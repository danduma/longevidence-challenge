#!/usr/bin/env python3
"""Independent replication and prespecified negative-control transfer analysis."""
from pathlib import Path
import hashlib
import json
import zipfile
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.sparse.linalg import spsolve

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "scout-data-sol" / "data"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
WT = "ENLYFQS"
AA = "ACDEFGHIKLMNPQRSTVWY"
CANON = f"[{AA}]{{7}}"
NEG = ["(1-13)-mCherry", "I C151A 30*31*", "no TEVp"]
PAIR = ["TEVp", "TEVs"]


def read(suffix):
    return pd.read_csv(DATA / f"protease_{suffix}.tsv", sep="\t", dtype={"TEVp": str})


def hamming(motif):
    return sum(a != b for a, b in zip(motif, WT))


def transfer(neg):
    records = []
    for panel, depth, max_mut in [("primary", 100, 7), ("singles_and_wt", 100, 1),
                                  ("all_depths", 0, 7)]:
        use = neg[neg.TEVs.str.fullmatch(CANON) & (neg["min"] >= depth)
                  & (neg.TEVs.map(hamming) <= max_mut)]
        table = use.pivot_table(index="TEVs", columns="TEVp", values="AUC").dropna()
        for name in NEG:
            held = table[name]
            other = table.drop(columns=name).mean(axis=1)
            high = held >= held.quantile(.9)
            called = other >= other.quantile(.9)
            records.append(dict(panel=panel, held=name, motifs=len(table),
                                rho=held.corr(other, method="spearman"), selected=called.sum(),
                                held_high=high.sum(), overlap=(called & high).sum(),
                                precision=high[called].mean(), base=high.mean()))
    pd.DataFrame(records).to_csv(OUT / "control_transfer.csv", index=False)


def cross_screen(neg):
    raw = (ROOT / "sources/screen3.tsv").read_bytes()
    blob = b"blob " + str(len(raw)).encode() + b"\x00" + raw
    assert hashlib.sha1(blob).hexdigest() == "a636c1a8dc9134fcaf098a666d988b06361e76df"
    old = pd.read_csv(ROOT / "sources/screen3.tsv", sep="\t")
    lib = old["library.ID"]
    use = (lib.isin(["libLH015","libLH016","libLH019"]) &
           ~(old.TEVs.eq(WT) & lib.isin(["libLH016","libLH019"])) &
           old.TEVp.isin(NEG) & old.TEVs.str.fullmatch(CANON) & (old["min"] >= 100))
    ref = neg[neg.TEVs.str.fullmatch(CANON) & (neg["min"] >= 100)].pivot(
        index="TEVs",columns="TEVp",values="AUC").dropna().mean(axis=1)
    rows = []
    for rep,group in old[use].groupby("rep"):
        mat = group.pivot_table(index="TEVs",columns="TEVp",values="AUC4").reindex(columns=NEG)
        obs = mat.dropna().mean(axis=1)
        panel = pd.concat([ref.rename("ref"),obs.rename("obs")],axis=1).dropna()
        selected = panel.ref >= panel.ref.quantile(.9)
        high = panel.obs >= panel.obs.quantile(.9)
        rows.append(dict(rep=rep,n=len(panel),rho=panel.ref.corr(panel.obs,method="spearman"),
                         selected=selected.sum(),high=high.sum(),overlap=(high & selected).sum()))
    pd.DataFrame(rows).to_csv(OUT/"cross_screen.csv",index=False)


def matched(controls, neg, info):
    pos = ["P6", "P5", "P4", "P3", "P2", "P1", "P1'"]
    v = read("vitro_relative").merge(read("vitro_absolute"), on=["TEVp","pos","TEVs_AA"],
                                     validate="one_to_one")
    v["TEVs"] = [WT[:i] + aa + WT[i+1:] for i, aa in zip(v.pos.map(pos.index), v.TEVs_AA)]
    rates = v.groupby(PAIR)[["vinirel","viniabs"]].mean()
    c = controls.set_index(PAIR)
    assert c.index.is_unique
    joined = rates.join(c[["AUC","min"]], how="left")
    info["validation"] = dict(original=len(v), unique=len(rates),
                              missing=[list(k) for k in joined[joined.AUC.isna()].index],
                              relative_missing=int(v.vinirel.isna().sum()),
                              wt_relative_range=v[v.TEVs==WT].vinirel.agg(["min","max"]).tolist(),
                              negative_relative=int((v.vinirel<0).sum()))
    joined = joined.dropna(subset="AUC").reset_index()
    joined["bg"] = joined.TEVs.map(neg.groupby("TEVs").AUC.mean())
    joined["subtracted"] = joined.AUC - joined.bg
    rows = []
    for name,g in joined.groupby("TEVp"):
        rows.append(dict(enzyme=name,n=len(g),min_depth=g["min"].min(),under100=(g["min"]<100).sum(),
                         raw=g.AUC.corr(g.vinirel, method="spearman"),
                         subtracted=g.subtracted.corr(g.vinirel,method="spearman"),
                         background=g.bg.corr(g.vinirel,method="spearman"),
                         absolute=g.AUC.corr(g.viniabs,method="spearman")))
    pd.DataFrame(rows).to_csv(OUT/"validation.csv",index=False)


def ap(labels, scores):
    # Exact grouped precision at each distinct score, weighted by new positives.
    bins = pd.DataFrame(dict(label=labels, score=scores)).groupby("score").label.agg(["sum","size"])
    descending = bins.iloc[::-1]
    cumulative = descending.cumsum()
    return (descending["sum"] * cumulative["sum"] / cumulative["size"]).sum() / labels.sum()


def count_audit(frame):
    flipped = frame[[f"F{i}" for i in range(1,6)]].to_numpy()
    totals = flipped + frame[[f"N{i}" for i in range(1,6)]].to_numpy()
    auc = np.trapezoid(flipped / totals, axis=1) / 4
    return dict(depth_disagreements=int((totals.min(axis=1)!=frame["min"]).sum()),
                max_auc_abs_error=float(np.max(np.abs(auc-frame.AUC))))


def design(d):
    # Encode only non-reference residues. Column order has no effect on ridge.
    def onehot(values):
        codes, unique = pd.factorize(values)
        return sparse.csr_matrix((np.ones(len(d)), (np.arange(len(d)), codes)),
                                 shape=(len(d),len(unique)))
    parts = [sparse.csr_matrix(np.ones((len(d),1))), onehot(d.TEVp)]
    for j, ref in enumerate(WT):
        letters = d.TEVs.str[j]
        parts.append(onehot(letters).multiply((letters!=ref).to_numpy()[:,None]))
    X = sparse.hstack(parts,format="csr")
    X.eliminate_zeros()
    return X


def replication(info, cutoff, matched_cutoff):
    with zipfile.ZipFile(DATA/"protease_DMS.zip") as ar, ar.open("Data_DMS.txt") as f:
        d = pd.read_csv(f,sep="\t")
    info["DMS_counts"] = count_audit(d)
    assert not d.duplicated(PAIR).any()
    info["filter_waterfall"] = {"source":len(d)}
    for label, mask in [("depth100",d["min"]>=100), ("enzyme_mut_count",d.n_mut_TEVp<=1),
                        ("canonical_motif",d.TEVs.str.fullmatch(CANON)),
                        ("wt_or_single",d.TEVp.str.fullmatch(r"wt|[A-Z]\d+[A-Z]"))]:
        d = d.loc[mask.loc[d.index]]
        info["filter_waterfall"][label] = len(d)
    assert d.TEVs.map(hamming).eq(d.n_mut_TEVs).all()
    info["enzyme_labels"] = d.TEVp.nunique()
    info["motifs"] = d.TEVs.nunique()
    info["motifs_by_distance"] = d.drop_duplicates("TEVs").n_mut_TEVs.value_counts().to_dict()
    X = design(d)
    y = d.AUC.to_numpy()
    rng = np.random.default_rng(20261126)
    sensitivity = []
    records = []
    for name, test in [("random",rng.uniform(size=len(d))<.2),("blocked",d.n_mut_TEVs.to_numpy()>1)]:
        train = d.loc[~test]
        hold = d.loc[test]
        T = X[~test]
        mu = train.AUC.mean()
        # Equivalent penalized least squares via direct normal equations.
        gram = (T.T @ T).tocsc() + sparse.eye(T.shape[1],format="csc")
        coef = spsolve(gram, T.T @ (y[~test]-mu))
        pred = X[test] @ coef + mu
        support = np.asarray(T.getnnz(axis=0))>0
        info[name] = dict(train=len(train),test=len(hold),
                          unseen_enzymes=len(set(hold.TEVp)-set(train.TEVp)),
                          unseen_motifs=len(set(hold.TEVs)-set(train.TEVs)),
                          unsupported_rows=int((X[test][:,~support].getnnz(axis=1)>0).sum()))
        matched_labels = y[test] > matched_cutoff
        sensitivity.append(dict(split=name, cutoff=matched_cutoff, positives=matched_labels.sum(),
                                prevalence=matched_labels.mean(), ridge_ap=ap(matched_labels,pred)))
        labels = y[test]>cutoff
        estimates = {"mean":np.full(test.sum(),mu), "ridge":pred}
        for key in PAIR:
            estimates[key] = hold[key].map(train.groupby(key).AUC.mean()).fillna(mu).to_numpy()
        for model,scores in estimates.items():
            records.append(dict(split=name,model=model,rmse=np.sqrt(np.mean((scores-y[test])**2)),
                                rho=stats.spearmanr(scores,y[test]).statistic if np.ptp(scores) else None,
                                positives=labels.sum(),prevalence=labels.mean(),ap=ap(labels,scores)))
        if name=="blocked":
            groups = pd.DataFrame(dict(motif=hold.TEVs,n=1, ridge=(y[test]-pred)**2,
                                       mean=(y[test]-mu)**2)).groupby("motif").sum()
            ix = rng.integers(len(groups), size=(3000,len(groups)))
            n, sq, null = groups.to_numpy()[ix].sum(axis=1).T
            info[name].update(ci95=np.quantile(np.sqrt(sq/n)-np.sqrt(null/n), [.025,.975]).tolist(),
                              n_range=groups.n.agg(["min","max"]).tolist(),
                              groups_ridge_worse=int((groups.ridge>groups["mean"]).sum()),
                              held_enzymes=hold.TEVp.nunique(),
                              min_prediction=float(pred.min()), negative_predictions=int((pred<0).sum()))
    pd.DataFrame(sensitivity).to_csv(OUT/"threshold.csv",index=False)
    pd.DataFrame(records).to_csv(OUT/"replication.csv",index=False)


def main():
    info = {}
    for row in pd.read_csv(DATA/"file_inventory.csv").itertuples(index=False,name=None):
        path = DATA/Path(row[0]).name
        if path.name.startswith("protease_"):
            assert hashlib.sha256(path.read_bytes()).hexdigest()==row[2]
    c = read("controls")
    neg = c[c.TEVp.isin(NEG)]
    info["control_counts"] = count_audit(c)
    info["controls"] = dict(total=len(c),negative=len(neg),min_depth=int(neg["min"].min()),
                            below100=int((neg["min"]<100).sum()),
                            noncanonical=int((~neg.TEVs.str.fullmatch(CANON)).sum()),
                            cutoff=float(neg.AUC.quantile(.99)),
                            canonical_cutoff=float(neg[neg.TEVs.str.fullmatch(CANON)].AUC.quantile(.99)),
                            exceed_noncanonical=int((~neg.TEVs.str.fullmatch(CANON) &
                                                     (neg.AUC>neg.AUC.quantile(.99))).sum()))
    transfer(neg)
    cross_screen(neg)
    matched(c,neg,info)
    eligible = c.TEVs.str.fullmatch(CANON) & (c["min"]>=100)
    four = c[eligible & c.TEVp.isin(NEG+["MBP"])].pivot(index="TEVs", columns="TEVp", values="AUC").dropna()
    reference = four[NEG].mean(axis=1)
    selected = reference >= reference.quantile(.9)
    high = four.MBP >= four.MBP.quantile(.9)
    info["supplement_MBP"] = dict(motifs=len(four), rho=four.MBP.corr(reference,method="spearman"),
                                  selected=int(selected.sum()), high=int(high.sum()),
                                  overlap=int((selected & high).sum()))
    elig_neg = c[eligible & c.TEVp.isin(NEG)]
    info["controls"]["matched_eligible"] = len(elig_neg)
    replication(info,info["controls"]["cutoff"],elig_neg.AUC.quantile(.99))
    (OUT/"audit.json").write_text(json.dumps(info,indent=2,allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
