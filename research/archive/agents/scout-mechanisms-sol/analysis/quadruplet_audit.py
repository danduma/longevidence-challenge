"""Exact audits of synthetic code maps; no cell or decoder is implemented.

Compare padding, a known parity-code idea, and a rotation-invariant dense map.
The latter is a mathematical control for joint position/information audits,
not a claimed novel code. All weighting assumptions are explicit.
"""
import collections
import csv
import itertools
import json
import math
import random
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
BASES='ACGU'
AA='ACDEFGHIKLMNPQRSTVWY'
WORDS=[''.join(x) for x in itertools.product(BASES,repeat=4)]


def rotations(w):return {w[i:]+w[:i] for i in range(4)}


def audit(code):
    sense=[w for w in WORDS if code[w] not in ['*','!']]
    changes=[]; erasures=[]; entropy=[]; ambiguous=[]
    for pos in range(4):
        miss=erased=0
        groups=collections.defaultdict(collections.Counter)
        for w in sense:
            groups[w[:pos]+w[pos+1:]][code[w]]+=1
            for base in BASES:
                if base==w[pos]:continue
                v=w[:pos]+base+w[pos+1:]
                miss+=code[v]!=code[w]
                erased+=code[v] in ['*','!']
        changes.append(miss/(3*len(sense)))
        erasures.append(erased/(3*len(sense)))
        h=0.0
        for counts in groups.values():
            n=sum(counts.values())
            h+=n/len(sense)*sum(-c/n*math.log2(c/n) for c in counts.values())
        entropy.append(h)
        ambiguous.append(sum(len(c)>1 for c in groups.values()))
    # Uniform, independent valid codon pairs. Check EVERY internal frame offset.
    accepted=[]
    for shift in [1,2,3]:
        count=sum(code[a[shift:]+b[:shift]] not in ['*','!']
                  for a in sense for b in sense)
        accepted.append(count/len(sense)**2)
    return {'sense_codons':len(sense),'amino_acids':len({code[w] for w in sense}),
            'position_change_including_erasure':changes,
            'position_erasure_probability':erasures,
            'conditional_amino_acid_entropy_given_other_3_bases_bits':entropy,
            'ambiguous_three_base_keys_by_deleted_position':ambiguous,
            'shifted_window_sense_acceptance_1_2_3':accepted,
            'weighting':'uniform over sense codons; uniform substitutions; independent uniform codon pairs'}


def main():
    # Deliberately padded artificial 20-AA code, not the standard genetic code.
    padded={w:AA[(BASES.index(w[0])*16+BASES.index(w[1])*4+BASES.index(w[2]))%20]
            for w in WORDS}
    parity=[w for w in WORDS if sum(BASES.index(c) for c in w)%4==0]
    sparse={w:'!' for w in WORDS}
    for w,a in zip(parity[:20],AA):sparse[w]=a
    # Aperiodic length-4 necklaces: 60 orbits of size 4; 16 periodic words STOP.
    reps=sorted({min(rotations(w)) for w in WORDS if len(rotations(w))==4})
    assert len(reps)==60
    oi={w:i for i,r in enumerate(reps) for w in rotations(r)}
    weights=collections.defaultdict(int)
    for w in oi:
        for pos in range(4):
            for b in BASES:
                v=w[:pos]+b+w[pos+1:]
                if b==w[pos] or v not in oi:continue
                a,c=sorted((oi[w],oi[v]))
                if a!=c:weights[a,c]+=1
    # Directed edges are counted twice above; either count gives same optimum.
    adj=[{} for _ in reps]
    for (a,b),v in weights.items():adj[a][b]=v;adj[b][a]=v
    rng=random.Random(20261126)
    labels=list(AA*3);rng.shuffle(labels)
    initial=list(labels)
    def score(ls):return sum(v for (a,b),v in weights.items() if ls[a]==ls[b])
    current=score(labels);best=current;bestlabels=list(labels)
    for step in range(120000):
        a,b=rng.sample(range(60),2)
        if labels[a]==labels[b]:continue
        la,lb=labels[a],labels[b]
        delta=sum(v*((lb==labels[k])-(la==labels[k])) for k,v in adj[a].items() if k!=b)
        delta+=sum(v*((la==labels[k])-(lb==labels[k])) for k,v in adj[b].items() if k!=a)
        temp=max(0.05,8*(1-step/120000)**2)
        if delta>=0 or rng.random()<math.exp(delta/temp):
            labels[a],labels[b]=lb,la;current+=delta
            if current>best:best=current;bestlabels=list(labels)
    assert score(bestlabels)==best
    def expand(ls):return {w:ls[oi[w]] if w in oi else '*' for w in WORDS}
    dense=expand(bestlabels)
    assert all(dense[w]==dense[w[1:]+w[:1]] for w in WORDS)
    assert collections.Counter(dense[w] for w in WORDS if dense[w]!='*')=={a:12 for a in AA}
    controls={'padded_artificial_20AA':audit(padded),
              'sparse_parity_20AA':audit(sparse),
              'rotation_dense_initial':audit(expand(initial)),
              'rotation_dense_optimized':audit(dense)}
    for k in ['rotation_dense_initial','rotation_dense_optimized']:
        for metric in ['position_change_including_erasure',
                       'conditional_amino_acid_entropy_given_other_3_bases_bits']:
            assert max(controls[k][metric])-min(controls[k][metric])<1e-12
    assert controls['sparse_parity_20AA']['position_change_including_erasure']==[1]*4
    assert controls['sparse_parity_20AA']['conditional_amino_acid_entropy_given_other_3_bases_bits']==[0]*4
    out={'seed':20261126,'annealing_steps':120000,'codes':controls,
         'single_error_correction_hamming_ball_size':13,
         '20_labels_times_ball_size':260,'available_words':256,
         'initial_synonymous_edge_objective':score(initial),'optimized_objective':best,
         'warning':'Heuristic code construction, not globally optimal; no translation, triplet exclusion, viability or novelty claim.'}
    (ROOT/'outputs/quadruplet_audit.json').write_text(json.dumps(out,indent=2)+'\n')
    with (ROOT/'outputs/quadruplet_dense_map.csv').open('w') as f:
        writer=csv.writer(f);writer.writerow(['codon','label'])
        writer.writerows(sorted(dense.items()))
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
