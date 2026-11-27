"""Independently reconstruct the graphs, purged folds, tuning and OOF metrics."""
from collections import Counter
from functools import cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

from Bio.Align import PairwiseAligner
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.stats import rankdata
from sklearn.feature_extraction import DictVectorizer
from sklearn.model_selection import GroupKFold
from sklearn.preprocessing import normalize

import prior_grouped as replay


def alignment_matrix(sequences):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score, aligner.mismatch_score = 2, -1
    aligner.open_gap_score, aligner.extend_gap_score = -5, -1

    @cache
    def measure(a, b):
        if a == b:
            return 1.
        path = aligner.align(a, b)[0].coordinates.T.tolist()
        matches, columns = 0, 0
        for (i, j), (u, v) in zip(path, path[1:]):
            columns += max(u - i, v - j)
            matches += sum(x == y for x, y in zip(a[i:u], b[j:v]))
        return matches / columns

    return np.array([[measure(*sorted((a, b))) for b in sequences] for a in sequences])


def graph(adjacency):
    return connected_components(csr_matrix(adjacency), directed=False)[1]


def equal(values):
    return np.equal.outer(values, values)


def rank(y, p):
    if np.ptp(y) == 0 or np.ptp(p) == 0:
        return None
    return float(np.corrcoef(rankdata(y), rankdata(p))[0, 1])


def refit(name, data, y, train, test, setting):
    if name == "embedding":
        return replay.embedding_fit(data, y, train, test, *setting)
    block = data[np.ix_(test, train)]
    target = y[train]
    if name == "trimer":
        normal = data[np.ix_(train, train)] + setting[0] * np.eye(len(train))
        pred = target.mean() + block @ np.linalg.solve(normal, target - target.mean())
    else:
        predictions = []
        for row in block:
            kth = sorted(row, reverse=True)[min(setting[0], len(row)) - 1]
            predictions.append(target[row >= kth].mean())
        pred = np.array(predictions)
    return pred.round(10)


def retained(candidate, held, adjacency, groups):
    rejected = {groups[i] for i in held}
    remaining = np.array([i for i in candidate if groups[i] not in rejected], dtype=int)
    assert not adjacency[np.ix_(remaining, held)].any()
    return remaining


