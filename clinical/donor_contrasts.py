"""Case-minus-control contrasts of donor means: exhaustive two-sided label permutations, Welch 95% intervals,
and Benjamini-Hochberg q-values within families. Groups with fewer than two donors get only the difference.
"""
import itertools

import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

TIE_TOLERANCE = 1e-12


def exact_welch(case, control):
    """Case-minus-control difference of donor means with exact permutation P, Welch 95% interval and Welch P."""
    a, b = np.asarray(case, float), np.asarray(control, float)
    out = dict(n_case=len(a), n_control=len(b), effect=float(a.mean() - b.mean()) if len(a) and len(b) else np.nan,
               ci_low=np.nan, ci_high=np.nan, p_exact_two_sided=np.nan, p_welch_two_sided=np.nan,
               n_permutations=0, status='insufficient_independent_donors')
    if min(len(a), len(b)) < 2 or not (np.isfinite(a).all() and np.isfinite(b).all()):
        return out
    v, delta = np.r_[a, b], out['effect']
    null = []
    for ix in itertools.combinations(range(len(v)), len(a)):
        m = np.zeros(len(v), bool)
        m[list(ix)] = True
        null.append(v[m].mean() - v[~m].mean())
    v1, v0 = a.var(ddof=1) / len(a), b.var(ddof=1) / len(b)
    se = np.sqrt(v1 + v0)
    df = (v1 + v0) ** 2 / (v1 * v1 / (len(a) - 1) + v0 * v0 / (len(b) - 1)) if se else np.inf
    half_width = stats.t.ppf(.975, df) * se
    p_welch = 2 * stats.t.sf(abs(delta / se), df) if se else (1. if delta == 0 else 0.)
    out.update(ci_low=delta - half_width, ci_high=delta + half_width, se=se, welch_df=df,
               p_exact_two_sided=float(np.mean(np.abs(null) >= abs(delta) - TIE_TOLERANCE)),
               p_welch_two_sided=p_welch, n_permutations=len(null), status='estimated_exploratory_small_n')
    return out


def bh(df, families, p='p_exact_two_sided', out='q_BH'):
    """Add BH q-values computed within each family; rows without a P value stay NA."""
    df[out] = np.nan
    tested = df[p].notna()
    df.loc[tested, out] = df[tested].groupby(families)[p].transform(lambda s: multipletests(s, method='fdr_bh')[1])
    return df
