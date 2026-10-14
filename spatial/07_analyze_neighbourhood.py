"""Somatic neighbour weights around HC-SPG beads in fixed distance annuli, related to MP3 (Fig 6g reports the Sertoli
OOF estimator). Focal beads: HC-SPG beads >= 8 median-NN units inside the convex hull. Three estimators per annulus:
high-minus-low MP3 tertile (unbalanced), the same on propensity-matched pairs, and the slope of neighbour weight on the
standardized out-of-fold MP3 residual (ridge cross-fitted over four quadrants, 5-NN-unit buffer). Intervals: 399
resamples of the original 6x6 spatial blocks; one random stream per array is shared by the estimators in table order.
Input: results/spatial/neighbourhood/<array>_bead_annotation.tsv.gz; output: neighbourhood/neighbour_effects.tsv.
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.spatial import ConvexHull, cKDTree
from scipy.spatial.distance import cdist
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, SplineTransformer, StandardScaler

from spatial_common import ARRAYS, HC_CLASSES, NEIGHBOURHOOD, TYPES

SEED = 20261004
HULL_MARGIN = 8          # median-NN units
RADIUS = 8               # local-density radius and outer annulus edge, median-NN units
BUFFER = 5               # training beads closer than this to the held-out quadrant are dropped
ANNULI = [(0, 2), (2, 4), (4, 6), (6, 8)]
NEIGHBOUR_TYPES = ['Sertoli', 'Peritubular_smooth_muscle', 'Leydig']
N_BLOCKS = 6             # 6 x 6 spatial blocks
N_BOOT = 399
MIN_NEIGHBOURS = 3
COVARIATES = ['logUMI', 'logdensity8', 'SPG', 'Sertoli', 'Peritubular_smooth_muscle', 'Leydig']
CALIPER = 0.2            # SD of the logit propensity
MAX_ABS_SMD = 0.1
MIN_PAIRS = 30
MIN_PAIR_BLOCKS = 6
MIN_UNBALANCED = 60
MIN_OOF_ROWS = 100
MIN_OOF_BLOCKS = 6
METHODS = ['original_unbalanced', 'matched_primary', 'OOF_continuous_secondary']


def bead_table(meta):
    xy = meta[['x_native', 'y_native']].to_numpy()
    tree = cKDTree(xy)
    mnn = np.median(tree.query(xy, k=2)[0][:, 1])
    hull = ConvexHull(xy)
    margin = np.min(-(xy @ hull.equations[:, :2].T + hull.equations[:, 2]), axis=1)
    meta['logdensity8'] = np.log1p(tree.query_ball_point(xy, RADIUS * mnn, return_length=True) - 1)
    bb = np.minimum(((xy - xy.min(0)) / (xy.max(0) - xy.min(0)) * N_BLOCKS).astype(int), N_BLOCKS - 1)
    meta['block'] = bb[:, 0] + N_BLOCKS * bb[:, 1]
    centre = np.median(xy, 0)
    meta['fold'] = (xy[:, 0] >= centre[0]).astype(int) + 2 * (xy[:, 1] >= centre[1]).astype(int)
    ix = np.flatnonzero(meta.HC_SPG.to_numpy() & (margin >= HULL_MARGIN * mnn))
    f = meta.iloc[ix].copy()
    cuts = f.MP_3.quantile([1 / 3, 2 / 3]).to_numpy()
    f['group'] = 'middle'
    f.loc[f.MP_3 <= cuts[0], 'group'] = 'low'
    f.loc[f.MP_3 >= cuts[1], 'group'] = 'high'
    f['global_index'] = ix
    return xy, tree, mnn, f


def oof_residuals(f, mnn):
    """MP3 minus its out-of-fold ridge prediction from depth, density, reported weights and spot class."""
    features = ['logUMI', 'logdensity8'] + TYPES + ['spot_class']
    f['prediction'] = np.nan
    for fold in range(4):
        test = f.fold == fold
        ti, tr = np.flatnonzero(test.to_numpy()), np.flatnonzero((~test).to_numpy())
        dist = cKDTree(f.iloc[ti][['x_native', 'y_native']]).query(f.iloc[tr][['x_native', 'y_native']], k=1)[0]
        tr = tr[dist > BUFFER * mnn]
        if len(tr) < 100 or len(ti) < 20:
            continue
        pre = ColumnTransformer([('depth', SplineTransformer(n_knots=5, degree=3, include_bias=False), ['logUMI']),
                                 ('continuous', StandardScaler(), ['logdensity8'] + TYPES),
                                 ('class', OneHotEncoder(handle_unknown='ignore'), ['spot_class'])])
        fit = make_pipeline(pre, StandardScaler(with_mean=False), Ridge(alpha=1, solver='lsqr'))
        fit.fit(f.iloc[tr][features], f.MP_3.iloc[tr])
        f.loc[f.index[ti], 'prediction'] = fit.predict(f.iloc[ti][features])
    f['OOF_residual'] = f.MP_3 - f.prediction


def smd(high, low, sd):
    """Standardized mean difference on the pooled pre-matching SD (0 for equal means with zero SD)."""
    if sd > 0:
        return (high - low) / sd
    return 0. if sd == 0 and high == low else np.nan


def match_tertiles(f):
    """High/low MP3 pairs: exact spot class, logit-propensity caliper, minimum-cost Euclidean matching on covariates."""
    hl = f[f.group.isin(['low', 'high'])].copy()
    X = hl[COVARIATES].to_numpy()
    sd = X.std(0, ddof=1)
    Z = (X - X.mean(0)) / np.where(sd > 0, sd, 1)
    t = (hl.group == 'high').astype(int).to_numpy()
    eta = LogisticRegression(C=1, max_iter=2000, random_state=SEED).fit(Z, t).decision_function(Z)
    caliper = CALIPER * eta.std(ddof=1)
    smd_sd = np.sqrt((X[t == 1].var(0, ddof=1) + X[t == 0].var(0, ddof=1)) / 2)   # pooled SD before matching
    pairs = []
    for sc in sorted(hl.spot_class.unique()):
        hi = np.flatnonzero((hl.group == 'high') & (hl.spot_class == sc))
        lo = np.flatnonzero((hl.group == 'low') & (hl.spot_class == sc))
        if not len(hi) or not len(lo):
            continue
        mn, mx = max(eta[hi].min(), eta[lo].min()), min(eta[hi].max(), eta[lo].max())
        h, l = hi[(eta[hi] >= mn) & (eta[hi] <= mx)], lo[(eta[lo] >= mn) & (eta[lo] <= mx)]
        if not len(h) or not len(l):
            continue
        costs = cdist(Z[h], Z[l])
        costs[np.abs(eta[h, None] - eta[l][None, :]) > caliper] = 1e6
        i, j = linear_sum_assignment(costs)
        pairs += [(hl.iloc[hh].barcode, hl.iloc[ll].barcode) for hh, ll, cost in zip(h[i], l[j], costs[i, j]) if cost < 1e6]
    pairs = pd.DataFrame(pairs, columns=['high_barcode', 'low_barcode'])
    matched = hl[hl.barcode.isin(set(pairs.high_barcode) | set(pairs.low_barcode))]
    smds = [smd(matched.loc[matched.group == 'high', c].mean(), matched.loc[matched.group == 'low', c].mean(), smd_sd[j])
            if (matched.group == 'high').any() and (matched.group == 'low').any() else np.nan for j, c in enumerate(COVARIATES)]
    blocks = [matched.loc[matched.group == g, 'block'].nunique() for g in ['high', 'low']]
    passed = bool(len(pairs) >= MIN_PAIRS and min(blocks) >= MIN_PAIR_BLOCKS and np.isfinite(smds).all()
                  and np.max(np.abs(smds)) <= MAX_ABS_SMD)
    return pairs, smd_sd, passed


def annulus_rows(f, meta, xy, tree, mnn):
    hc = meta.spot_class.isin(HC_CLASSES).to_numpy()
    rows = []
    for row in f.itertuples():
        i = int(row.global_index)
        nb = np.array(tree.query_ball_point(xy[i], RADIUS * mnn))
        dd = np.linalg.norm(xy[nb] - xy[i], axis=1) / mnn
        for lo, hi in ANNULI:
            js = nb[(dd > lo) & (dd <= hi) & hc[nb]]
            if len(js) < MIN_NEIGHBOURS:
                continue
            for ct in NEIGHBOUR_TYPES:
                rows.append(dict(barcode=row.barcode, group=row.group, block=int(row.block), mid=(lo + hi) / 2, type=ct,
                                 weight=float(meta.iloc[js][ct].mean()), OOF_residual=row.OOF_residual))
    return pd.DataFrame(rows)


def block_bootstrap(ss, method, scale, rng):
    blocks = ss.block.unique()
    draw = rng.integers(0, len(blocks), size=(N_BOOT, len(blocks)))
    if method.startswith('OOF'):
        zz = ss.assign(xx=ss.OOF_residual / scale)
        zz = zz.assign(xy=zz.xx * zz.weight, x2=zz.xx ** 2)
        bb = zz.groupby('block').agg(n=('weight', 'size'), sx=('xx', 'sum'), sy=('weight', 'sum'), sxy=('xy', 'sum'),
                                     sx2=('x2', 'sum')).reindex(blocks)
        n, sx, sy, sxy, sx2 = bb.to_numpy()[draw].sum(1).T
        den = sx2 - sx * sx / n
        return np.divide(sxy - sx * sy / n, den, out=np.full(N_BOOT, np.nan), where=den > 0)
    bb = ss.groupby(['block', 'group']).weight.agg(['sum', 'count']).unstack('group').reindex(blocks).fillna(0)
    n1, n0 = bb['count']['high'].to_numpy()[draw].sum(1), bb['count']['low'].to_numpy()[draw].sum(1)
    s1, s0 = bb['sum']['high'].to_numpy()[draw].sum(1), bb['sum']['low'].to_numpy()[draw].sum(1)
    return (np.divide(s1, n1, out=np.full(N_BOOT, np.nan), where=n1 > 0)
            - np.divide(s0, n0, out=np.full(N_BOOT, np.nan), where=n0 > 0))


def point_estimate(ss, method, scale):
    if method.startswith('OOF'):
        xx, yy = ss.OOF_residual.to_numpy() / scale, ss.weight.to_numpy()
        return float(np.cov(xx, yy, ddof=0)[0, 1] / np.var(xx)) if np.var(xx) > 0 else np.nan
    h, l = ss[ss.group == 'high'].weight, ss[ss.group == 'low'].weight
    return h.mean() - l.mean() if len(h) and len(l) else np.nan


def neighbour_effects(array):
    meta = pd.read_csv(NEIGHBOURHOOD / f'{array}_bead_annotation.tsv.gz', sep='\t')
    xy, tree, mnn, f = bead_table(meta)
    oof_residuals(f, mnn)
    pairs, smd_sd, passed = match_tertiles(f)
    nr = annulus_rows(f, meta, xy, tree, mnn)
    rng = np.random.default_rng(SEED)
    estimates = []
    for (ct, mid), sub in nr.groupby(['type', 'mid']):
        for method in METHODS:
            extra = {}
            if method == 'original_unbalanced':
                ss = sub[sub.group.isin(['high', 'low'])]
                estimable = len(ss) >= MIN_UNBALANCED
            elif method == 'matched_primary':
                present = set(sub.barcode)
                rp = pairs[pairs.high_barcode.isin(present) & pairs.low_barcode.isin(present)]
                keep = set(rp.high_barcode) | set(rp.low_barcode)
                ss = sub[sub.barcode.isin(keep)]
                cov = f[f.barcode.isin(keep)]
                hs, ls = cov[cov.group == 'high'], cov[cov.group == 'low']
                smds = [smd(hs[c].mean(), ls[c].mean(), smd_sd[j]) if len(rp) else np.nan for j, c in enumerate(COVARIATES)]
                estimable = bool(passed and len(rp) >= MIN_PAIRS and hs.block.nunique() >= MIN_PAIR_BLOCKS
                                 and ls.block.nunique() >= MIN_PAIR_BLOCKS and np.isfinite(smds).all()
                                 and np.max(np.abs(smds)) <= MAX_ABS_SMD)
                extra = dict(n_complete_pairs=len(rp), ring_gate_pass=estimable,
                             ring_max_abs_SMD=float(np.nanmax(np.abs(smds))) if np.isfinite(smds).any() else np.nan)
            else:
                ss = sub.dropna(subset=['OOF_residual'])
                estimable = len(ss) >= MIN_OOF_ROWS and ss.block.nunique() >= MIN_OOF_BLOCKS and ss.OOF_residual.std() > 0
            out = dict(array=array, type=ct, mid=mid, method=method, n=len(ss), n_blocks=ss.block.nunique(),
                       estimable=bool(estimable), matching_pass=passed, effect=np.nan, low=np.nan, high=np.nan, **extra)
            if estimable:
                scale = ss.OOF_residual.std(ddof=0) if method.startswith('OOF') else 1
                effect = point_estimate(ss, method, scale)
                boots = block_bootstrap(ss, method, scale, rng)
                out.update(effect=effect, low=float(np.nanquantile(boots, .025)), high=float(np.nanquantile(boots, .975)),
                           n_valid_bootstraps=int(np.isfinite(boots).sum()), scale_residual_SD=float(scale))
            estimates.append(out)
    print(array, len(f), 'focal HC-SPG beads,', len(pairs), 'matched pairs, balance', 'passed' if passed else 'failed')
    return pd.DataFrame(estimates)


pd.concat([neighbour_effects(a) for a in ARRAYS], ignore_index=True).to_csv(
    NEIGHBOURHOOD / 'neighbour_effects.tsv', sep='\t', index=False)
