"""Five matched non-member controls per target from standardized WT node features (log1p mean CP10k, detection, log1p
variance, log1p outdegree, log1p outstrength): non-prior nodes within 1 SD on every feature, nearest five by Euclidean
distance (gene-name tie break), balanced if every five-control mean difference is <= 0.5 SD. Targets: 49 MP3 members plus
focal targets (primary), the five focal targets (seed17, down80_seed29).
Inputs: results/knockout/networks/<config>/<key>/WT_node_features.tsv, inputs/<key>/all_gene_features.tsv.gz
Outputs: results/knockout/matching/<config>/{target_eligibility,matched_controls,perturbation_jobs}.tsv
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd

from ko_common import (FOCAL_TARGETS, KO, N_CONTROLS, NETWORK_CONFIGS, SAMPLE_KEYS, all_targets, fixed_gene_sets,
                       prior_genes, read_tsv, write_tsv)

FEATURES = ['log1p_mean_CP10k', 'detection_fraction', 'log1p_variance_logCP10k', 'log1p_outdegree', 'log1p_outstrength']
CALIPER_SD = 1.0
CENTROID_SD = 0.5
TOL = 1e-12


def feature_z(nodes):
    X = np.column_stack([np.log1p(nodes.mean_CP10k), nodes.detection_fraction, np.log1p(nodes.variance_logCP10k),
                         np.log1p(nodes.outdegree), np.log1p(nodes.outstrength)])
    sd = X.std(0, ddof=1)
    assert np.all(sd > 0)
    return (X - X.mean(0)) / sd


def match_donor(cfg, key, targets, excluded, mp3):
    nodes = read_tsv(KO / 'networks' / cfg / key / 'WT_node_features.tsv').set_index('gene')
    features = read_tsv(KO / 'inputs' / key / 'all_gene_features.tsv.gz')
    measured = set(features.loc[~features.mapping_ambiguous, 'gene'])
    names = nodes.index.to_numpy()
    Z = feature_z(nodes)
    pool = np.flatnonzero(~nodes.index.isin(excluded) & (nodes.variance_logCP10k.to_numpy() > 0)
                          & (nodes.outdegree.to_numpy() > 0))
    quality, matches, perturb = [], [], set()
    for g in targets:
        in_net = g in nodes.index
        ok = in_net and nodes.at[g, 'outdegree'] > 0
        why = 'PASS' if ok else 'NO_WT_OUTGOING_EDGES' if in_net else 'NOT_IN_NETWORK' if g in measured else 'NOT_IN_ASSAY'
        n_caliper, picked, balanced, centroid = 0, [], False, np.nan
        if ok:
            perturb.add(g)
            i = np.flatnonzero(names == g)[0]
            delta = Z[pool] - Z[i]
            dist = np.sqrt(np.sum(delta ** 2, 1))
            passed = np.flatnonzero(np.max(abs(delta), axis=1) <= CALIPER_SD + TOL)
            n_caliper = len(passed)
            if n_caliper >= N_CONTROLS:
                picked = pool[sorted(passed, key=lambda k: (dist[k], names[pool[k]]))[:N_CONTROLS]]
                mean_delta = (Z[picked] - Z[i]).mean(0)
                centroid = float(np.max(abs(mean_delta)))
                balanced = centroid <= CENTROID_SD + TOL
                for rank, j in enumerate(picked, 1):
                    row = dict(config=cfg, sample_key=key, target=g, control=names[j], match_rank=rank,
                               distance=float(np.linalg.norm(Z[j] - Z[i])), match_pass=bool(balanced),
                               max_abs_centroid_SMD=centroid)
                    for k, feature in enumerate(FEATURES):
                        row['SMD_' + feature] = float(Z[j, k] - Z[i, k])
                        row['mean5_SMD_' + feature] = float(mean_delta[k])
                    matches.append(row)
                if balanced:
                    perturb.update(names[picked])
        reason = ('PASS' if balanced else why if not ok else
                  'FEWER_THAN_FIVE_CALIPER_CONTROLS' if n_caliper < N_CONTROLS else 'FIVE_CONTROL_CENTROID_IMBALANCE')
        quality.append(dict(config=cfg, sample_key=key, target=g, raw_present=g in measured, in_network=in_net,
                            KO_eligible=bool(ok), KO_reason=why, n_caliper_controls=n_caliper,
                            n_selected_controls=len(picked), control_matching_eligible=bool(balanced),
                            matching_reason=reason, max_abs_centroid_SMD=centroid))
    jobs = [dict(config=cfg, sample_key=key, perturbed_gene=g,
                 role='systematic_target' if g in mp3 else 'candidate_reference' if g in FOCAL_TARGETS else 'matched_control')
            for g in sorted(perturb)]
    return quality, matches, jobs


sets = fixed_gene_sets()
excluded = prior_genes(sets)
for cfg in NETWORK_CONFIGS:
    targets = all_targets(sets) if cfg == 'primary' else FOCAL_TARGETS
    parts = [match_donor(cfg, key, targets, excluded, set(sets['MP_3'])) for key in SAMPLE_KEYS]
    folder = KO / 'matching' / cfg
    folder.mkdir(parents=True, exist_ok=True)
    for i, name in enumerate(['target_eligibility.tsv', 'matched_controls.tsv', 'perturbation_jobs.tsv']):
        write_tsv(pd.DataFrame(sum((p[i] for p in parts), [])), folder / name)
    print(cfg, sum(len(p[2]) for p in parts), 'knockouts,',
          sum(q['control_matching_eligible'] for p in parts for q in p[0]), 'matched targets')
