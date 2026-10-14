"""Shared controls for the same-node Full-SPG / U-only networks of donors A1 and A2 (targets FGFR3 and COMMD3).
Features standardized within each network: log1p mean CP10k, detection, log1p outdegree, log1p outstrength.
The five controls are the nearest non-prior nodes by mean distance across both networks (gene-name tie break) and are
used in both networks; per-feature caliper and centroid balance are recorded as diagnostics only.
Inputs: results/knockout/networks/<paired config>/<key>/WT_node_features.tsv
Outputs: results/knockout/matching/<paired config>/{target_eligibility,matched_controls,perturbation_jobs}.tsv
"""
import sys; from pathlib import Path; sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
import pandas as pd

from ko_common import KO, N_CONTROLS, PAIRED_CONFIGS, PAIRED_KEYS, PAIRED_TARGETS, fixed_gene_sets, prior_genes, read_tsv

CALIPER_SD = 1.0
CENTROID_SD = 0.5


def feature_z(nodes):
    X = np.column_stack([np.log1p(nodes.mean_CP10k), nodes.detection_fraction, np.log1p(nodes.outdegree),
                         np.log1p(nodes.outstrength)])
    sd = X.std(0, ddof=1)
    assert np.all(sd > 0)
    return pd.DataFrame((X - X.mean(0)) / sd, index=nodes.index)


excluded = prior_genes(fixed_gene_sets())
matches, quality, jobs = [], [], []
for key in PAIRED_KEYS:
    nodes = {c: read_tsv(KO / 'networks' / c / key / 'WT_node_features.tsv').set_index('gene') for c in PAIRED_CONFIGS}
    zs = {c: feature_z(nodes[c]) for c in PAIRED_CONFIGS}
    first = nodes[PAIRED_CONFIGS[0]]
    assert first.index.equals(nodes[PAIRED_CONFIGS[1]].index)
    pool = sorted(g for g in first.index if g not in excluded and all(nodes[c].at[g, 'outdegree'] > 0 for c in PAIRED_CONFIGS))
    perturb = {c: set() for c in PAIRED_CONFIGS}
    for target in PAIRED_TARGETS:
        ok = all(target in nodes[c].index and nodes[c].at[target, 'outdegree'] > 0 for c in PAIRED_CONFIGS)
        balanced = ok and len(pool) >= N_CONTROLS
        controls = []
        if balanced:
            distances = {c: np.sqrt(((zs[c].loc[pool] - zs[c].loc[target]) ** 2).sum(1)) for c in PAIRED_CONFIGS}
            joint = sum(distances.values()) / len(PAIRED_CONFIGS)
            controls = sorted(pool, key=lambda g: (joint.loc[g], g))[:N_CONTROLS]
        for c in PAIRED_CONFIGS:
            in_net = target in nodes[c].index
            eligible = in_net and nodes[c].at[target, 'outdegree'] > 0
            if eligible:
                perturb[c].add(target)
            if balanced:
                perturb[c].update(controls)
            for rank, g in enumerate(controls, 1):
                delta = zs[c].loc[g] - zs[c].loc[target]
                centroid = (zs[c].loc[controls] - zs[c].loc[target]).mean(0)
                row = dict(config=c, sample_key=key, target=target, control=g, match_rank=rank, distance=distances[c].at[g],
                           joint_distance=joint.at[g], match_pass=True,
                           individual_caliper_pass=bool(abs(delta).max() <= CALIPER_SD),
                           centroid_diagnostic_pass=bool(abs(centroid).max() <= CENTROID_SD),
                           max_abs_centroid_SMD=float(abs(centroid).max()))
                for f, v in zip(['mean', 'detection', 'degree', 'strength'], delta):
                    row[f'SMD_{f}'] = v
                matches.append(row)
            quality.append(dict(config=c, sample_key=key, target=target, raw_present=True, in_network=in_net,
                                KO_eligible=eligible, KO_reason='PASS' if eligible else 'NO_TARGET_NETWORK_EDGES',
                                n_caliper_controls=len(pool), n_selected_controls=len(controls),
                                control_matching_eligible=balanced,
                                matching_reason='PASS' if balanced else 'NO_JOINT_TARGET_OR_FIVE_CONTROLS'))
    for c in PAIRED_CONFIGS:
        jobs += [dict(config=c, sample_key=key, perturbed_gene=g,
                      role='candidate_reference' if g in PAIRED_TARGETS else 'matched_control') for g in sorted(perturb[c])]

for c in PAIRED_CONFIGS:
    folder = KO / 'matching' / c
    folder.mkdir(parents=True, exist_ok=True)
    for data, name in [(matches, 'matched_controls.tsv'), (quality, 'target_eligibility.tsv'), (jobs, 'perturbation_jobs.tsv')]:
        df = pd.DataFrame(data)
        df[df.config.eq(c)].to_csv(folder / name, sep='\t', index=False)
print('paired', len(jobs), 'knockouts')
