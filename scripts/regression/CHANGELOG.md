# Regression sandbox reference log

One entry per accepted reference run (`run_chain.py accept`), newest last. A run
is accepted only after the owner has looked at its differences from the previous
reference and said they are intended. The reference itself is
`/exp/icarus/data/users/micarrig/regression/reference/<tier>/current`; compare new
runs against it with `compare.py`.

## 2026-11-18 15:40 | full | 6b91c98d394c

- Reference: `/exp/icarus/data/users/micarrig/regression/reference/full/20261118-154015_6b91c98d394c`
- Code: 6b91c98; configs from `/exp/icarus/data/users/micarrig/regression/src_6b91c98d394c`
- Approved: initial reference: committed HEAD 6b91c98 (sandbox commit) with its own TOMLs (--configs-from-ref); no change to approve
- Diff none (first reference of this tier)

## 2026-11-18 16:00 | full | 6b91c98d394c-dirty-ppfxon

- Reference: `/exp/icarus/data/users/micarrig/regression/reference/full/20261118-160048_6b91c98d394c-dirty-ppfxon`
- Code: 6b91c98 + uncommitted PPFX changes (code.diff, code_untracked.txt in the reference, saved just after the run; see code_note.txt); configs from `/exp/icarus/app/users/micarrig/nuESpine/medulla`
- Approved: PPFX CV weight: true_ppfx_cv_weight branch in all selection trees; phase 1 PPFX-weighted (ppfx_from_truth); approved by owner 2026-11-18
- Diff vs 20261118-154015_6b91c98d394c: 23 files, 243 differences: RUNS DIFFER

```
239 x branch true_ppfx_cv_weight only in B
DIFF detsys_results/induction_high__reco_ele_beam_open_angle: 1 bin values differ (first bin 14: np.float64(0.013777941446432199) vs np.float64(0.01337811074613779))
DIFF detsys_results/induction_high__reco_ele_beam_open_angle: 1 bin errors differ (first bin 14: np.float64(0.007954698202976531) vs np.float64(0.007724117638637887))
DIFF detsys_results/induction_high__reco_leading_electron_energy: 3 bin values differ (first bin 1: np.float64(0.004592647148810733) vs np.float64(0.004511293476370932))
DIFF detsys_results/induction_high__reco_leading_electron_energy: 3 bin errors differ (first bin 1: np.float64(0.004592647148810733) vs np.float64(0.004511293476370932))
```
