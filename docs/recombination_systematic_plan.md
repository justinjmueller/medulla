# Recombination systematic: findings and plan

Nov 3, 2026 · Mike Carrigan

The EMB recombination variation acts on electron energy almost entirely as a 4.4% charge-scale shift, and the plan is to replace that unconstrained shift with a π⁰-constrained shower energy-scale systematic (σ ≈ 1.6%) plus a small recombination residual.

## Findings

The EMB variation lowers reconstructed electron shower energy by 4.4% without changing the resolution. From the matched electron resolution fits: CV peak at +0.6%, variation at −3.8%, both with σ = 0.0334 and the same tails.

**Mechanism.** The CV and data use the modified box (MB) model. Reco divides shower charge by one flat R ≈ 0.68, the MB value at MIP and nominal field (SPINE `RecombinationCalibrator` applies its MIP factor to every non-track object). The variation simulates the ellipsoidal modified box (EMB) model, which yields about 4% less charge per MeV at shower-like dE/dx and angles. Reco is unchanged, so showers come out about 4% low.

**Matched-event test.** Scaling CV shower energy by 0.956 reproduces the variation's N_var/N_CV to about 1% in every reco energy bin, the first bin aside (+24% vs +21%). Matching on generator events moved the first bin from about +7% to +21%, so unmatched events cluster near threshold.

**Conclusion.** For electron energy, the recombination variation is a charge-scale shift. The non-scale part (dE/dx and angle dependence) is at the 1% level, about 3% in the first bin.

## Why the bands spike at 300 MeV and 1.5 GeV

The spikes come from the spectrum's slope at the bin edges, not from statistics. The detector splines are reco-binned weights in `reco_leading_electron_energy` (edges 300, 390, 490, 600, 730, 870, 1030, 1240, 1500, 2000, 5000 MeV) × opening angle. A fractional shift s moves s·E·n(E) events across an edge at E, where n is events per MeV, so a bin [a, b] changes by:

```latex
\frac{\Delta N}{N} \approx s\,\frac{b\,n(b) - a\,n(a)}{N}
```

| Bin | Reco-level change | Cause |
| --- | --- | --- |
| 300–390 MeV | about +21% | Sits on the selection threshold and a steep rise in the spectrum |
| 1.5–2 GeV | about −10% | Density falls from about 0.8 to 0.3 events/MeV across the bin: about 53 events leave at 1.5 GeV, 26 arrive from 2 GeV, −9.5% predicted |
| 2–5 GeV | about −11.5% | About 26 events leave at 2 GeV, almost none arrive from above 5 GeV |

Counts per bin look flat (about 285, 283, 258) only because the bins widen.

**Cross section.** The 1.5–2 GeV bin reaches about 18% in the cross section, against about 10% at reco level. In the GUNDAM forward fit, the dial lowers both the signal templates and the background there by about 10%. Data are fixed, so the signal parameter rises by about 10%, plus 10% × (1 − p)/p to fill the background gap: about 10%/p in total. A purity near 55–60% would give 18%. This is a candidate explanation, not yet tested, and it is not obviously consistent with the 2–5 GeV bin (12%).

## Which energy factors a variation should recompute

None of the factors in the shower energy formula should change for the EMB variation, so the 4.4% is a genuine prediction of the current chain. Rule: recompute in a variation only what data would recompute; anything derived from CV MC stays fixed.

```latex
E_{\mathrm{shower}} = W_{\mathrm{ion}} \cdot C_{\mathrm{cal}} \cdot C_{\mathrm{adj}} \cdot \frac{1}{R} \sum_i e^{t_{\mathrm{drift}}/\tau} Q_i
```

| Factor | Where it comes from | In the EMB variation |
| --- | --- | --- |
| C_cal ≈ 75.0 ± 1.1 e⁻/ADC | Joint fit of gain and recombination to stopping cosmic muons and protons in data (arXiv:2407.12969) | Keep. EMB changes how many electrons survive, not the electronics; the joint fit would return the simulated gain |
| R ≈ 0.68 | Flat MB value at MIP and nominal field | Keep. Fixed configuration, same as data |
| C_adj | CV MC truth: double-sided Crystal Ball means of reco/true per true-energy bin, averaged, inverted (π⁰ tech note, Appendix B). Covers subthreshold and clustering losses | Keep. MC-derived, unrelated to recombination |
| π⁰ shower scale (not in the analysis now) | Per-sample factor 135 MeV/c² ÷ fitted π⁰ peak, in data and simulation separately (π⁰ tech note, Sec. 2.2.3) | Would be refit per sample if adopted, since data go through it too |

**Gain.** A gain fitted with the recombination model held fixed would absorb about 3.5% and look different in the variation. That is not how ICARUS calibrates, so recalibrating the variation's gain would hide a real effect.

## Data constraint from the π⁰ mass peak

Data and CV MC photon showers agree to 0.8 ± 1.4%, so data disfavour a full 4.4% shower-scale shift at about 2 to 2.5σ. Numbers from the π⁰ tech note, Figure 2 (before the shower-scale adjustment; BNB Run 2, 1.62×10¹⁹ POT):

| Sample | Fitted π⁰ peak (MeV/c²) | Average error (MeV/c²) |
| --- | --- | --- |
| Data | 130.7 (+0.9 / −1.6) | ~1.25 |
| Simulation (CV, MB) | 131.7 (+1.1 / −1.4) | ~1.25 |

- **Gap:** 1.0 MeV/c², or 0.8%.
- **Error on the gap:** the two fit errors in quadrature, 1.3% to 1.5% depending on which asymmetric sides are used; 1.4% taken.
- **If data behaved like EMB:** the data peak would sit near 126 MeV/c².

**Proposed σ.** Gap and error in quadrature, which is standard and defensible. Using the 0.8% gap alone would be too aggressive.

```latex
\sigma = \sqrt{0.8^2 + 1.4^2} \approx 1.6\%
```

**What it constrains.** The net shower energy scale (gain, recombination, lifetime, clustering, C_adj), not recombination alone. That matches what the observable needs.

**Caveats.** BNB not NuMI; photon showers of a few hundred MeV, not electrons up to 5 GeV; fit errors are statistical only; the data/simulation labels in Figure 2 still need confirming. Redoing the fit on NuMI with this analysis's π⁰ selection addresses most of these.

## Plan

Split the recombination effect into a scale part and a residual, then replace only the scale part. The two factors multiply back to the original, so nothing is counted twice. N_shift is the CV with shower quantities scaled by s = 0.956.

```latex
\frac{N_{\mathrm{var}}}{N_{\mathrm{CV}}} = \underbrace{\frac{N_{\mathrm{var}}}{N_{\mathrm{shift}}}}_{\text{residual: keep}} \times \underbrace{\frac{N_{\mathrm{shift}}}{N_{\mathrm{CV}}}}_{\text{scale: replace}}
```

**Systematics after the change**

- **Shower energy scale (π⁰-constrained):** σ from the π⁰ data/MC comparison (about 1.6%), applied to shower quantities only. Built from CV, so it carries no spline noise from a separate sample.
- **Recombination model (residual):** the full residual as 1σ. Carries the dE/dx and angle effects, tracks and dE/dx-based particle ID.

**Steps**

1. Confirm the approach with the detector-systematics conveners. Ask the π⁰ analysis authors whether they refit their Sec. 2.2.3 factor for each detector variation.
2. Run the π⁰ selection on NuMI data, CV and the EMB variation, and fit each peak with a double-sided Crystal Ball. Include fit-model systematics: range, tail parameters, background shape.
    - Data vs CV gives the constraint and σ.
    - The variation's peak shift shows the probe is sensitive to recombination.
    - Comparing the π⁰ shift with the 4.4% electron shift checks that photons stand in for electrons.
3. Measure s per sample family (main and nue-enhanced `Recomb_Model`) from matched events.
4. Add a per-sample shower scale to the medulla selection. Use it to make:
    - the scale-removed recombination sample (factor 1/s);
    - CV at 1 ± σ for the new dial.

    Apply it to every shower calorimetric quantity: energies, the shower part of visible energy, and shower dE/dx if used. It must run inside the selection so the 300 MeV threshold and dE/dx cuts see it. Tracks stay untouched.
5. Replace `recombination` in the `[variations]` keys of both systematics TOMLs with the residual, add the shower-scale variation, and rerun phase 1, phase 2, merge B and stage 6 under a new project version.
6. In GUNDAM, remove the old recombination dial and add the two new ones with independent unit priors.
7. Validate:
    - closure: residual × 4.4% scale ratio reproduces the original N_var/N_CV in every bin;
    - fake-data fits with the EMB variation and with CV shifted by ±σ;
    - the signal-only vs background-only `applyCondition` split, to test the 10%/p mechanism;
    - bands before and after.

**Avoiding double counting**

- Keep the other shower-shifting variations (lifetime, YZ, noise) as they are; apply the π⁰ constraint only to recombination's scale part.
- Do not add the 1.5% gain uncertainty separately; the π⁰ comparison already includes the gain.
- Check whether an existing calorimetric or energy-scale systematic already shifts shower energy. If so, the new dial replaces it for showers.
- Make sure the NuMI π⁰ sample is orthogonal to the νe selection (NC π⁰ is a νe background).

## Expected impact and open questions

Going from a 4.4% to a 1.6% scale should cut the scale part of the recombination band to about 0.36 of its current size, assuming the band is linear in the shift (check this first with CV shifted by 1.6%).

| Bin | Now | Scale part at σ = 1.6% |
| --- | --- | --- |
| 300–390 MeV (reco) | ~21% | ~8% |
| 1.5–2 GeV (reco) | ~10% | ~3.6% |
| 2–5 GeV (reco) | ~11.5% | ~4% |
| 1.5–2 GeV (cross section) | ~18% | ~6.5% |

The residual adds in quadrature: about 3% in the first bin, about 1% elsewhere at reco level, possibly more in the cross section through the same background amplification. That puts recombination well below GENIE (~11%) and flux (~7%).

**Quick preview from existing samples.** Rebuild the band from CV shifted by 1.6%, combine it in quadrature with the matched residual, and plot it against the current recombination band.

**Open questions**

- [ ] Does the NuMI π⁰ fit agree as well as BNB? A significant gap (say 3 ± 1%) would mean correcting the central value with the per-sample π⁰ factor rather than widening σ.
- [ ] Did the π⁰ analysis refit its shower-scale factor per detector variation?
- [ ] Will the detector-systematics conveners accept a data-constrained σ?
- [ ] Why do unmatched events cluster near the 300 MeV threshold in the variation samples? This affects every detector variation.
- [ ] Does a flat σ hold from 300 MeV to 5 GeV? The reco/true response trend is about 1% across that range; a π⁰ fit in photon-energy bins would help.
- [ ] Is the 18% in the cross-section 1.5–2 GeV bin from background amplification (the 10%/p mechanism), or from the 2D energy × angle spline cells? Compare the reco-level band after phase 2 with the direct ratio.

## Sources

- [Calibration and simulation of ionization signal and electronics noise in the ICARUS LArTPC (arXiv:2407.12969)](https://arxiv.org/pdf/2407.12969): gain from a joint gain and recombination fit, 75.0 ± 1.1 e⁻/ADC; nominal MC uses the angle-independent ArgoNeuT modified box.
- [SPINE on PyPI](https://pypi.org/project/spine/) (v1.4.4, `spine/calib/recombination.py`): flat MIP recombination factor for non-track objects.
- ICARUS νμ CC π⁰ tech note (internal): Eq. 13 and Appendix B (C_adj), Sec. 2.2.3 and Figure 2 (π⁰ shower scale), Table 6 (recombination 3.6% on the π⁰ rate).
- This analysis: matched-event N_var/N_CV vs CV scaled by 0.956, electron resolution fits, `NuMI_nue_gundam.toml` spline binning.

EMB and MB parameter values used in rough estimates (ArgoNeuT α = 0.930, β = 0.212; ICARUS EMB α = 0.904, β = 0.204, R = 1.25) are approximate and should be checked against the simulation fcl and arXiv:2407.12969.
