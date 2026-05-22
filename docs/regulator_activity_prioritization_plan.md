# Regulator-Activity Prioritization — Pre-registration (v10)

**Status:** pre-registration · May 2026 · **Companion:** `docs/v10_route_assessment_and_plan.md`

## Purpose

Upgrade v9 from a state map around an established phenotype into an
**activity-anchored candidate-mechanism map**. Central question: *which
regulatory and signaling activity axes plausibly organize the recurring
DCT/NCC-low RNA state around the established spaceflight NCC/DCT phenotype?*

This is **computational mechanism prioritization, not causal discovery**. All
outputs use constrained wording: "candidate upstream organizer", "activity
anchor", "negative boundary" — never "mechanism" or "causal axis". Suggestive,
clearly-labelled results and ranked candidate lists are acceptable and intended
deliverables; the requirement is honest labelling of confidence, not certainty.

## Layers (trimmed from the original six)

**Layer A — kinase activity (KSEA).** Classic Casado et al. (2013) z-score
enrichment on the OSD-462 genome-wide phosphosite flight effects. Curated
renal core kinase-substrate net ships in-repo (`data/external/kinase_substrate/`);
a PhosphoSitePlus-format table can be supplied for the full panel.
*Positive control (pre-registered):* WNK and SPAK_OSR1 must return
`inferred_activity_down` at p < 0.05. **Outcome: passed** — SPAK_OSR1
z = −6.31 (p = 2.8×10⁻¹⁰), WNK z = −4.12 (p = 3.7×10⁻⁵).

**Layer B — TF / pathway activity.** decoupler ULM with PROGENy (pathway) and
DoRothEA/CollecTRI (TF) mouse priors, run **within each cohort** (RRRM-2 ISS-T,
RRRM-2 LAR, OSD-513, OSD-253, OSD-462 RNA) — never by raw cross-study pooling.
Cross-cohort recurrence is classified per regulator. Priors are network-fetched
at run time.

**Dropped from the original v10 plan:** perturbation-signature matching
(LINCS = human cell lines vs mouse DCT — generic, fishable, low defensible
yield) and the histology "audit" (no tissue/raw images available — it is
literature review, kept as a Discussion paragraph citing Siew et al., not a
layer). Rationale in `docs/v10_route_assessment_and_plan.md` §5.

## Integration

One mechanism-evidence ranking table grading candidate axes (WNK/SPAK/NCC;
S1P/adhesion/integrin; TGF-β/ECM; TLR/innate/macrophage; stress/preservation;
hypoxia/NRF2) on RNA recurrence, OSD-462 phospho support, and protein-abundance
boundary, with evidence grades: `activity_anchor`,
`candidate_upstream_organizer`, `candidate_context_axis`,
`technical_or_biological_caution`, `negative_boundary`. Failed/null axes are
retained, not dropped.

## Honest scope

- The kinase layer's top result (WNK-SPAK/OSR1 down) reproduces the Cosmic
  Kidney Disease NCC finding and is the positive control, not a discovery.
- Any kinase axis beyond WNK-SPAK/OSR1 is exploratory nomination.
- TF/pathway activity is RNA-level inference (regulator nomination), not
  measured protein activity.
- Mouse-site numbering in the curated kinase net is approximate for the WNK→
  SPAK/OSR1 sites; the NCC N-terminal SPAK/OSR1 cluster is well established.

## Code and reproducibility

- Module: `src/multiomics/regulator_activity.py`
- Orchestrator: `scripts/regulator_activity/run_regulator_activity.py`
- Tests: `tests/test_regulator_activity.py` (8 tests, offline machinery)
- Curated net: `data/external/kinase_substrate/renal_kinase_substrate_core.tsv`
- Layer A is fully offline and verified. Layer B requires one-time network
  access to fetch decoupler priors; the orchestrator prints the run recipe and
  skips gracefully if priors are unreachable.
