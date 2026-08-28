# Cross-mission clinically anchored renal tissue programs

**Status:** frozen retrospective targeted analysis  
**Lock date:** 2026-08-11  
**Machine-readable specification:** `config/clinical_renal_axes_cross_mission.yaml`

## Decision

This analysis is worth running. It asks a narrower and more translational
question than “do kidneys rewire,” while retaining a biological rather than a
methods-paper endpoint:

> Which clinically anchored renal tissue programs show directional engagement
> during terminal spaceflight, and which recur across independent mouse
> missions?

The study measures kidney-tissue RNA. It does **not** measure urine or blood
biomarkers, GFR, albuminuria, histologic fibrosis, morphology, or renal cancer.
Accordingly, a positive result can nominate prospective assays; it cannot claim
that urinary KIM-1, NGAL, DKK3, MMP-7, collagen fragments, or albumin actually
changed.

## Why the earlier recurrence numbers are not the answer

Two existing analyses overlap this question but are not valid confirmation.
The v11 recurrence code meta-analyzed genes and then Stouffer-combined gene
statistics, treating correlated genes as inferential replication. The Stage-0
effect-size analysis used animal-level scores but omitted OSD-462 and OSD-771,
included a live-return cohort with terminal cohorts, used small-sample-poor
meta-analysis, and pooled an invalid OSD-253 control comparison during its QC
gate. The often-quoted matrix `p=7e-4` and fibrosis `g=0.80` are therefore prior
design signals, not confirmatory results.

## Frozen four-axis family

All scores are oriented so positive values indicate the declared adverse tissue
program direction.

1. **Glomerular barrier identity loss.** Lower `Nphs1`, `Nphs2`, `Synpo`,
   `Ptpro`, `Magi2`, and `Wt1`. `Podxl` and `Cd2ap` enter only an expanded
   sensitivity because they are less podocyte-restricted.
2. **Tubular epithelial injury induction.** Higher `Havcr1` and `Lcn2`. Both
   genes must be observable and both are reported individually. `Il18` is
   secondary because epithelial and inflammatory sources vary.
3. **Fibrosis/maladaptive remodeling.** An equal-weight composite of (a)
   higher `Dkk3`, `Mmp7`, `Timp1`, and `Ccl2` plus lower `Egf`, and (b) higher
   collagen/ECM deposition genes. Both subdomains must retain the pooled
   direction.
4. **Distal transport identity loss.** An equal-weight DCT and ASDN transcript
   identity reference. This is transcript abundance, not transport flux or
   WNK-SPAK activity. It is multiplicity-corrected with the biological axes,
   but a distal-only result cannot make the paper a go.

The panels are external to the flight expression contrasts, but the analysis is
not preregistered: prior repository work had already revealed an ECM-like signal
and a heterogeneous distal program. This fact must remain visible in the paper.

## Primary mission family

The terminal synthesis comprises OSD-102, OSD-163, OSD-253, OSD-462, and the
ISS-terminal arm of OSD-771. OSD-253 uses only C57BL/6J animals and estimates
the approximately 25- and 75-day contrasts separately before combining them
within mission. OSD-462 uses wild-type total-RNA data; mRNA and UPX are paired
technical sensitivities. OSD-771 estimates young and old effects separately and
combines them within mission.

OSD-513 and the OSD-771 live-animal-return arm are recovery moderators, not
independent terminal replications. Age, duration, RNA preparation, and control
scenario never count as extra missions.

OSD-253 has an irreducible design caveat: its original ground controls differ
from flight in light-treatment annotation, while the white-light rerun controls
were sequenced in a different run/read-length context. The original same-run
ground control is primary and the white-light rerun is sensitivity. The exact
matched sample cells must pass the coverage diagnostic before inclusion.

## Primary estimator

For each mission, eligible genes are z-scored across the complete frozen
analysis sample set, assigned the frozen biological sign, and averaged within
subdomain. Composite axes average subdomains equally. Flight-versus-ground
effects are Hedges `g`. Age and duration strata are combined inside their
mission before cross-mission analysis.

Mission effects are synthesized using REML random effects, modified
Hartung-Knapp uncertainty, a prediction interval, and explicit heterogeneity.
The primary p-values come from at least 100,000 joint blocked label
permutations. Each permutation preserves OSD-253 duration and OSD-771 age,
recomputes all mission and meta statistics, and records the maximum absolute
statistic over exactly four axes. This creates one family-wise error-controlled
test without pretending that genes are independent cohorts.

### Label-blind lock amendment

Before calculating any new axis effect, the gene-observability audit showed
that the initially written `CPM >= 1` rule made HAVCR1/LCN2 unavailable in two
missions and removed most low-abundance repair sentinels despite finite OSDR
processed measurements. The rule was changed once, before examining group
directions, to `CPM >= 0.1 in at least half of either arm`. At these library
sizes that still requires roughly 5–15 expected counts. Features remaining near
zero, most notably `Mmp7` in several missions, are excluded mission by mission
and cannot be replaced. This amendment is machine-recorded in the YAML lock.

## Required adversarial checks

- exact contrast-level gene-body coverage and technical-balance audit;
- leave-one-mission and leave-one-gene analyses;
- mean versus median scores and a common-gene-intersection score;
- separate repair and ECM subdomain estimates;
- OSD-462 total-RNA versus mRNA and UPX;
- OSD-253 original versus rerun ground control;
- OSD-771 terminal versus live-return and young versus old;
- flight-versus-vivarium and basal-versus-ground context checks where valid;
- cell-composition adjustment as a different-estimand sensitivity;
- exploratory OSD-462 podocyte phosphoproteomics only if barrier RNA recurs
  without OSD-462.

## Interpretation and stop rule

An axis earns “recurrent tissue engagement” only if it clears the predeclared
family-wise test, small-sample interval, direction, heterogeneity, influence,
and score-definition gates in the YAML specification. The project is a paper
go only if the barrier, tubular-injury, or fibrosis/remodeling axis passes.

If no biological axis passes, this branch is retired. We will not promote an
individual gene, a recovery-only contrast, an uncorrected p-value, or whichever
alternative score looks best.
