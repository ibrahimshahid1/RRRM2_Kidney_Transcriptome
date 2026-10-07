# Local runbook: set up this project on a workstation and execute the research plan

**For:** a local Claude Code session on the owner's computer.
**Executes:** `docs/RESEARCH_PLAN_2026-10.md` (v1.1). Read this runbook fully first, then the plan.
**Written:** 2026-10-07, by the cloud session that produced G1/G2 and the plan.

## Why a local session

The cloud session ran in an ephemeral container. Nothing under `data/` survived it.

| Cloud limitation | What this runbook does locally |
|---|---|
| `data/` is gitignored and was lost on every container reset | Data persist on your disk; fetch once |
| The Ensembl REST gene map drifts with each Ensembl release | The exact gene map, tiers, KEGG library and atlas pseudobulk are **committed** in `resources/frozen/2026-10-06/`; restore them, don't rebuild |
| No R / DESeq2 | Install R + DESeq2 (needed for A6 heart VST) |
| Large downloads (OSD-102 proteomics, 3 GB) were impractical | Fetch locally, gated (§L5) |
| Some hosts were blocked (www.ncbi.nlm.nih.gov, pmc, bioRxiv) | Normal internet locally |
| Usage-limit interruptions lost agent work | Commit after every stage; at most 2–3 agents in flight |
| Outreach (emails) and the second-person ledger audit need a human | The session drafts; **the owner sends and audits** |

## Ground rules (non-negotiable)

1. **Never modify** frozen or locked artifacts:
   - `config/clinical_axes_podocyte_sensitivities.yaml`
   - `config/clinical_axes_recovery_persistence.yaml`
   - `config/gene_map_reconstruction.yaml`
   - `resources/frozen/**`
   - the G1/G2 results doc
   - any registration file after it is committed

   New analyses add files; they never rewrite frozen rows.
2. **Blinding.** The coordinator (you, the main session) is unblinded: you have read the plan's §2.2.
   Therefore:
   - **You do not write G3 gates and you do not write or run the A1–A4 code.** Fresh subagents do,
     given only `docs/registrations/G3_BLIND_BRIEF.md` (§L3).
   - Never paste observed numbers into a subagent prompt for A0 or A1–A4.
   - The gate-writer's YAML is adopted **verbatim**. The owner may reject it only for internal
     inconsistency, in writing, which triggers a new writer, not an edit.
3. **Predictions before unseen data.**
   - `scripts/fetch_inputs.py` refuses heart, brain and proteomics downloads until
     `docs/registrations/A6_PREDICTIONS.md` / `A8_PREDICTIONS.md` are committed.
   - A5 claims and predictions are committed before scoring.
4. **Commit and push after every stage.** Work on a branch named `local/plan-g3-<yyyymmdd>`.
   **Never merge to `main` without the owner's explicit approval.**
5. **Stop and report** at every checkpoint marked ⛔. Never push past a failed gate silently.
6. **Agent budget** (`docs/strategy/2026-10-07/README.md`):
   - at most 2–3 subagents at once;
   - Sonnet for reading, downloading and well-specified coding; Opus for adversarial review and
     statistical design;
   - every subagent gets an explicit file list, a word-limited return, and writes its output to a
     committed file.
7. **Untrusted downloads.** Keep downloaded files in their own directories under `data/`. Run Python
   that reads them with `python -I`. Never execute anything that came from a download.

## Phase L0: machine setup (⛔ report when done)

1. **OS.** On Windows, use **WSL2 (Ubuntu)**; paths and shell scripts assume POSIX. macOS and Linux
   work natively.
2. **Disk.** At least 10 GB free; more if OSD-102 proteomics is fetched.
3. **Repository.**
   - First time: `git clone https://github.com/ibrahimshahid1/RRRM2_Kidney_Transcriptome.git`, then
     `cd RRRM2_Kidney_Transcriptome && git checkout main && git pull`.
   - Confirm that `git log --oneline -1` is at or after `36a4eea` ("Research plan v1.1").
4. **Python 3.11**, isolated environment (conda/mamba or `python3.11 -m venv .venv`):
   ```bash
   pip install -r requirements-analysis.txt
   pip install -e .          # installs the `rrrm2` CLI; or use `python -m rrrm2`
   ```
5. **R ≥ 4.3 + DESeq2** (needed only in L5/A6). Install with conda (`r-base`,
   `bioconductor-deseq2`) or CRAN plus `BiocManager::install("DESeq2")`. Record `sessionInfo()` in
   `docs/registrations/R_SESSION.txt`.
6. **Tests.** Run `python -m pytest -q`. Expected baseline failures are **11 legacy tests** needing
   archived results or optional packages:
   - `test_external_aging_axis`
   - `test_osd462_stage0` (1)
   - `test_osd462_stage0_reporting` (2)
   - `test_regulator_activity`
   - `test_rrrm2_registry[network-rewiring]` (needs local `METHODOLOGY.md`/`FILE_SUMMARY.md`)
   - `test_v11_h2_enrichment` (3)
   - `test_v11_rna_protein_propagation` (2)

   Anything else failing → ⛔ stop.

## Phase L1: inputs (⛔ report when done)

1. **Restore the frozen resources:**
   ```bash
   python scripts/restore_frozen_resources.py      # must print "12/12 resources verified"
   ```
   This restores:
   - the Ensembl-116 gene map (`98e622bc…`);
   - the gene universe;
   - `KEGG_2019_Mouse.json` (`5a5cff30…`);
   - `frozen_compartment_tiers.tsv` (`f080a0d9…`) with its counts and manifest;
   - the atlas compartment expression, pseudobulk counts and sample metadata.
2. **Build the reconstructed baseline map:**
   ```bash
   python scripts/clinical_axes/build_reconstructed_id_map.py   # must end with sha256 d36462df…
   ```
3. **Fetch the primary OSDR inputs** (94 verified files; about 410 MB without the supplementary
   MultiQC bundles):
   ```bash
   python scripts/fetch_inputs.py                  # retries transient 5xx; extracts ISA bundles
   ```
   It must end with `done: 44/44 rows ok` (supplementary rows are skipped by default). OSD-771's
   `s_/a_/i_` tables land in `data/raw/metadata/`.
4. **Optional provenance check** (not required, because the tiers are restored):
   - Download the Mouse Kidney Atlas (`https://zenodo.org/api/records/17395591/files/mka.h5ad/content`,
     485 MB; MD5 must be `6a2e53086f100acffe99c054d3e8b7e9`) to
     `data/external/single_cell_atlases/mouse_kidney_atlas/mka.h5ad`.
   - Run `python -m rrrm2 run clinical-axes --stage atlas-pseudobulk --only`, then `marker-tiers`.
   - The tiers must reproduce `f080a0d9…`.
5. Run `python -m rrrm2 doctor`. clinical-axes must report no missing required inputs.

## Phase L2: reproduce the frozen results before anything new (⛔ hard gate)

```bash
RUN=data/results/clinical-axes/$(date -u +%Y%m%dT%H%M%SZ)_local_repro
for s in cross-mission compartment-context compartment-context-median compartment-context-common \
         podocyte-disjoint podocyte-scaffold-specificity podocyte-sensitivities atlas-loso-markers \
         atlas-loso-markers-idmap recovery-persistence recovery-persistence-idmap; do
  python -m rrrm2 run clinical-axes --run-dir "$RUN" --stage "$s" --only || break
done
```

Compare against the 2026-10-06 run. Point estimates and CIs must match to 1e-4. Permutation p-values
must match to 0.002 with identical package versions, or 3 Monte-Carlo SE otherwise.

| Output | Reference |
|---|---|
| `compartment_context_meta_results.tsv`, podocyte__high_specificity | g 0.68897 (0.041563, 1.336377); max-T 0.01917; I² 0 |
| same file, structural control | 0.5493 (−0.0947, 1.1933); max-T 0.1292 |
| `primary_meta_results.tsv` (4 axes) | −0.7163 / 0.5695 / 0.3115 / −0.1535 |
| `podocyte_scaffold_meta.tsv`, HC3 | P0 0.1983 (0.0230, 0.3736); P1 0.1123 (−0.0306, 0.2553); S1 0.0312 |
| `podocyte_scaffold_permutation.tsv`, P1 | legacy 0.0605; joint_label 0.0388; freedman_lane_joint 0.0415 |
| `podocyte_sensitivities/podocyte_sensitivity_gates.json` | overall ROBUST_DIRECTION; median 0.6311; common 0.6926 |
| `atlas_loso/atlas_loso_summary.tsv` | 8/8 gate pass; minimum retention 0.854 (Miao21) |
| `recovery_persistence/persistence_contrasts.tsv` | podocyte interaction −0.2806, FL p 0.0783, INCONCLUSIVE; structural_disjoint −0.3495, FL p 0.0088, INCONCLUSIVE |

The results doc lists sha256 prefixes of the 2026-10-06 outputs. Byte-identical hashes are a bonus,
not a requirement. **Any numeric mismatch → ⛔ stop and report; do not proceed to L3.**

## Phase L3: registration week (plan §5, A0; ⛔ owner checkpoints)

Create `docs/registrations/` and commit each item separately:

1. **`G3_BLIND_BRIEF.md`.** Copy plan Appendix A **verbatim** and nothing else.
2. **A0 gate-writer.** Spawn a *fresh* subagent (Opus) with this instruction: "Read only
   `docs/registrations/G3_BLIND_BRIEF.md` and the code paths it lists. Do not open `docs/` otherwise,
   `data/results/`, or `docs/strategy/`. Write `config/clinical_axes_g3_registration.yaml` as
   specified. Return only the file path."
   - Commit the YAML unchanged.
   - ⛔ Show it to the owner. The owner may reject only for internal inconsistency, which means
     re-spawning a new writer.
3. **A5 ledger** (coordinator plus a Sonnet reader, scoped to 15 searches): enumerate the claims (plan
   §6/A5) and the frozen mapping rule into `A5_LEDGER.md`.
   - Directional predictions come from the claims themselves, not from data.
   - ⛔ The owner (or a colleague) audits the claim → gene-set mapping. Commit after the audit.
4. **`A6_PREDICTIONS.md` and `A8_PREDICTIONS.md`.** Write directional predictions for H6a/H6b,
   the handling panels, and H8 per plan §6. Commit both. This unlocks the gated downloads.
5. **Outreach drafts** in `docs/outreach/` (the owner sends them):
   - the morphometry request (plan A9);
   - the GEO/author report of the GSE295428 duplicate upload (evidence:
     `docs/strategy/2026-10-07/GSE295428_FEASIBILITY.md`);
   - the OSDR query about RRRM-1 bulk kidney mRNA.
6. **A1–A4 executor.** Spawn a *fresh* subagent (Opus for A1's τ² model; Sonnet acceptable for A2–A4
   once A1 is specified). Give it only the blind brief, the registration YAML and the listed code
   paths.
   - It implements A1–A4 as new modules or stages, following the repo's conventions:
     - a preregistration guard (refuse an uncommitted YAML);
     - `--gene-map` dual runs;
     - synthetic tests;
     - stages registered in `config/revisions/clinical-axes.yaml` with stage-level `requires` and
       their own output folders.
   - It must **not** run on real data yet.
   - Commit, then write `docs/registrations/G3_HASHLOCK.json` (sha256 of the registration YAML and
     every new code file, plus the git HEAD). Commit.

## Phase L4: execute A1–A4 (⛔ report results)

1. The **executor** (the same blind subagent, or a new blind one) runs the hash-locked stages on real
   data into `data/results/clinical-axes/<UTC>_g3`. It logs every run, failed ones included, in
   `docs/registrations/G3_RUN_LOG.md`.
2. Only after completion, the coordinator opens the outputs, applies the registered labels
   **mechanically**, and determines the plan §7.2 row (R1–R4).
3. Write `docs/G3_RESULTS_<date>.md` with:
   - the provenance (lock commit, hashes, the known-outcome disclosure from plan §5.3);
   - per-set H1 labels, the H2 percentile, and the H4 table under both gene maps;
   - the §7.2 row.
4. One **Opus adversarial reviewer** reads the results doc against the registration YAML and the raw
   outputs. Fix the issues it raises, then commit.
5. ⛔ Report to the owner. Ask before merging to `main`.

## Phase L5: unseen-data analyses (one at a time; commit after each)

- **A5 ledger scoring.** Score the claims with the G1 machinery plus A1 calibration, two families with
  Holm, and verdicts per the registration. Write `docs/A5_LEDGER_RESULTS_<date>.md`.
- **A6 heart.**
  1. Fetch: `python scripts/fetch_inputs.py --manifest resources/plan_datasets.tsv --group a6_heart`
     (gated).
  2. Compute VST with DESeq2 `vst(blind=TRUE)` with GeneLab settings.
  3. **Validation gate first:** on one kidney mission that has a GeneLab VST, recompute VST from its
     counts and require per-gene Pearson ≥ 0.99 *and* 49-set mission g within ±0.02. Otherwise ⛔.
  4. Run H6a, H6b and the handling panels on the ISS-T arm only (LAR heart is right ventricle).
  5. Brain (`--group a6_brain`) is descriptive only.
- **A7 OSD-457.**
  - Fetch: `--group a7_mhu3` (not gated).
  - Load the 12 kidney samples of 192 through the project loaders.
  - Wild-type primary, Nrf2-KO separate. Supplement only.
- **A8 protein.**
  1. Fetch: `--group a8_osd163` (gated).
  2. Confirm the sample-ID linkage to the OSD-163 RNA animals.
  3. Coverage gate: ≥ 10 podocyte HS and ≥ 50 structural proteins, each in ≥ 80% of samples.
  4. Compare with the *same mission's* RNA effect using the CI rules.
  5. Commit `docs/registrations/A8_OSD163_RESULT.md`. Only if it is "informative", fetch `--group
     a8_osd102` (3 GB). ⛔ Ask the owner before that download.

## Phase L6: write-up and freeze (plan §9–§13)

Follow plan §9 (manuscript plan), §13 (housekeeping) and the week 6–7 freeze:
1. environment lock, CI, Docker cleanup and a README refresh;
2. tag `v1.0-manuscript`;
3. internal hostile review;
4. the Casaletto package.

⛔ Every merge to `main` and every external communication requires the owner's approval.

## Known gotchas

- **Do not rebuild the gene map from Ensembl REST.** It drifts. Restore it from `resources/frozen/`.
- **`gseapy`'s Enrichr call uses plain http.** The repo loader now fetches over HTTPS, but the KEGG
  file is restored anyway.
- **OSDR occasionally returns 502/503.** `fetch_inputs.py` retries; re-run it, since it skips files
  already present.
- **OSD-580 has no GeneLab VST** and uses the `GLDS-573_*` prefix. OSD-561/562/613 use
  `GLDS-556/557/589_*`.
- **OSD-771 dissection dates partly overlap across groups**; only the ERCC mix is confounded with
  group. Basal animals were dissected about 2 months earlier.
- **GSE295428** "Kidney, LAR" is a duplicate spleen upload. Do not analyse it.
- **The CPM eligibility rule is arm-aware.** Call it that, not "label-blind".
- **Do not install unpinned scanpy or anndata** (they pull numpy 2).

## Session hand-off (if this session ends mid-way)

Before stopping:
1. Commit and push.
2. Update `docs/registrations/STATUS.md` with the phase reached, the last commit, the next action and
   any open owner decision.

A new session resumes by reading `STATUS.md`, then this runbook.
