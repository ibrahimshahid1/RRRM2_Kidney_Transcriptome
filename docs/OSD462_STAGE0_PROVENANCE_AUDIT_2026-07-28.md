# OSD-462 Stage 0 provenance audit

Stage 0 now reads the OSD-462 MS workbooks and ISA metadata directly and
produces an exact sample/channel table, multiplex-level raw-file inventory,
residue-resolved NCC/SPAK audit, QC table, and SHA-256 manifest. The manifest
also records a deterministic `git status --porcelain` listing and digest.

Run:

```bash
venv/bin/python scripts/osd462/08_stage0_provenance_audit.py
PYTHONPATH=. venv/bin/pytest -q tests/test_osd462_stage0.py tests/test_regulator_activity.py tests/test_osd462_anchor.py
```

Outputs are in `data/results/run_20260728_osd462_stage0/`.

## Blocking scientific corrections

1. The detailed assay protocol specifies TMTpro reagent, the +304.207-Da tag
   mass, and TMTpro impurity correction. A separate legacy investigation
   description says iTRAQ. These evidence sources are emitted in separate
   fields. The defensible assay name is **TMTpro isobaric-tag proteomics**,
   with the legacy inconsistency disclosed.
2. Both modalities contain 30 WT female B6129SF2/J left-kidney samples: two
   15-plexes with five baseline, five flight, and five ground samples each.
   The primary contrast is 10 flight versus 10 ground animals. The broader
   RR-10 p21-null arm is not part of these MS workbooks.
3. Condition is perfectly aligned with tag blocks in both plexes and the map
   is reused without a label swap. Thus condition cannot be separated from a
   systematic reporter-tag effect by the current design alone.
4. NCC workbook positions 65 and 68 are **Y65 and Y68**, not canonical
   SPAK/OSR1 sites. The strong `65;68` feature is `Y65;Y68`. Mixed features
   such as `53;65` and `58;65` are `T53;Y65` and `T58;Y65`; their signals
   cannot be attributed to the canonical threonine alone.
5. The strict mouse NCC sites are T53/T58/S71 (PMID:24039833). Only an exact
   T53 rollup is present. The strict SPAK WNK-target sites are T243/S383
   (PMCID:PMC2944316). Only an exact S383 rollup is present.
6. The T53 and S383 rollups both display co-modified peptide sequences
   (`T53+Y65` and `S382+S383`). They are residue-indexed site features, not
   isolated single-site occupancy measurements. Therefore the current
   assay-qualified isolated-canonical feature table contains zero rows.
   The audit emits source sheet/row, accession, motif, AScore, redundancy,
   modified and unmodified peptide sequence, motif-anchored absolute
   phosphoform labels, peptide depth, and reporter completeness for tracing.
7. The protein `Samp6-10` raw-file sheet repeats aliases `tc-885_11` and
   `tc-885_12`, although the acquisition IDs themselves are unique.

The corrected literature reference is
`data/external/kinase_substrate/renal_kinase_substrate_core.tsv`. Runtime site
lists and legacy plotting/compendium scripts were corrected so Y65/Y68, SPAK
S382, and position-only T53/S383 matches are not counted as isolated canonical
evidence.

No manuscript file was edited. The next reanalysis should treat isolated
canonical-site coverage as absent, report T53- and S383-indexed co-modified
features only as context, and empirically assess reporter-tag confounding.
