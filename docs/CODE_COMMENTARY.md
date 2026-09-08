# Archived code commentary

Verbose comment blocks removed from the source tree, preserved verbatim and
organized by filename. Each entry records the exact line range the block
occupied, the one-line comment that replaced it, and the original text.

**226 blocks** from **98 files**; 70 carry a claim boundary,
prohibition, or caveat and are marked `guardrail: true`. Their replacements
keep the prohibition explicit — a generic summary would invert several of them.

Where a rule provably cannot fit on one line, a **Preservation note** records
what the replacement points at rather than states. This document is the
relocation target for that text: it is the written record, and the source
one-liner is the pointer.

Line ranges are the pre-removal positions, verified byte-identical against the
tree immediately before the replacements were applied. They locate the block in
the history, not in the current file.

## Contents

- [`src/clinical_axes/data.py`](#srcclinicalaxesdatapy) — 2 block(s)
- [`src/clinical_axes/matching.py`](#srcclinicalaxesmatchingpy) — 1 block(s)
- [`src/clinical_axes/osd656_context.py`](#srcclinicalaxesosd656contextpy) — 3 block(s)
- [`src/clinical_axes/statistics.py`](#srcclinicalaxesstatisticspy) — 6 block(s)
- [`src/common.py`](#srccommonpy) — 5 block(s)
- [`src/data/build_id_map.py`](#srcdatabuildidmappy) — 2 block(s)
- [`src/enrichment/biological_grounding.py`](#srcenrichmentbiologicalgroundingpy) — 1 block(s)
- [`src/enrichment/gene_set_loader.py`](#srcenrichmentgenesetloaderpy) — 6 block(s)
- [`src/enrichment/regression_enrichment.py`](#srcenrichmentregressionenrichmentpy) — 1 block(s)
- [`src/grey60/adversarial.py`](#srcgrey60adversarialpy) — 5 block(s)
- [`src/markers/discover_dct.py`](#srcmarkersdiscoverdctpy) — 5 block(s)
- [`src/markers/discover_markers.py`](#srcmarkersdiscovermarkerspy) — 2 block(s)
- [`src/multiomics/celltype_panels.py`](#srcmultiomicscelltypepanelspy) — 4 block(s)
- [`src/multiomics/osd462_anchor.py`](#srcmultiomicsosd462anchorpy) — 7 block(s)
- [`src/multiomics/osd462_stage0.py`](#srcmultiomicsosd462stage0py) — 5 block(s)
- [`src/multiomics/phenotype_anchor.py`](#srcmultiomicsphenotypeanchorpy) — 5 block(s)
- [`src/multiomics/regulator_activity.py`](#srcmultiomicsregulatoractivitypy) — 5 block(s)
- [`src/networks/alternative_methods.py`](#srcnetworksalternativemethodspy) — 1 block(s)
- [`src/networks/contrast_vectors.py`](#srcnetworkscontrastvectorspy) — 10 block(s)
- [`src/networks/cross_osdr_projection.py`](#srcnetworkscrossosdrprojectionpy) — 1 block(s)
- [`src/networks/edge_regression.py`](#srcnetworksedgeregressionpy) — 2 block(s)
- [`src/networks/embeddings.py`](#srcnetworksembeddingspy) — 3 block(s)
- [`src/networks/external_aging_axis.py`](#srcnetworksexternalagingaxispy) — 1 block(s)
- [`src/networks/lioness.py`](#srcnetworkslionesspy) — 3 block(s)
- [`src/networks/mechanism_axis.py`](#srcnetworksmechanismaxispy) — 1 block(s)
- [`src/networks/procrustes.py`](#srcnetworksprocrustespy) — 2 block(s)
- [`src/networks/shared_topology.py`](#srcnetworkssharedtopologypy) — 5 block(s)
- [`src/networks/stability_test.py`](#srcnetworksstabilitytestpy) — 2 block(s)
- [`src/networks/wgcna_followup.py`](#srcnetworkswgcnafollowuppy) — 1 block(s)
- [`src/preprocessing/data_alignment.py`](#srcpreprocessingdataalignmentpy) — 4 block(s)
- [`src/preprocessing/deconvolution_sanity.py`](#srcpreprocessingdeconvolutionsanitypy) — 1 block(s)
- [`src/run_all_phases.py`](#srcrunallphasespy) — 5 block(s)
- [`src/statistics/direct_coexpression_test.py`](#srcstatisticsdirectcoexpressiontestpy) — 1 block(s)
- [`src/statistics/full_pipeline_permutation.py`](#srcstatisticsfullpipelinepermutationpy) — 1 block(s)
- [`src/statistics/full_regression.py`](#srcstatisticsfullregressionpy) — 2 block(s)
- [`src/statistics/interaction_metrics.py`](#srcstatisticsinteractionmetricspy) — 1 block(s)
- [`src/statistics/permutation_bootstrap.py`](#srcstatisticspermutationbootstrappy) — 1 block(s)
- [`src/statistics/silent_shifters.py`](#srcstatisticssilentshifterspy) — 1 block(s)
- [`src/subtype_reference/reference_builder.py`](#srcsubtypereferencereferencebuilderpy) — 3 block(s)
- [`src/v11/aldosterone_axis.py`](#srcv11aldosteroneaxispy) — 3 block(s)
- [`src/v11/cmap_screen.py`](#srcv11cmapscreenpy) — 2 block(s)
- [`src/v11/core_analysis.py`](#srcv11coreanalysispy) — 4 block(s)
- [`src/v11/dct_continuous_gradient.py`](#srcv11dctcontinuousgradientpy) — 2 block(s)
- [`src/v11/h2_composition_aware_phospho.py`](#srcv11h2compositionawarephosphopy) — 1 block(s)
- [`src/v11/human_concordance.py`](#srcv11humanconcordancepy) — 3 block(s)
- [`src/v11/kinome_atlas_ksea.py`](#srcv11kinomeatlaskseapy) — 7 block(s)
- [`src/v11/observability_audit.py`](#srcv11observabilityauditpy) — 6 block(s)
- [`src/v11/publication_figures.py`](#srcv11publicationfigurespy) — 1 block(s)
- [`src/v11/recurrence_meta.py`](#srcv11recurrencemetapy) — 6 block(s)
- [`src/v13/compartment_adversarial_audit.py`](#srcv13compartmentadversarialauditpy) — 1 block(s)
- [`src/v13/continuous_phospho_inference.py`](#srcv13continuousphosphoinferencepy) — 12 block(s)
- [`src/v13/reporting.py`](#srcv13reportingpy) — 1 block(s)
- [`src/validation/continuous_target.py`](#srcvalidationcontinuoustargetpy) — 1 block(s)
- [`src/validation/cross_validation.py`](#srcvalidationcrossvalidationpy) — 3 block(s)
- [`src/validation/enhanced_cv.py`](#srcvalidationenhancedcvpy) — 1 block(s)
- [`src/validation/external_replication.py`](#srcvalidationexternalreplicationpy) — 1 block(s)
- [`src/validation/multi_study_pool.py`](#srcvalidationmultistudypoolpy) — 1 block(s)
- [`src/validation/sample_features.py`](#srcvalidationsamplefeaturespy) — 5 block(s)
- [`src/validation/wgcna_external_validation.py`](#srcvalidationwgcnaexternalvalidationpy) — 1 block(s)
- [`src/visualization/network_diagnostics.py`](#srcvisualizationnetworkdiagnosticspy) — 1 block(s)
- [`scripts/align_metadata_to_counts.py`](#scriptsalignmetadatatocountspy) — 4 block(s)
- [`scripts/clinical_axes/run_compartment_context.py`](#scriptsclinicalaxesruncompartmentcontextpy) — 1 block(s)
- [`scripts/clinical_axes/run_podocyte_gene_coherence.py`](#scriptsclinicalaxesrunpodocytegenecoherencepy) — 1 block(s)
- [`scripts/clinical_axes/run_podocyte_podxl_nid1_sensitivity.py`](#scriptsclinicalaxesrunpodocytepodxlnid1sensitivitypy) — 1 block(s)
- [`scripts/clinical_axes/run_sensitivities.py`](#scriptsclinicalaxesrunsensitivitiespy) — 1 block(s)
- [`scripts/clinical_axes/run_strict_podocyte_matching_audit.py`](#scriptsclinicalaxesrunstrictpodocytematchingauditpy) — 1 block(s)
- [`scripts/map_top_hits.py`](#scriptsmaptophitspy) — 1 block(s)
- [`scripts/module_convergence_fisher.py`](#scriptsmoduleconvergencefisherpy) — 1 block(s)
- [`scripts/osd462/05_plot_dashboard.py`](#scriptsosd46205plotdashboardpy) — 1 block(s)
- [`scripts/osd462/09_stage0_manuscript_reporting.py`](#scriptsosd46209stage0manuscriptreportingpy) — 1 block(s)
- [`scripts/osd462/_common.py`](#scriptsosd462commonpy) — 1 block(s)
- [`scripts/phase3_node_rewiring_from_deltas.py`](#scriptsphase3noderewiringfromdeltaspy) — 2 block(s)
- [`scripts/phase3_procrustes_rewiring.py`](#scriptsphase3procrustesrewiringpy) — 3 block(s)
- [`scripts/phase4_anchor_qc_report.py`](#scriptsphase4anchorqcreportpy) — 1 block(s)
- [`scripts/phase5_derive_interaction_persistence.py`](#scriptsphase5deriveinteractionpersistencepy) — 1 block(s)
- [`scripts/phase7_grounding_fast.py`](#scriptsphase7groundingfastpy) — 1 block(s)
- [`scripts/pick_anchors.py`](#scriptspickanchorspy) — 2 block(s)
- [`scripts/plot_skeleton_diagnostics.py`](#scriptsplotskeletondiagnosticspy) — 1 block(s)
- [`scripts/regulator_activity/run_phenotype_anchor.py`](#scriptsregulatoractivityrunphenotypeanchorpy) — 1 block(s)
- [`scripts/regulator_activity/run_regulator_activity.py`](#scriptsregulatoractivityrunregulatoractivitypy) — 1 block(s)
- [`scripts/run_full_pipeline.py`](#scriptsrunfullpipelinepy) — 1 block(s)
- [`scripts/run_phase1_networks.py`](#scriptsrunphase1networkspy) — 1 block(s)
- [`scripts/run_phase2_pipeline.py`](#scriptsrunphase2pipelinepy) — 1 block(s)
- [`scripts/run_phase3_pipeline.py`](#scriptsrunphase3pipelinepy) — 1 block(s)
- [`scripts/run_phase4_7_pipeline.py`](#scriptsrunphase47pipelinepy) — 1 block(s)
- [`scripts/sf_classifier_component_decomp.py`](#scriptssfclassifiercomponentdecomppy) — 1 block(s)
- [`scripts/sf_classifier_negctrl_figure.py`](#scriptssfclassifiernegctrlfigurepy) — 1 block(s)
- [`scripts/spaceflight_kidney_classifier.py`](#scriptsspaceflightkidneyclassifierpy) — 1 block(s)
- [`scripts/stage0/axis_effect_size_anchor.py`](#scriptsstage0axiseffectsizeanchorpy) — 1 block(s)
- [`scripts/stage0/coverage_confounding_gate.py`](#scriptsstage0coverageconfoundinggatepy) — 1 block(s)
- [`scripts/stage0/protocol_inventory.py`](#scriptsstage0protocolinventorypy) — 1 block(s)
- [`scripts/subtype_reference/02a_extract_gse228367_raw_pseudobulk.py`](#scriptssubtypereference02aextractgse228367rawpseudobulkpy) — 1 block(s)
- [`scripts/v11/02_run_v11_core_analysis.py`](#scriptsv1102runv11coreanalysispy) — 1 block(s)
- [`scripts/v13/intensity_confound_audit.py`](#scriptsv13intensityconfoundauditpy) — 2 block(s)
- [`scripts/v13/layer_block_shift_comparison.py`](#scriptsv13layerblockshiftcomparisonpy) — 1 block(s)
- [`scripts/v13/reporter_position_diagnostic.py`](#scriptsv13reporterpositiondiagnosticpy) — 1 block(s)
- [`tests/test_continuous_phospho_inference.py`](#teststestcontinuousphosphoinferencepy) — 1 block(s)
- [`tests/test_osd462_anchor.py`](#teststestosd462anchorpy) — 1 block(s)

---

## `src/clinical_axes/data.py`

### Block 1 — module docstring, lines 1-6 — guardrail: true

**Replacement one-liner:**

```python
"""Load frozen cross-mission kidney RNA contrasts; no effect estimate selects any sample or gene."""
```

**Archived original:**

```python
"""Load and harmonize the frozen cross-mission kidney RNA contrasts.

All public functions in this module are label-transparent: sample membership is
derived only from the frozen configuration and repository metadata.  No effect
estimate is used to select a sample or gene.
"""
```

### Block 2 — comment run, lines 77-79 — guardrail: false

**Replacement one-liner:**

```python
        # OSDR DE-table annotations extend the 14k universe; the curated map wins on ID collisions.
```

**Archived original:**

```python
        # The OSDR differential-expression table supplies annotations for genes
        # outside the historical 14k repository universe.  The curated map
        # remains authoritative when both sources contain an ID.
```


---

## `src/clinical_axes/matching.py`

### Block 1 — function docstring, lines 59-66 — guardrail: true

**Replacement one-liner:**

```python
    """Draw unique matched panels within a balance caliper; no outcome or treatment label used."""
```

**Archived original:**

```python
    """Draw unique matched panels and enforce aggregate covariate balance.

    Each target draws from its ``pool_size`` closest candidates. Candidate
    genes cannot repeat within a panel. A draw is retained only when the
    absolute difference between every target and candidate panel-mean
    covariate is no larger than ``balance_caliper``. Inputs should already be
    standardized; no outcome or treatment label is used here.
    """
```


---

## `src/clinical_axes/osd656_context.py`

### Block 1 — module docstring, lines 1-8 — guardrail: true

**Replacement one-liner:**

```python
"""Descriptive OSD-656 urine context for frozen renal axes; never validation or meta-analysis."""
```

**Archived original:**

```python
"""Descriptive OSD-656 urine context for the frozen renal tissue axes.

OSD-656 contains post-flight urine measurements from four Inspiration4 crew
members.  This module intersects that assay with the *already frozen* mouse
kidney tissue-axis genes and preserves subjects and recovery timepoints.  It
does not perform a hypothesis test and its outputs must not be used as
validation or added to the mouse cross-mission meta-analysis.
"""
```

### Block 2 — function docstring, lines 26-29 — guardrail: true

**Replacement one-liner:**

```python
    """Flatten frozen primary-family subdomain genes only; no sensitivity or secondary markers."""
```

**Archived original:**

```python
    """Flatten only the frozen primary-family subdomain genes.

    Sensitivity additions and secondary markers are intentionally excluded.
    """
```

### Block 3 — function docstring, lines 129-134 — guardrail: true

**Replacement one-liner:**

```python
    """Pair recovery values to the subject's mean preflight baseline; timepoints never pooled."""
```

**Archived original:**

```python
    """Pair each recovery value to that subject's mean preflight baseline.

    Multiple preflight collections are first averaged within subject and
    analyte.  Recovery timepoints remain separate; they are never pooled or
    treated as independent subjects.
    """
```


---

## `src/clinical_axes/statistics.py`

### Block 1 — module docstring, lines 1-18 — guardrail: true

**Replacement one-liner:**

```python
"""Label-blind axis scoring and meta-inference primitives; the caller owns gene sign direction."""
```

**Archived original:**

```python
"""Statistical primitives for cross-mission renal tissue-axis analyses.

This module deliberately contains no repository paths, cohort names, gene
panels, or biological interpretation.  It provides the label-blind scoring
and inferential operations used by the clinical-axis runner:

* signed gene-wise z scores, with optional equal-weight subdomains;
* Hedges' g within an exchangeability stratum;
* inverse-variance fixed-effect pooling of strata within a mission;
* REML random-effects pooling across missions with modified
  Hartung--Knapp uncertainty and a prediction interval; and
* whole-pipeline blocked label permutations with max-|T| family-wise error
  control across a frozen set of axes.

Positive scores and effects always mean movement in the direction encoded by
the supplied gene signs.  The caller, not this module, owns that biological
direction convention.
"""
```

### Block 2 — function docstring, lines 127-133 — guardrail: true

**Replacement one-liner:**

```python
    """Z-standardize each gene across analysis samples; missing values stay missing, not neutral."""
```

**Archived original:**

```python
    """Z-standardize each gene across the supplied analysis samples.

    Missing observations remain missing and are not converted to a neutral
    value.  A gene with at least two finite observations but zero variance is
    assigned zero for its finite observations.  A gene with fewer than two
    finite values cannot be standardized and is left entirely missing.
    """
```

### Block 3 — function docstring, lines 241-252 — guardrail: false

**Replacement one-liner:**

```python
    """Signed label-blind per-sample axis score, equal-weighting subdomains before averaging."""
```

**Archived original:**

```python
    """Compute a signed, label-blind tissue-axis score for every sample.

    Genes are z-standardized separately across the samples in ``expression``
    before their frozen +/- direction is applied.  Without ``subdomains``, all
    observable genes receive equal weight.  With ``subdomains``, genes are
    first aggregated inside each subdomain and the subdomain scores are then
    averaged, ensuring (for example) that a large ECM list cannot outweigh a
    smaller maladaptive-repair list simply because it contains more genes.

    ``method='median'`` changes the within-domain gene aggregation; equal
    subdomain weighting remains an arithmetic mean across domains.
    """
```

### Block 4 — function docstring, lines 313-319 — guardrail: true

**Replacement one-liner:**

```python
    """Return Hedges' g and its sampling variance; non-finite values are rejected, not dropped."""
```

**Archived original:**

```python
    """Return Hedges' g and its conventional sampling variance.

    The small-sample correction is ``J = 1 - 3/(4*df - 1)`` and the sampling
    variance is ``(n1+n0)/(n1*n0) + g^2/(2*df)``.  Non-finite observations are
    rejected rather than silently removed because exchangeability-block sample
    counts are part of the frozen design.
    """
```

### Block 5 — function docstring, lines 454-461 — guardrail: true

**Replacement one-liner:**

```python
    """REML pooling with modified Hartung-Knapp scale floored at one; PI only for k >= 3."""
```

**Archived original:**

```python
    """REML random-effects meta-analysis with modified HK uncertainty.

    The Hartung--Knapp scale is floored at one (the modified HK rule), so a
    very homogeneous small set of missions cannot yield a standard error below
    the conventional random-effects standard error.  The prediction interval
    uses ``t_(k-2) * sqrt(tau^2 + SE_mKH^2)`` and is reported only for at least
    three missions.
    """
```

### Block 6 — function docstring, lines 760-772 — guardrail: false

**Replacement one-liner:**

```python
    """Whole-pipeline blocked max-|T| permutation over the frozen axis family, plus-one p-values."""
```

**Archived original:**

```python
    """Run a whole-pipeline, exchangeability-blocked max-|T| permutation.

    Treatment labels are shuffled independently inside each mission/stratum
    block while retaining its observed treatment count.  Every allocation is
    transformed into stratum-level Hedges g, fixed-effect mission estimates,
    and a cross-mission REML/mKH t statistic.  The maximum absolute t statistic
    across the supplied axes controls the single frozen axis family.

    The plus-one Monte Carlo correction is used for both unadjusted and max-T
    p-values.  Gene standardization and scoring must be completed before this
    function; because they are label-blind, they need not be recalculated under
    each permutation.
    """
```


---

## `src/common.py`

### Block 1 — module docstring, lines 2-10 — guardrail: false

**Replacement one-liner:**

```python
"""Shared RRRM-2 helpers: repo root, sample-column detection, label normalization, BH-FDR."""
```

**Archived original:**

```python
"""
Shared utilities for the RRRM-2 pipeline.

Centralizes functions that were previously copy-pasted across multiple modules:
  - REPO_ROOT: repository root path
  - find_sample_col: detect sample identifier column in metadata
  - normalize_labels: canonical Age/Arm/EnvGroup label normalization
  - bh_fdr: Benjamini-Hochberg FDR correction
"""
```

### Block 2 — function docstring, lines 25-29 — guardrail: false

**Replacement one-liner:**

```python
    """Find the OSD-771 sample identifier column in metadata, falling back to the first column."""
```

**Archived original:**

```python
    """Find the sample identifier column in metadata.

    Checks for common column names used across OSD-771 metadata files.
    Falls back to the first column if none match.
    """
```

### Block 3 — function docstring, lines 37-48 — guardrail: false

**Replacement one-liner:**

```python
    """Normalize to canonical Age (YNG/OLD), Arm (ISS-T/LAR), EnvGroup (FLT/GC/VIV/BSL) labels."""
```

**Archived original:**

```python
    """Normalize Age, Arm, and EnvGroup labels to canonical forms.

    Canonical forms:
        Age:      YNG, OLD
        Arm:      ISS-T, LAR
        EnvGroup: FLT, GC, VIV, BSL

    This is the single authoritative normalization used across all pipeline
    phases.  Previous versions in full_regression.py applied .str.upper()
    before replacement, which caused HGC to remain unmapped — that bug is
    fixed here.
    """
```

### Block 4 — function docstring, lines 75-84 — guardrail: false

**Replacement one-liner:**

```python
    """Benjamini-Hochberg FDR correction; returns q-values clipped to [0, 1], shaped like input."""
```

**Archived original:**

```python
    """Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    p : array-like of raw p-values

    Returns
    -------
    q : ndarray of adjusted p-values (same shape as input), clipped to [0, 1]
    """
```

### Block 5 — function docstring, lines 166-170 — guardrail: true

**Replacement one-liner:**

```python
    """Resolve configured symbols/IDs to Ensembl IDs and panel status; no rewiring stat is used."""
```

**Archived original:**

```python
    """Resolve configured gene symbols/IDs to Ensembl IDs and panel status.

    No observed rewiring statistic is consulted here.  Ambiguous symbols are kept
    as separate resolved rows; callers can decide whether ambiguity is acceptable.
    """
```


---

## `src/data/build_id_map.py`

### Block 1 — function docstring, lines 23-26 — guardrail: false

**Replacement one-liner:**

```python
    """Batch POST to Ensembl /lookup/id; returns {ensembl_id: {display_name, biotype, ...}}."""
```

**Archived original:**

```python
    """
    Batch POST to /lookup/id.
    Returns dict: {ensembl_id: {display_name, biotype, description, ...}}
    """
```

### Block 2 — function docstring, lines 68-71 — guardrail: false

**Replacement one-liner:**

```python
    """Look up a gene symbol via Ensembl xrefs; returns the matching Ensembl gene IDs."""
```

**Archived original:**

```python
    """
    Look up a gene symbol → Ensembl gene IDs via xrefs.
    Returns list of Ensembl gene IDs matching the symbol.
    """
```


---

## `src/enrichment/biological_grounding.py`

### Block 1 — module docstring, lines 2-20 — guardrail: false

**Replacement one-liner:**

```python
"""Phase 7: Fisher enrichment of pre-registered gene sets among top-decile rewiring genes."""
```

**Archived original:**

```python
"""
Phase 7: Fast biological grounding (poster-friendly).

1) Test enrichment of pre-registered gene sets among high-Δ genes (top decile)
   using Fisher's exact test.

Optional:
2) Cluster a reference embedding (k-means) and test cluster enrichment of high-Δ genes.

Inputs:
  - rewiring agg table (has gene + rewiring_mean)
  - optional embedding npy for module clustering
  - optional gene mapping TSV (Ensembl -> Symbol) for readability

Outputs:
  data/results/phase7_grounding/
    gene_set_enrichment.tsv
    cluster_enrichment.tsv (if embedding provided)
"""
```


---

## `src/enrichment/gene_set_loader.py`

### Block 1 — module docstring, lines 2-18 — guardrail: false

**Replacement one-liner:**

```python
"""Load Enrichr, .gmt, and curated-YAML gene sets, resolving symbols via the Ensembl id_map."""
```

**Archived original:**

```python
"""
Gene set loader for enrichment analysis.

Fetches gene set collections from Enrichr (via gseapy), caches locally
for reproducibility, and supports offline .gmt files.

Symbol resolution uses the Ensembl-backed id_map.json (built by
src.data.build_id_map via the Ensembl REST API) instead of the
old title-casing heuristic.

Usage:
    from src.enrichment.gene_set_loader import load_gene_sets

    sets, libs = load_gene_sets()
    # sets : dict[set_name → list[mouse_symbol]]
    # libs : dict[set_name → library_name]
"""
```

### Block 2 — function docstring, lines 45-49 — guardrail: false

**Replacement one-liner:**

```python
    """Load the Ensembl->Symbol map as a lowercase-symbol to canonical-mouse-symbol lookup."""
```

**Archived original:**

```python
    """Load the Ensembl→Symbol JSON map and build a case-insensitive
    symbol → canonical symbol lookup.

    Returns dict: lowercase_symbol → canonical_mouse_symbol
    """
```

### Block 3 — function docstring, lines 75-79 — guardrail: false

**Replacement one-liner:**

```python
    """Resolve a gene symbol to canonical mouse form via the id_map, falling back to the input."""
```

**Archived original:**

```python
    """Resolve a gene symbol to its canonical mouse form using the
    Ensembl-backed id_map.

    Falls back to the original string if not found.
    """
```

### Block 4 — function docstring, lines 92-96 — guardrail: false

**Replacement one-liner:**

```python
    """Load curated gene sets from config/gene_sets.yaml, resolving symbols via the id_map."""
```

**Archived original:**

```python
    """Load curated gene sets from config/gene_sets.yaml and resolve
    each symbol against the Ensembl-backed id_map.

    Returns: dict[set_name → list[resolved_mouse_symbol]]
    """
```

### Block 5 — function docstring, lines 167-173 — guardrail: false

**Replacement one-liner:**

```python
    """Load Phase 1.5b per-segment marker panels as 'segment_markers::SEGMENT' gene sets."""
```

**Archived original:**

```python
    """Load data-driven marker panels from Phase 1.5b (discover_markers.py).

    Each segment has a <SEGMENT>_marker_panel.txt file containing one gene
    per line.  These are loaded as 'segment_markers::SEGMENT' gene sets.

    Also checks the run-specific results directory.
    """
```

### Block 6 — function docstring, lines 267-289 — guardrail: false

**Replacement one-liner:**

```python
    """Load size-filtered Enrichr/.gmt/YAML gene sets; returns (gene_sets, set_to_library)."""
```

**Archived original:**

```python
    """
    Load gene sets from Enrichr libraries, .gmt files, and/or curated YAML.

    All symbols are resolved to canonical mouse form using the Ensembl-backed
    id_map.json (if available).  Falls back to raw strings if unavailable.

    Parameters
    ----------
    libraries : list of Enrichr library names (None → DEFAULT_LIBRARIES)
    gmt_files : list of paths to .gmt files
    cache_dir : directory for caching Enrichr downloads
    id_map_path : path to id_map.json (None → default location)
    curated_yaml : path to curated gene sets YAML (None → config/gene_sets.yaml)
    min_size, max_size : gene-set size filter (inclusive)
    include_curated : also include curated sets from YAML

    Returns
    -------
    gene_sets : dict[str, list[str]]
        set_name → list of mouse gene symbols
    set_to_library : dict[str, str]
        set_name → source library name
    """
```


---

## `src/enrichment/regression_enrichment.py`

### Block 1 — module docstring, lines 2-14 — guardrail: false

**Replacement one-liner:**

```python
"""Phase 7 supplement: gene-set enrichment over Phase 6 regression-significant genes."""
```

**Archived original:**

```python
"""
Phase 7 supplement: Run gene-set enrichment on regression-significant genes
(from Phase 6 full_regression) instead of rewiring-top-decile.

This tests: "Are genes with significant arm×flight or age×flight interaction
enriched for DCT/NCC-WNK or other predefined pathways?"

Usage:
    python -m src.enrichment.regression_enrichment \
        --reg_dir data/results/<run>/phase6_regression \
        --outdir data/results/<run>/phase7_regression_enrichment \
        --map data/processed/resources/id_map.tsv
"""
```


---

## `src/grey60/adversarial.py`

### Block 1 — module docstring, lines 1-5 — guardrail: false

**Replacement one-liner:**

```python
"""Core statistics for the frozen Grey60 adversarial reanalysis; the runner owns all I/O."""
```

**Archived original:**

```python
"""Core statistics for the frozen Grey60 adversarial reanalysis.

The functions here contain no repository path logic.  The runner is
responsible for loading the frozen inputs and writing manifests.
"""
```

### Block 2 — function docstring, lines 84-88 — guardrail: false

**Replacement one-liner:**

```python
    """Fixed-weight mean gene-z score normalized by the sum of absolute weights (scale-free)."""
```

**Archived original:**

```python
    """Return a fixed-weight mean gene-z score.

    Weights retain their frozen sign and are normalized by the sum of their
    absolute values, so an arbitrary scale change cannot alter the score.
    """
```

### Block 3 — function docstring, lines 193-197 — guardrail: false

**Replacement one-liner:**

```python
    """Contrasts over eight cell means ordered ISS-Y, ISS-O, LAR-Y, LAR-O, each as GC then FLT."""
```

**Archived original:**

```python
    """Contrasts over eight cell means.

    Cell order:
    ISS-Y GC/F, ISS-O GC/F, LAR-Y GC/F, LAR-O GC/F.
    """
```

### Block 4 — function docstring, lines 292-297 — guardrail: false

**Replacement one-liner:**

```python
    """Blocked max-|t| null over 11 contrasts x modules, permuted within four Age x Arm strata."""
```

**Archived original:**

```python
    """Generate a blocked max-|t| null across 11 contrasts x modules.

    Labels are permuted independently inside the four 5/5 Age x Arm strata.
    The implementation uses saturated cell means and a pooled 32-df residual
    variance, exactly matching the historical factorial OLS model.
    """
```

### Block 5 — function docstring, lines 382-388 — guardrail: true

**Replacement one-liner:**

```python
    """REML synthesis with HK scale; ``modified`` floors it at one so CIs cannot narrow."""
```

**Archived original:**

```python
    """Random-effects REML synthesis with Hartung-Knapp uncertainty.

    When ``modified`` is true, the Hartung-Knapp scale factor is floored at
    one. This prevents unusually homogeneous small meta-analyses from
    producing intervals narrower than the conventional random-effects
    interval.
    """
```


---

## `src/markers/discover_dct.py`

### Block 1 — function docstring, lines 27-34 — guardrail: false

**Replacement one-liner:**

```python
    """Build design matrix [intercept | CLR(DCT) | CLR(other) | tech | cells]; returns (X, idx)."""
```

**Archived original:**

```python
    """
    Build the full design matrix X [n_samples × p]:
      intercept | CLR(DCT) | CLR(other) | tech covariates | cell dummies

    Returns:
        X : (n, p) float64 design matrix
        dct_col_idx : column index of CLR(DCT) in X
    """
```

### Block 2 — function docstring, lines 108-118 — guardrail: false

**Replacement one-liner:**

```python
    """Vectorised OLS over all genes; returns beta_dct, t_stat, p_value, partial_r2 per gene."""
```

**Archived original:**

```python
    """
    Vectorised OLS: solve X @ B = Y in one shot.

    Args:
        X:   (n, p) design matrix
        Y:   (n, G) expression matrix (samples × genes)
        dct_col: index of CLR(DCT) column in X

    Returns:
        DataFrame with columns: beta_dct, t_stat, p_value, partial_r2
    """
```

### Block 3 — function docstring, lines 160-169 — guardrail: false

**Replacement one-liner:**

```python
    """Flag genes whose |beta_DCT| ranks top-N among segment coefficients (DCT is a minority)."""
```

**Archived original:**

```python
    """
    For each gene, check if |β_DCT| is among the top-N segment coefficients.

    DCT is a minority cell type in bulk kidney, so requiring it to be THE
    largest beta is too stringent (PT/TAL_LOH dominate). Top-3 ensures
    DCT markers are reasonably specific without being impossible.

    Returns:
        Boolean array (G,): True if β_DCT ranks in top-N segments.
    """
```

### Block 4 — function docstring, lines 192-197 — guardrail: false

**Replacement one-liner:**

```python
    """Stratified within-cell bootstrap fraction where beta_DCT > 0 and q < alpha."""
```

**Archived original:**

```python
    """
    Stratified bootstrap within experimental cells.
    Counts gene as passing when BOTH β_DCT > 0 AND q < alpha,
    matching the actual marker selection criterion.
    Returns: (G,) array of fraction of bootstraps where gene passes.
    """
```

### Block 5 — function docstring, lines 247-255 — guardrail: false

**Replacement one-liner:**

```python
    """Marginal Pearson bootstrap fallback used when the full OLS model is over-parameterised."""
```

**Archived original:**

```python
    """
    Fallback bootstrap using marginal Pearson correlation instead of OLS.
    Used when the full OLS model is over-parameterised (≈0 genes pass q<α).

    For each bootstrap resample, compute Pearson r(gene, DCT) and test
    whether r > 0 with BH-corrected p < α.

    Returns: (G,) array of fraction of bootstraps where gene passes.
    """
```


---

## `src/markers/discover_markers.py`

### Block 1 — function docstring, lines 30-33 — guardrail: false

**Replacement one-liner:**

```python
    """Build design matrix with the target segment as coefficient; returns (X, col_idx, names)."""
```

**Archived original:**

```python
    """Build design matrix with target segment as the coefficient of interest.

    Returns (X, target_col_idx, column_names).
    """
```

### Block 2 — function docstring, lines 199-202 — guardrail: false

**Replacement one-liner:**

```python
    """Discover marker genes for one segment; returns (results_df, panel_genes)."""
```

**Archived original:**

```python
    """Discover marker genes for a single segment.

    Returns (results_df, panel_genes).
    """
```


---

## `src/multiomics/celltype_panels.py`

### Block 1 — comment run, lines 29-30 — guardrail: false

**Replacement one-liner:**

```python
# Compartment marker proxies for the memo Section 8.2.3 split, reusing project mechanism sets.
```

**Archived original:**

```python
# Compartment scores used by the memo Section 8.2.3 split (re-uses project
# mechanism sets where available; symbols below are the marker proxies here).
```

### Block 2 — function docstring, lines 38-42 — guardrail: false

**Replacement one-liner:**

```python
    """Mean flight stat per panel; panels with < ``min_genes`` mapped report a NaN mean."""
```

**Archived original:**

```python
    """Mean flight stat per panel from a (gene symbol, stat) table.

    Returns one row per panel with the mean stat, n mapped genes, and the genes
    found. Panels with < ``min_genes`` mapped are still reported with NaN mean.
    """
```

### Block 3 — function docstring, lines 69-73 — guardrail: false

**Replacement one-liner:**

```python
    """Per-sample mean-z panel scores (panels x samples) from an ENSMUSG-indexed VST matrix."""
```

**Archived original:**

```python
    """Per-sample mean-z panel scores from an ENSMUSG-indexed VST matrix.

    Returns panels x samples. ``vst`` rows are ENSMUSG ids; ``sym_to_ens`` maps
    a lowercase symbol to its ENSMUSG id set.
    """
```

### Block 4 — function docstring, lines 87-90 — guardrail: false

**Replacement one-liner:**

```python
    """Classify bulk-RNA ambiguity per memo Section 8.3 using drop/stable mean-z thresholds."""
```

**Archived original:**

```python
    """Classify the bulk-RNA ambiguity per memo Section 8.3.

    ``drop``/``stable`` are flight-effect thresholds on the mean-z panel score.
    """
```


---

## `src/multiomics/osd462_anchor.py`

### Block 1 — function docstring, lines 40-43 — guardrail: false

**Replacement one-liner:**

```python
    """Map a row-3 machine header and row-2 label to channel metadata; None for metadata columns."""
```

**Archived original:**

```python
    """Map a row-3 machine header + row-2 sample label to channel metadata.

    Returns ``None`` for non-quantitative columns (metadata columns).
    """
```

### Block 2 — function docstring, lines 74-90 — guardrail: false

**Replacement one-liner:**

```python
    """Parse one TMT workbook sheet into a :class:`TmtTable` (``header_row`` is 1-based, 3 here)."""
```

**Archived original:**

```python
    """Parse one TMT workbook sheet into a :class:`TmtTable`.

    Parameters
    ----------
    path, sheet
        Workbook path and sheet name.
    gene_col
        Name of the gene-symbol column in the row-``header_row`` header.
    peptide_cols
        Mapping ``{PLEX1: <peptide col name>, PLEX2: <peptide col name>}``.
    id_col
        Optional protein-id column to retain.
    extra_meta_cols
        Additional metadata column names to retain (e.g. site position).
    header_row
        1-based index of the machine-header row (3 for these workbooks).
    """
```

### Block 3 — function docstring, lines 198-210 — guardrail: false

**Replacement one-liner:**

```python
    """Per-row Flight - Ground effect averaged across plexes, optionally channel-median centered."""
```

**Archived original:**

```python
    """Per-row Flight - Ground flight effect, estimated within each plex.

    For each plex the effect is ``mean(log2 FL) - mean(log2 GC)`` over the
    available channels; the two plex estimates are then averaged.  Optional
    per-channel med=ian centering (``channel_center``) performs standard TMT
    sample-loading normalization within each plex; because the differencing is
    within plex, both the per-plex normalization constant and (after centering)
    per-channel loading cancel.

    Returns one row per input row with the flight effect, per-plex effects,
    plex coverage, peptide counts, and a mean-abundance column used for
    matched-null binning.
    """
```

### Block 4 — function docstring, lines 260-268 — guardrail: false

**Replacement one-liner:**

```python
    """Per-row FL - GC effect, SE, CI and p from a per-site ``y ~ flight + plex`` linear model."""
```

**Archived original:**

```python
    """Per-row FL - GC effect with CI from a plex-adjusted linear model.

    For each row, log2 scaled FL/GC channels (both plexes) are fit with
    ``y ~ flight + plex``; the flight coefficient is the plex-corrected
    FL - GC effect.  By default, per-channel median centering performs the same
    within-plex TMT loading normalization used for protein effects.  Returns
    effect, SE, 95% CI, two-sided p, and channel counts.  Suitable for the
    small phosphosite tables where a per-site interval is wanted.
    """
```

### Block 5 — function docstring, lines 347-352 — guardrail: false

**Replacement one-liner:**

```python
    """Collapse protein rows to one row per gene by peptide-weighted mean of the flight effect."""
```

**Archived original:**

```python
    """Collapse multiple protein rows per gene symbol into one gene-level row.

    Many-to-one collisions (isoforms / shared symbols) are resolved by a
    peptide-weighted mean of the flight effect; peptide counts are summed and
    the number of collapsed rows is logged in ``n_protein_rows``.
    """
```

### Block 6 — function docstring, lines 452-459 — guardrail: false

**Replacement one-liner:**

```python
    """Stratum-matched random-gene-set null: each draw resamples the target's per-stratum counts."""
```

**Archived original:**

```python
    """Stratum-matched random-gene-set null for a gene-set statistic.

    ``pool`` is the background of all genes (with effects + match columns);
    ``target_mask`` selects the gene set; ``stat_fn`` maps a sub-frame of
    ``pool`` to a scalar; ``strata`` is the matched-sampling stratum per pool
    row.  Each null draw samples, within each stratum, the same number of genes
    the target set has there, then recomputes ``stat_fn``.
    """
```

### Block 7 — function docstring, lines 512-517 — guardrail: false

**Replacement one-liner:**

```python
    """Mean gene effect per pathway; pathways with < ``min_genes`` finite members are dropped."""
```

**Archived original:**

```python
    """Mean gene effect per pathway -> ordered pathway vector + coverage table.

    ``gene_effect`` is indexed by the same id space as the gene-set members
    (e.g. ENSMUSG).  Pathways with fewer than ``min_genes`` mapped, finite
    members are dropped from the vector.
    """
```


---

## `src/multiomics/osd462_stage0.py`

### Block 1 — module docstring, lines 1-12 — guardrail: true

**Replacement one-liner:**

```python
"""Stage 0 OSD-462 provenance from source workbooks and ISA metadata, never analysis outputs."""
```

**Archived original:**

```python
"""Reproducible Stage 0 provenance checks for the OSD-462 MS assays.

This module deliberately separates three objects that were previously mixed in
the analysis narrative:

* biological samples and their reporter channels;
* multiplex-level raw LC-MS acquisitions (which are not sample-specific);
* workbook phosphosite features and their residue identities.

The functions read the source workbooks and ISA metadata directly.  They do
not use manuscript tables or previously generated analysis outputs.
"""
```

### Block 2 — function docstring, lines 257-265 — guardrail: false

**Replacement one-liner:**

```python
    """Build the 60-row sample/modality design and the 96-row multiplex raw-acquisition table."""
```

**Archived original:**

```python
    """Build the exact sample design and multiplex-level raw inventory.

    Returns
    -------
    design
        One row per biological sample per MS modality (60 rows).
    raw_inventory
        One row per multiplex-level raw LC-MS acquisition (96 rows).
    """
```

### Block 3 — comment run, lines 338-339 — guardrail: true

**Replacement one-liner:**

```python
        # Keep protocol evidence and legacy metadata separate; never collapse into one label.
```

**Archived original:**

```python
        # Preserve protocol evidence and contradictory legacy metadata as
        # separate fields; do not collapse them into one ambiguous assay label.
```

### Block 4 — function docstring, lines 489-496 — guardrail: false

**Replacement one-liner:**

```python
    """Map every ``#`` in a workbook peptide sequence to an absolute site by 13-aa motif anchor."""
```

**Archived original:**

```python
    """Map every ``#`` in a workbook peptide sequence to an absolute site.

    The workbook encodes phosphorylation by placing ``#`` after the modified
    residue.  A 13-aa localization motif anchors its center to the reported
    position.  Longest-overlap alignment between each motif and the peptide
    determines the peptide's absolute start; independent component anchors
    must agree.
    """
```

### Block 5 — function docstring, lines 1074-1079 — guardrail: true

**Replacement one-liner:**

```python
    """Isolated canonical-site assay rows; intentionally empty for the current OSD-462 workbook."""
```

**Archived original:**

```python
    """Return assay rows qualified as isolated canonical-site features.

    The filter requires a strict residue-aware literature match, a single-site
    rollup, and exactly one ``#`` phosphomodification in the reported peptide
    sequence.  It is intentionally empty for the current OSD-462 workbook.
    """
```


---

## `src/multiomics/phenotype_anchor.py`

### Block 1 — comment run, lines 13-15 — guardrail: true

**Replacement one-liner:**

```python
# Literature site definition only; OSD-462 has no isolated canonical phosphoform at these sites.
```

**Archived original:**

```python
# Literature definition.  This is deliberately separate from OSD-462 assay
# qualification: the source workbook does not contain an isolated canonical
# phosphoform for any of these positions.
```

### Block 2 — comment run, lines 21-23 — guardrail: true

**Replacement one-liner:**

```python
# No OSD-462 row is an isolated canonical site: T53 co-occurs with pY65, S383 with pS382.
```

**Archived original:**

```python
# No OSD-462 row qualifies as an isolated canonical-site measurement after
# residue and peptide-phosphoform provenance are enforced.  The T53-indexed
# sequence also contains pY65; the S383-indexed sequence also contains pS382.
```

### Block 3 — comment run, lines 26-27 — guardrail: true

**Replacement one-liner:**

```python
# Showable as residue-indexed co-modified features only; never scored as isolated canonical sites.
```

**Archived original:**

```python
# These rows may be shown as supportive, residue-indexed *co-modified*
# features.  They must not be scored or labeled as isolated canonical sites.
```

### Block 4 — comment run, lines 33-35 — guardrail: true

**Replacement one-liner:**

```python
# Legacy position-only names resolve to the empty qualified set, so callers fail closed.
```

**Archived original:**

```python
# Backward-compatible names intentionally resolve to the empty qualified set.
# This makes legacy position-only code fail closed instead of silently treating
# T53/Y65 or S382/S383 phosphoforms as isolated regulatory-site evidence.
```

### Block 5 — function docstring, lines 88-91 — guardrail: false

**Replacement one-liner:**

```python
    """Mean z-scored value per sample over the requested features; missing features skipped."""
```

**Archived original:**

```python
    """Mean z-scored value across the requested features, per sample (column).

    ``values``: features x samples. Missing features are skipped (logged by
    caller via the returned coverage in n)."""
```


---

## `src/multiomics/regulator_activity.py`

### Block 1 — function docstring, lines 32-38 — guardrail: false

**Replacement one-liner:**

```python
    """Load a kinase-substrate table; requires kinase, substrate_gene, substrate_site columns."""
```

**Archived original:**

```python
    """Load a kinase-substrate table.

    Required columns: ``kinase``, ``substrate_gene``, ``substrate_site``.
    Optional columns (e.g. ``evidence``) are preserved. A PhosphoSitePlus
    ``Kinase_Substrate_Dataset`` export can be remapped to these column names
    upstream; the curated core panel ships in this format directly.
    """
```

### Block 2 — function docstring, lines 60-69 — guardrail: false

**Replacement one-liner:**

```python
    """Casado (2013) KSEA z-score; positive z = substrates more phosphorylated in flight."""
```

**Archived original:**

```python
    """Kinase-substrate enrichment analysis (Casado et al. 2013 z-score).

    For kinase ``k`` with quantified substrate-site set ``S``::

        z = (mean(effect_S) - mean(effect_all)) * sqrt(|S|) / sd(effect_all)

    and a two-sided normal p-value. Positive ``z`` => substrates collectively
    *more* phosphorylated in flight (inferred higher kinase activity); negative
    ``z`` => collectively *less* (inferred lower activity).
    """
```

### Block 3 — function docstring, lines 121-126 — guardrail: true

**Replacement one-liner:**

```python
    """Evaluate a predeclared down-direction control; the caller must first enforce provenance."""
```

**Archived original:**

```python
    """Evaluate a predeclared down-direction control *after* site qualification.

    This helper does not establish substrate identity. The caller must first
    enforce residue and phosphoform provenance; an unscored or under-covered
    control correctly returns ``False``.
    """
```

### Block 4 — function docstring, lines 142-151 — guardrail: false

**Replacement one-liner:**

```python
    """decoupler ULM activity from a contrasts x genes matrix and source/target/weight prior."""
```

**Archived original:**

```python
    """Univariate-linear-model activity inference via decoupler.

    ``gene_effects``: contrasts x genes matrix (one row per cohort contrast,
    columns are gene identifiers, values are flight-effect statistics).
    ``net``: long-format prior with columns ``source``, ``target``, ``weight``
    (PROGENy or DoRothEA/CollecTRI as returned by ``decoupler.op``).

    Returns a tidy table: contrast, source (TF/pathway), activity score, p, q.
    decoupler is imported lazily so this module is usable without it installed.
    """
```

### Block 5 — function docstring, lines 181-185 — guardrail: false

**Replacement one-liner:**

```python
    """Classify cross-cohort recurrence from per-cohort activity; requires a consistent sign."""
```

**Archived original:**

```python
    """Classify a regulator's cross-cohort recurrence from per-cohort activity.

    ``threshold`` is on the absolute activity score (decoupler ULM scores are
    approximately z-scaled). Recurrence requires a consistent sign.
    """
```


---

## `src/networks/alternative_methods.py`

### Block 1 — module docstring, lines 2-12 — guardrail: true

**Replacement one-liner:**

```python
"""Phase 8 sensitivity networks: SSN perturbation; CSN/scLink report unavailable, not LIONESS."""
```

**Archived original:**

```python
"""
Alternative sample-specific network methods for Phase 8 sensitivity.

Implemented:
  * SSN-style leave-one-sample perturbation of edge correlations.

Guarded wrappers:
  * CSN and scLink are exposed as explicit modes but require external
    implementations/dependencies; the benchmark records them as unavailable
    rather than silently substituting LIONESS.
"""
```


---

## `src/networks/contrast_vectors.py`

### Block 1 — function docstring, lines 47-55 — guardrail: false

**Replacement one-liner:**

```python
    """Compute beta, cos, rho and ||R|| for one (A_FLT, A_GC) pair, with optional weights."""
```

**Archived original:**

```python
    """Compute beta, cos, rho, ||R|| for one (A_FLT, A_GC) pair.

    Parameters
    ----------
    a_flt, a_gc : np.ndarray
        Flight and ground-control aging vectors. Must share length.
    weights : np.ndarray or None
        Optional non-negative per-feature precision weights (Guardrail C).
    """
```

### Block 2 — function docstring, lines 103-108 — guardrail: false

**Replacement one-liner:**

```python
    """Return R = A_FLT - beta * A_GC in the original basis; beta is fit in the weighted space."""
```

**Archived original:**

```python
    """Return the residual vector R = A_FLT - beta * A_GC.

    Note: beta is computed in the (possibly weighted) inner product space, but
    the residual is returned in the original coordinate basis so callers can
    interpret per-feature R values directly.
    """
```

### Block 3 — function docstring, lines 155-159 — guardrail: false

**Replacement one-liner:**

```python
    """Per-category bootstrap fractions; headline label only if it clears ``headline_fraction``."""
```

**Archived original:**

```python
    """Return per-category fractions and the headline assignment if any.

    The headline category is assigned only when its mass passes
    ``headline_fraction`` of the bootstrap distribution (§4.4).
    """
```

### Block 4 — function docstring, lines 192-196 — guardrail: false

**Replacement one-liner:**

```python
    """Manuscript-safe bootstrap category label; reports the mix unless one category dominates."""
```

**Archived original:**

```python
    """Return the manuscript-safe bootstrap category label (§4.4).

    A single headline category is used only when the configured fraction
    threshold is met. Otherwise the label explicitly reports the category mix.
    """
```

### Block 5 — function docstring, lines 213-217 — guardrail: false

**Replacement one-liner:**

```python
    """Resample indices with replacement within each unique stratum value."""
```

**Archived original:**

```python
    """Resample with replacement within each unique stratum value.

    Returns an array of indices into the original sample table with the same
    length as ``strata``.
    """
```

### Block 6 — function docstring, lines 258-263 — guardrail: false

**Replacement one-liner:**

```python
    """Return (A_GC, A_FLT) Old-minus-Young within-cell means, optionally restricted by arm."""
```

**Archived original:**

```python
    """Return (A_GC, A_FLT) = (Old-Young) in GC and FLT subsets.

    ``feature_matrix`` is samples x features. Aggregation is the within-cell
    mean. The optional ``arm_mask`` restricts to a single arm before computing
    the contrast.
    """
```

### Block 7 — function docstring, lines 291-295 — guardrail: false

**Replacement one-liner:**

```python
    """Per-feature precision weights 1 / (var_k + eps) from (B, K) bootstrap vectors."""
```

**Archived original:**

```python
    """Variance-based per-feature precision weights from bootstrap replicates.

    bootstrap_vectors is (B, K). Returns length-K weight vector w_k = 1 / (var_k + eps)
    where eps is the ``floor_percentile`` of the variance vector.
    """
```

### Block 8 — function docstring, lines 329-335 — guardrail: false

**Replacement one-liner:**

```python
    """Summarize a bootstrap decomposition into the pre-registered per-statistic artifact table."""
```

**Archived original:**

```python
    """Summarize a bootstrap decomposition distribution for output artifacts.

    The returned table matches the pre-registered artifact shape in
    agents_instruction.md §2.4: one row per statistic with point estimate,
    bootstrap median and percentile interval, interpretation-category fractions,
    and optional empirical permutation p-values.
    """
```

### Block 9 — function docstring, lines 385-390 — guardrail: false

**Replacement one-liner:**

```python
    """Bootstrap (beta, cos, rho); ``vector_builder(indices)`` rebuilds (A_GC, A_FLT) per draw."""
```

**Archived original:**

```python
    """Bootstrap the (beta, cos, rho) statistics.

    ``vector_builder(indices)`` must return ``(A_GC, A_FLT)`` recomputed on the
    resampled rows. ``strata`` is a per-sample stratum (e.g. concatenation of
    Age, Arm, EnvGroup labels) governing the stratified resample.
    """
```

### Block 10 — function docstring, lines 428-432 — guardrail: true

**Replacement one-liner:**

```python
    """Permute Old/Young within strata and recompute; ``strata`` must exclude the permuted label."""
```

**Archived original:**

```python
    """Permute Old/Young labels within supplied strata and recompute stats.

    ``strata`` should exclude the permuted label itself. If omitted, the legacy
    behavior stratifies by EnvGroup plus the supplied arm mask.
    """
```


---

## `src/networks/cross_osdr_projection.py`

### Block 1 — function docstring, lines 127-135 — guardrail: false

**Replacement one-liner:**

```python
    """FLT/GC label-permutation null for cosine alignment; the reference vector is held fixed."""
```

**Archived original:**

```python
    """Build a FLT/GC label-permutation null for cosine alignment.

    ``external_vector_builder(None)`` must return the observed external
    flight vector. ``external_vector_builder(permuted_labels)`` must rebuild
    the external vector after assigning the provided permuted FLT/GC labels to
    the same samples. The reference vector is held fixed because the question
    is whether the external FLT/GC contrast aligns with the RRRM-2 direction
    more than expected under exchangeable external labels.
    """
```


---

## `src/networks/edge_regression.py`

### Block 1 — module docstring, lines 2-18 — guardrail: true

**Replacement one-liner:**

```python
"""Phase 2 A2-A3: limma edge regression on LIONESS weights; covariates only if unresidualized."""
```

**Archived original:**

```python
"""
Phase 2 Step A2-A3: Edge-wise Regression + Predicted Networks.

Fits limma edge-wise regression on LIONESS edge weights. In the primary
pipeline the expression input has already been residualized for technical
covariates, deconvolution, and selected SVs, so the edge model uses only
Age × Arm × EnvGroup cell means. Nuisance covariates are allowed only when the
edge weights were built from a non-residualized expression source.

Outputs:
  - Contrast effects (Δz) for rewiring analysis
  - Predicted condition-specific networks (z_hat)
  - limma topTable results per contrast

Usage:
    python -m src.networks.edge_regression
"""
```

### Block 2 — comment run, lines 309-310 — guardrail: false

**Replacement one-liner:**

```python
        # Save the fit2$coefficients edge contrast under the legacy filename for compatibility.
```

**Archived original:**

```python
        # Save edge-weight contrast from fit2$coefficients (legacy filename
        # keeps downstream compatibility).
```


---

## `src/networks/embeddings.py`

### Block 1 — module docstring, lines 2-19 — guardrail: false

**Replacement one-liner:**

```python
"""Phase 3.2: multi-seed PecanPy node2vec embeddings on Phase 2 predicted edge-weight networks."""
```

**Archived original:**

```python
"""
Phase 3.2: PecanPy node2vec embeddings on predicted edge-weight networks.

Builds biased random walks + trains embeddings using PecanPy for speed.

Graph:
  - fixed topology from Phase 2 skeleton (edge_i, edge_j)
  - primary signed mode embeds positive and negative edge-weight channels
    separately, then concatenates them per gene

Targets:
  - by default embed only FLT/GC predicted networks via --patterns
  - multi-seed for robustness; saves embeddings per seed + stability stats

Usage:
  python -m src.networks.embeddings
  python -m src.networks.embeddings --num_seeds 2 --num_walks 20 --walk_length 40  # quick test
"""
```

### Block 2 — function docstring, lines 41-45 — guardrail: false

**Replacement one-liner:**

```python
    """Write the tab-separated weighted edgelist PecanPy expects, using integer node IDs."""
```

**Archived original:**

```python
    """
    Write weighted edgelist expected by PecanPy:
      node_u <tab> node_v <tab> weight
    We use integer node IDs as strings to keep things compact.
    """
```

### Block 3 — function docstring, lines 59-62 — guardrail: false

**Replacement one-liner:**

```python
    """Call PecanPy ``g.embed`` with only the keyword arguments the installed version supports."""
```

**Archived original:**

```python
    """
    Call PecanPy g.embed(...) but only with args supported by the installed version.
    This prevents crashes like: unexpected keyword argument 'workers' or 'seed'.
    """
```


---

## `src/networks/external_aging_axis.py`

### Block 1 — function docstring, lines 117-123 — guardrail: false

**Replacement one-liner:**

```python
    """Export the TMS kidney female old-minus-young donor-pseudobulk log2(CPM+1) aging axis."""
```

**Archived original:**

```python
    """Create the TMS kidney female old-vs-young aging-axis TSV.

    Cells are first summed to donor-level pseudobulk counts, then converted to
    log2(CPM + 1). The exported effect is the old-donor mean minus the young-
    donor mean for each gene. Age ranges are inclusive when a two-value range is
    supplied; a single-value list such as ``[3]`` selects exactly 3-month donors.
    """
```


---

## `src/networks/lioness.py`

### Block 1 — module docstring, lines 2-14 — guardrail: true

**Replacement one-liner:**

```python
"""Phase 2 A1: LIONESS weights on skeleton E, rank-normalized rather than Fisher-z transformed."""
```

**Archived original:**

```python
"""
Phase 2 Step A1: Compute LIONESS-style sample-specific edge weights on skeleton E.

Uses the residualized, non-cell-standardized Rtech for sample-specific edge
weights. The primary output is a raw LIONESS correlation contribution that is
rank-normalized edge-wise across samples before regression. Fisher z is applied
only to pooled/leave-one-out correlation estimates in the explicit
``z_contribution`` sensitivity mode; sample-specific LIONESS weights are not
treated as correlations and are not Fisher transformed in the default mode.

Usage:
    python -m src.networks.lioness --lioness-transform raw_ranknorm
"""
```

### Block 2 — function docstring, lines 50-54 — guardrail: false

**Replacement one-liner:**

```python
    """Rank-normalize each edge across samples to normal scores, the default regression scale."""
```

**Archived original:**

```python
    """Rank-normalize each edge across samples to normal scores.

    This is the default regression scale for LIONESS contributions because the
    sample-specific values are linear contributions, not bounded correlations.
    """
```

### Block 3 — function docstring, lines 95-101 — guardrail: true

**Replacement one-liner:**

```python
    """LIONESS edge weights for genes x samples; ``z_contribution`` is a sensitivity mode only."""
```

**Archived original:**

```python
    """Compute LIONESS edge weights for X (genes x samples).

    ``raw*`` modes apply the LIONESS linear formula on Pearson correlations and
    then optionally normalize sample-specific contributions. ``z_contribution``
    applies the linear formula after Fisher-z-transforming the pooled and LOO
    correlations; it is retained only as a sensitivity analysis.
    """
```


---

## `src/networks/mechanism_axis.py`

### Block 1 — function docstring, lines 118-123 — guardrail: false

**Replacement one-liner:**

```python
    """Resolve symbol/Ensembl YAML sets to scored Ensembl IDs, keeping ``protected_genes``."""
```

**Archived original:**

```python
    """Resolve symbol/Ensembl YAML sets to scored Ensembl IDs.

    ``deoverlap_gene_ids`` is used for candidate-upstream mechanism scores so
    axis-defining genes can be reported and removed.  Symbols listed in a set's
    ``protected_genes`` field remain scored even if they overlap the axis.
    """
```


---

## `src/networks/procrustes.py`

### Block 1 — module docstring, lines 2-13 — guardrail: false

**Replacement one-liner:**

```python
"""Phase 3.3: multi-seed Procrustes alignment and 1 - cosine rewiring over configured anchors."""
```

**Archived original:**

```python
"""
Phase 3.3: Procrustes alignment + cosine-distance rewiring (multi-seed)

- Loads configured anchors from config/anchor_genes.yaml
- Resolves symbols to Ensembl IDs through id_map.tsv
- Aligns embeddings per seed via orthogonal Procrustes
- Computes rewiring = 1 − cosine_similarity
- Reports seed mean/std + rank variance

Usage:
    python -m src.networks.procrustes
"""
```

### Block 2 — function docstring, lines 88-92 — guardrail: false

**Replacement one-liner:**

```python
    """Rotate B onto A by orthogonal Procrustes; returns (B_aligned, R)."""
```

**Archived original:**

```python
    """
    Rotate B to best match A using orthogonal Procrustes:
      find R minimizing ||B R - A||_F
      return (B_aligned, R)
    """
```


---

## `src/networks/shared_topology.py`

### Block 1 — module docstring, lines 2-18 — guardrail: false

**Replacement one-liner:**

```python
"""Phase 2.1-2.3: build the shared sparse skeleton E by Ledoit-Wolf top-k partial correlation."""
```

**Archived original:**

```python
"""
Phase 2 Step 2.1-2.3: Build Shared Sparse Skeleton E

Cell-standardize expression within each (Age × Arm × EnvGroup) cell,
then build a sparse partial correlation network using Ledoit-Wolf
shrinkage and top-k neighbors per gene.

The skeleton E is fixed for all downstream sample-specific weighting.

Usage:
    python -m src.networks.shared_topology --max_genes 2500 --topk 80

With biotype filtering (recommended):
    python -m src.networks.shared_topology --max_genes 2500 --topk 80 \\
        --id_map data/processed/resources/id_map.tsv \\
        --biotype_filter protein_coding
"""
```

### Block 2 — function docstring, lines 45-48 — guardrail: false

**Replacement one-liner:**

```python
    """Load ensembl_gene_id, mgi_symbol, biotype from id_map.tsv."""
```

**Archived original:**

```python
    """Load Ensembl biotype annotations from id_map.tsv.

    Returns DataFrame with columns: ensembl_gene_id, mgi_symbol, biotype
    """
```

### Block 3 — function docstring, lines 130-146 — guardrail: false

**Replacement one-liner:**

```python
    """Select up to ``max_genes`` HVGs after biotype/noise filtering, force-included genes first."""
```

**Archived original:**

```python
    """Select top genes by variance, with biotype filtering and force-include.

    Pipeline:
        1. Drop ERCC spike-ins
        2. (Optional) Filter to allowed biotypes via id_map.tsv
        3. (Optional) Exclude Gm-prefix and other noise symbols
        4. Force-include genes get priority slots
        5. HVGs fill remaining slots so total ≤ max_genes

    Args:
        rtech: Expression matrix (genes x samples), index = Ensembl IDs
        max_genes: Maximum total genes in panel
        force_include: Ensembl IDs to force-include (bypass biotype filter)
        biotype_map: DataFrame with ensembl_gene_id, mgi_symbol, biotype
        allowed_biotypes: List of allowed biotypes (e.g. ["protein_coding"])
        exclude_noise_symbols: If True, drop Gm\\d+, Rik, etc. from HVG pool
    """
```

### Block 4 — function docstring, lines 219-231 — guardrail: false

**Replacement one-liner:**

```python
    """Standardize expression within each experimental cell; returns n_samples x n_genes."""
```

**Archived original:**

```python
    """
    Standardize within each experimental cell (defined by cell_cols).
    
    Args:
        rtech_gxs: genes x samples DataFrame
        meta: metadata with index = sample IDs
        cell_cols: columns defining experimental cells
        eps: small constant for numerical stability
        sd_floor: minimum SD to avoid division issues
    
    Returns:
        Z: (n_samples x n_genes) cell-standardized matrix
    """
```

### Block 5 — function docstring, lines 264-268 — guardrail: false

**Replacement one-liner:**

```python
    """Build the skeleton from the union of each gene's top-k neighbors (~G*k edges, not G*k/2)."""
```

**Archived original:**

```python
    """
    Build skeleton by taking top-k neighbors per gene.
    
    The union of all top-k neighbors gives ~G*k edges (not G*k/2).
    """
```


---

## `src/networks/stability_test.py`

### Block 1 — function docstring, lines 30-35 — guardrail: false

**Replacement one-liner:**

```python
    """Read stability gates from config/contrast_vector_framework.yaml into median/lower keys."""
```

**Archived original:**

```python
    """Extract stability gates from config/contrast_vector_framework.yaml.

    The YAML uses explicit keys (``stability_median`` / ``stability_lower``)
    while the in-memory gate application uses the shorter ``median`` /
    ``lower`` names.
    """
```

### Block 2 — function docstring, lines 128-132 — guardrail: false

**Replacement one-liner:**

```python
    """Bootstrap-vs-full angular stability of A^a_GC (``agc_builder(None)`` = full sample)."""
```

**Archived original:**

```python
    """Compute bootstrap-vs-full-sample angular stability of A^a_GC.

    ``agc_builder(indices)`` should rebuild A^a_GC on the resampled rows.
    Calling it with ``None`` returns the full-sample (no-resample) vector.
    """
```


---

## `src/networks/wgcna_followup.py`

### Block 1 — module docstring, lines 2-5 — guardrail: false

**Replacement one-liner:**

```python
"""WGCNA follow-up: contrasts, hub genes, eigengene plots, GO/KEGG enrichment, integrated table."""
```

**Archived original:**

```python
"""
WGCNA Manuscript Follow-up: Simple contrasts, hub genes, eigengene plots,
GO/KEGG enrichment, and final integrated table.
"""
```


---

## `src/preprocessing/data_alignment.py`

### Block 1 — function docstring, lines 26-30 — guardrail: false

**Replacement one-liner:**

```python
    """Read a GeneLab counts CSV into a DataFrame indexed by gene_id with sample columns."""
```

**Archived original:**

```python
    """
    Reads a GeneLab counts CSV shaped like:
      gene_id, sample1, sample2, ...
    Returns a DataFrame indexed by gene_id with sample columns.
    """
```

### Block 2 — function docstring, lines 62-65 — guardrail: false

**Replacement one-liner:**

```python
    """Align metadata to the canonical sample ordering, raising if any sample is missing."""
```

**Archived original:**

```python
    """
    Align metadata to the canonical sample ordering.
    Raises error if any samples are missing from metadata.
    """
```

### Block 3 — function docstring, lines 75-78 — guardrail: false

**Replacement one-liner:**

```python
    """Align count-matrix columns to the canonical sample order, raising if any is missing."""
```

**Archived original:**

```python
    """
    Align count matrix columns to the canonical sample ordering.
    Raises error if any samples are missing from count matrix.
    """
```

### Block 4 — function docstring, lines 90-95 — guardrail: false

**Replacement one-liner:**

```python
    """Parse EnvGroup, Arm, Age, AnimalCode out of an RRRM2 sample name."""
```

**Archived original:**

```python
    """
    Optional: extract simple design factors from the RRRM2 sample name.
    Example: RRRM2_R-KDN_BSL_ISS-T_YNG_BY1
    
    Returns dict with keys: EnvGroup, Arm, Age, AnimalCode
    """
```


---

## `src/preprocessing/deconvolution_sanity.py`

### Block 1 — module docstring, lines 1-4 — guardrail: false

**Replacement one-liner:**

```python
"""Validate Phase 0 deconvolution by correlating CLR cell fractions with canonical marker VST."""
```

**Archived original:**

```python
"""
Script: src/preprocessing/deconvolution_sanity.py
Purpose: Validate deconvolution (Phase 0) by correlating cell type fractions (CLR) with canonical marker expression (VST).
"""
```


---

## `src/run_all_phases.py`

### Block 1 — function docstring, lines 37-47 — guardrail: false

**Replacement one-liner:**

```python
    """Resolve an artifact from the current run, the latest prior run, then the legacy path."""
```

**Archived original:**

```python
    """Resolve an artifact from the current run or the most recent prior run.

    Search order:
      1. Current run's results dir  (data/results/<current_run>/<relpath>)
      2. Most recent prior run       (data/results/<latest_run>/<relpath>)
      3. Legacy global location      (data/processed/<relpath>)

    *phase_subdir* is unused but kept for future sub-folder logic.

    Returns the first existing Path, or None.
    """
```

### Block 2 — function docstring, lines 71-83 — guardrail: false

**Replacement one-liner:**

```python
    """Create the versioned data/results/<run_id>/ phase directories and write run_metadata.json."""
```

**Archived original:**

```python
    """Initialize versioned output directories and save run metadata.

    All pipeline outputs live under a single run directory:
        data/results/<run_id>/
            deconvolution/     ← Phase 0
            phase1_residuals/  ← Phase 1
            dct_markers/       ← Phase 1.5
            networks/          ← Phase 2
            phase3_embeddings/ ← Phase 3
            phase3_rewiring/   ← Phase 3
            ...                ← Phase 5-9
            run_metadata.json
    """
```

### Block 3 — function docstring, lines 170-175 — guardrail: false

**Replacement one-liner:**

```python
    """Build or rebuild the Ensembl->Symbol id map over all expressed genes and curated sets."""
```

**Archived original:**

```python
    """Build/rebuild the Ensembl→Symbol ID map if needed.

    Uses the full Rtech gene list (all expressed genes) to ensure pathway
    genes outside the network skeleton are still mappable.  Also auto-resolves
    curated gene symbols from config/gene_sets.yaml.
    """
```

### Block 4 — function docstring, lines 878-882 — guardrail: false

**Replacement one-liner:**

```python
    """Project RRRM-2 WGCNA modules into external OSD cohorts; permutation-test FLT vs GC."""
```

**Archived original:**

```python
    """WGCNA module-projection external validation in OSD cohorts.

    Projects RRRM-2 WGCNA module gene sets into external cohorts and tests
    FLT vs GC module score shifts via permutation.
    """
```

### Block 5 — function docstring, lines 1007-1013 — guardrail: false

**Replacement one-liner:**

```python
    """Phase 8b: stratum-specific LOO-CV with an expression-only baseline and a shuffle null."""
```

**Archived original:**

```python
    """Phase 8b: Enhanced Predictive Validation

    Three improvements over Phase 8:
      1. Stratum-specific LOO-CV (n=10 per stratum, 5 FLT + 5 GC)
      2. Expression-based classification as control baseline
      3. Permutation null distribution (1,000 shuffles) with p-values
    """
```


---

## `src/statistics/direct_coexpression_test.py`

### Block 1 — function docstring, lines 81-84 — guardrail: false

**Replacement one-liner:**

```python
    """Build the gene-score null by shuffling FLT/GC labels; returns an (n_perms, n_genes) array."""
```

**Archived original:**

```python
    """Build gene-score null distribution by shuffling FLT/GC labels.

    Returns (n_perms, n_genes) array of null gene scores.
    """
```


---

## `src/statistics/full_pipeline_permutation.py`

### Block 1 — module docstring, lines 2-13 — guardrail: true

**Replacement one-liner:**

```python
"""Full-pipeline permutation manifest for pre-registered top hits; no run unless --execute."""
```

**Archived original:**

```python
"""
Full-pipeline permutation manifest for selected top-hit genes.

This is deliberately separate from src.statistics.permutation_bootstrap. The
fast Phase 6 test calibrates edge-sum node rewiring; this driver is reserved for
the expensive statistic that matches Phase 3 node2vec/Procrustes cosine-distance
rewiring. It permutes labels, reruns the requested pipeline commands into
isolated directories, and records a manifest for auditable execution.

By default the module writes a reviewable manifest only. Use --execute after a
pre-registered top-hit file and command template have been reviewed.
"""
```


---

## `src/statistics/full_regression.py`

### Block 1 — module docstring, lines 2-32 — guardrail: true

**Replacement one-liner:**

```python
"""Phase 6: 2x2x2 factorial edge regression; signed gene stats by permutation, not Stouffer."""
```

**Archived original:**

```python
"""
Phase 6: Edge-level regression using ALL 80 samples (full factorial model).

Fits a FULL 2×2×2 factorial model per edge to properly control for all
structure (Age, Arm, Flight) and extract clean effect estimates.

Model per edge (FLT vs GC only):
  w_e ~ Age + Arm + Flight + Age:Flight + Arm:Flight + Age:Arm + Age:Arm:Flight

Key outputs:
  - Flight main effect (FLT vs GC, controlling for Age and Arm)
  - Age×Flight interaction (does flight effect differ by age?)
  - Arm×Flight interaction (does flight effect differ by arm?)
  - Age×Arm×Flight (3-way interaction)

Then aggregate incident edge t-statistics into signed gene-level statistics and
calibrate them by within-Age×Arm label permutation. Stouffer/Brown-style
combination is not used for primary inference because unsigned edge p-values
discard direction and edge dependence is substantial on a shared skeleton.

Inputs:
  - data/processed/networks/phase2/lioness_edges.npy (N × E edge weights)
  - data/processed/networks/phase2/edge_i.npy, edge_j.npy
  - data/processed/networks/phase2/phase2_genes.txt
  - data/processed/networks/phase2/lioness_samples.txt
  - data/processed/phase1_residuals/meta_phase1.tsv.gz

Outputs:
  - data/results/phase6_regression/gene_<effect>.tsv for each effect
  - data/results/phase6_regression/edge_regression_stats.tsv (optional)
"""
```

### Block 2 — function docstring, lines 197-207 — guardrail: false

**Replacement one-liner:**

```python
    """Fit OLS across all edges at once; returns per-term beta/p arrays keyed like 'p_flight'."""
```

**Archived original:**

```python
    """
    Fit OLS regression for all edges simultaneously.
    
    Args:
        W: Edge weight matrix (N_samples × E_edges)
        X: Design matrix (N_samples × P_terms)
        term_indices: Dict mapping term names to column indices in X
    
    Returns:
        Dict with keys like 'p_flight', 'p_age_flight', 'beta_flight', etc.
    """
```


---

## `src/statistics/interaction_metrics.py`

### Block 1 — module docstring, lines 2-14 — guardrail: false

**Replacement one-liner:**

```python
"""Phase 5 derived: interaction and persistence/recovery metrics from Phase 3.3 rewiring."""
```

**Archived original:**

```python
"""
Phase 5 (derived): interaction + persistence/recovery metrics from Phase 3.3 outputs.

Inputs (from Phase 3.3):
  data/results/phase3_rewiring/*_rewiring_agg.tsv

Outputs:
  data/results/phase5_derived/
    ISS_T_interaction.tsv
    LAR_interaction.tsv
    ISS_minus_LAR_YNG_persistence.tsv
    ISS_minus_LAR_OLD_persistence.tsv
"""
```


---

## `src/statistics/permutation_bootstrap.py`

### Block 1 — module docstring, lines 2-16 — guardrail: true

**Replacement one-liner:**

```python
"""Phase 6: fast edge-sum node-rewiring tests; not inference for Phase 3 cosine rewiring."""
```

**Archived original:**

```python
"""
Phase 6: fast uncertainty tests for edge-sum node rewiring.

This module works in LIONESS edge-weight space.  For each FLT-control contrast
it computes a per-edge mean difference, then a per-gene statistic equal to the
sum of absolute incident edge differences.  These p-values are therefore
edge-sum node-rewiring tests. They are not direct inference for the Phase 3
node2vec/Procrustes cosine-distance rewiring statistic.

Focused testing is limited to statistically valid pre-specified modes:
  * --candidate-genes: BH only over an external, pre-registered candidate file.
  * --hierarchical-fdr: Benjamini-Bogomolov-style two-stage FDR over configured
    gene families/pathways. Overlapping families are additionally reported with
    a BY column because dependence is not fully eliminated by this procedure.
"""
```


---

## `src/statistics/silent_shifters.py`

### Block 1 — module docstring, lines 2-15 — guardrail: true

**Replacement one-liner:**

```python
"""Phase 5: DE-aware silent-shifter tiers; missing DE errors unless --allow_missing_de."""
```

**Archived original:**

```python
"""
Phase 5: DE-aware silent shifter generation.

Definitions:
  * candidate rewired genes: top rewiring quantile, regardless of DE.
  * DE-supported genes: high rewiring with differential-expression support.
  * strict silent shifters: high rewiring and bounded/small mean-expression
    change (shrunken |log2FC| < threshold, most of the 95% CI lies inside a
    small-effect interval, and DE FDR is not significant).
  * supported strict subset: strict silent shifters with Phase 6 support.

Missing DE is an error by default. Older rewiring-only behavior can be requested
only with --allow_missing_de and is labelled exploratory.
"""
```


---

## `src/subtype_reference/reference_builder.py`

### Block 1 — module docstring, lines 2-10 — guardrail: true

**Replacement one-liner:**

```python
"""Build frozen distal-nephron signatures from reference-only inputs; gates non-evaluable."""
```

**Archived original:**

```python
"""Build frozen distal-nephron signatures without consulting flight results.

The module consumes reference-only differential-expression and whole-kidney
expression summaries.  It deliberately has no dependency on the OSD-462
phosphoproteomic code.  Final data-derived signatures are emitted only when
both the independent distal-nephron validation and whole-kidney specificity
inputs are present; otherwise the output records a non-evaluable gate rather
than silently promoting discovery-only genes.
"""
```

### Block 2 — function docstring, lines 139-146 — guardrail: false

**Replacement one-liner:**

```python
    """Split expression breadth from the stricter distal-specificity flag frozen DCT sets use."""
```

**Archived original:**

```python
    """Separate true expression breadth from the signature-specificity filter.

    ``broadly_expressed`` requires expression in the distal target and in the
    configured number of unrelated compartments. The compatibility column
    ``non_distal_specific_or_broad`` additionally captures a single unrelated
    compartment whose expression is comparable with the distal target; frozen
    DCT signature construction continues to use that stricter flag.
    """
```

### Block 3 — function docstring, lines 496-502 — guardrail: true

**Replacement one-liner:**

```python
    """GSE150338-only DCT2/CNT signature; no GSE228367 or discovery membership is consulted."""
```

**Archived original:**

```python
    """Derive a GSE150338-only DCT2/CNT signature without GSE228367 membership.

    The fine-subtype and microdissected-segment inputs are both from
    GSE150338. The whole-kidney atlas is used only as the predeclared combined
    distal-specificity/breadth exclusion. No discovery-table membership is
    consulted.
    """
```


---

## `src/v11/aldosterone_axis.py`

### Block 1 — function docstring, lines 36-43 — guardrail: false

**Replacement one-liner:**

```python
    """Direction-corrected mean of a signed panel in one cohort, with a Wilcoxon p vs 0."""
```

**Archived original:**

```python
    """Direction-corrected mean of a directional panel within one cohort.

    ``gene_stats``: Series gene -> flight-effect statistic (z-like).
    ``signs``: Series gene -> {+1, -1, 0}; only nonzero-sign genes contribute.
    Returns the signed axis effect (positive = aldosterone program up), the
    unsigned mean for transparency, a within-panel Wilcoxon p vs 0, and the
    list of contributing genes.
    """
```

### Block 2 — function docstring, lines 72-79 — guardrail: false

**Replacement one-liner:**

```python
    """Gene-label permutation null for the panel mean; two-sided and one-sided suppression p."""
```

**Archived original:**

```python
    """Gene-label permutation null for the direction-corrected panel mean.

    Draws ``n_perm`` random gene sets of the same size from the cohort's gene
    universe, assigns them the observed panel sign vector, and recomputes the
    direction-corrected mean. Returns a two-sided competitive p and a one-sided
    *suppression* p (observed <= null), the relevant tail for the manuscript's
    distal-nephron suppression prediction.
    """
```

### Block 3 — function docstring, lines 100-105 — guardrail: false

**Replacement one-liner:**

```python
    """Pool per-cohort axis effects into a Stouffer Z, recurrence class, and sign test."""
```

**Archived original:**

```python
    """Pool per-cohort axis effects into a cross-cohort verdict.

    Reuses ``signed_stouffer_z`` for a combined Z/p and ``recurrence_class`` for
    the up/down/mixed label, and adds a directional sign-test (how many cohorts
    are negative, i.e. consistent with predicted suppression).
    """
```


---

## `src/v11/cmap_screen.py`

### Block 1 — function docstring, lines 49-55 — guardrail: false

**Replacement one-liner:**

```python
    """Load an optional mouse->human symbol map; an absent file falls back to uppercase matching."""
```

**Archived original:**

```python
    """Load an optional mouse->human symbol map.

    Expected columns are flexible but should contain mouse and human symbol
    names.  If no path is supplied or the file is absent, an empty map is
    returned and the query builder falls back to conservative uppercase symbol
    matching.
    """
```

### Block 2 — function docstring, lines 159-163 — guardrail: false

**Replacement one-liner:**

```python
    """Approximate signed CMap score; positive mimics the flight signature, negative reverses."""
```

**Archived original:**

```python
    """Approximate signed CMap score for signature rows.

    Positive means the signature mimics the mouse flight meta-signature
    (up-query genes high and down-query genes low); negative means reversal.
    """
```


---

## `src/v11/core_analysis.py`

### Block 1 — function docstring, lines 863-869 — guardrail: false

**Replacement one-liner:**

```python
    """Pick the most responsive phosphosite per parent gene (p value, then most negative effect)."""
```

**Archived original:**

```python
    """Select one phosphosite row per parent gene for row-dependence sensitivity.

    The representative row is the most statistically responsive site on the
    parent gene, using phosphosite p value as the primary key and more negative
    flight effect as the tie-breaker. This is distinct from
    ``is_single_site``, which only excludes composite/multi-position rows.
    """
```

### Block 2 — function docstring, lines 1157-1162 — guardrail: false

**Replacement one-liner:**

```python
    """Permute parent-gene subtype-prior flags within site-count strata, holding rows fixed."""
```

**Archived original:**

```python
    """Shuffle subtype-prior flags among parent genes within site-count strata.

    The row-level suppressed/not-suppressed pattern is held fixed; only the
    parent-gene subtype-prior label is permuted within bins of quantified site
    count. This preserves parent-gene site density in the null.
    """
```

### Block 3 — function docstring, lines 1510-1515 — guardrail: true

**Replacement one-liner:**

```python
    """Label a DCT1-vs-effect Spearman by observed sign; positive rho never reads as suppression."""
```

**Archived original:**

```python
    """Sign-aware interpretation of a DCT1-score vs phosphosite-effect Spearman.

    A continuous suppression gradient predicts NEGATIVE rho (more negative effect
    at higher DCT1 prior). The label reports the OBSERVED sign so a positive rho is
    never silently read as supportive of suppression.
    """
```

### Block 4 — comment run, lines 1786-1788 — guardrail: true

**Replacement one-liner:**

```python
    # Stage 0 invalidated the mediation outcome: T53/S383 rows are co-modified. Fail closed.
```

**Archived original:**

```python
    # Stage 0 invalidated the outcome used by the historical mediation:
    # position-indexed T53 and S383 rows are co-modified phosphoforms, and no
    # isolated canonical NCC/SPAK feature qualifies. Fail closed.
```


---

## `src/v11/dct_continuous_gradient.py`

### Block 1 — function docstring, lines 49-55 — guardrail: false

**Replacement one-liner:**

```python
    """OLS slope of ``outcome`` on ``coord`` plus optional covariates, in the units of ``coord``."""
```

**Archived original:**

```python
    """OLS slope of ``outcome`` on ``coord`` (+ optional covariates).

    The slope is reported in the units of ``coord`` as passed in: standardize
    ``coord`` before calling for a per-SD slope. A planted positive trend is
    recovered as a positive slope (sign-faithful), which is what the unit test
    checks.
    """
```

### Block 2 — function docstring, lines 136-143 — guardrail: false

**Replacement one-liner:**

```python
    """Natural-cubic-spline non-linearity test: nested partial F of spline(coord) vs linear."""
```

**Archived original:**

```python
    """Natural-cubic-spline non-linearity test (nested F vs the linear contrast).

    Fits a reduced model (const + coord [+ covars]) and a full model
    (const + natural-cubic-spline(coord, df) [+ covars]). The natural cubic
    spline space contains all linear functions, so the reduced model is nested
    in the full one and a standard partial F-test is valid. A small ``f_p``
    means the gradient departs from a straight line (non-monotone or curved).
    """
```


---

## `src/v11/h2_composition_aware_phospho.py`

### Block 1 — comment run, lines 537-538 — guardrail: false

**Replacement one-liner:**

```python
        # Within-site demeaning absorbs the phosphosite baseline (fixed-effect intercept analogue).
```

**Archived original:**

```python
        # Absorb phosphosite baseline by within-site demeaning. This is the
        # scalable fixed-effect analogue of a random phosphosite intercept.
```


---

## `src/v11/human_concordance.py`

### Block 1 — function docstring, lines 348-356 — guardrail: true

**Replacement one-liner:**

```python
    """Collapse scored analytes to physiological axes; per-analyte trials inflate the sign test."""
```

**Archived original:**

```python
    """Collapse scored analytes to independent physiological axes.

    The scored machine-readable + figure analytes are not statistically
    independent: several index the same physiology (e.g. 24 h urine volume and
    AQP2 both report water balance).  Treating each as its own Bernoulli trial
    inflates the sign test, so we collapse scored rows to their axis and treat
    each axis as a single trial.  An axis is concordant only if *every* scored
    analyte on it is concordant.
    """
```

### Block 2 — function docstring, lines 429-434 — guardrail: true

**Replacement one-liner:**

```python
    """Long-format OSD-656 urine inflammation panel; not inflight kidney proteomics."""
```

**Archived original:**

```python
    """Return long-format OSD-656 urine inflammation panel data.

    The submitted workbook is an Inspiration4 urine Multiplex/NULISAseq-style
    inflammation panel with preflight and recovery samples.  It is not treated as
    inflight kidney proteomics.
    """
```

### Block 3 — function docstring, lines 552-556 — guardrail: true

**Replacement one-liner:**

```python
    """Catalog OSD-656 files as recovery context only; never enters the Twins sign test."""
```

**Archived original:**

```python
    """Catalog optional OSD-656 files and summarize detectable pre/recovery markers.

    The submitted result workbook is summarized as recovery/inflammation context
    only.  It never enters the primary Twins concordance sign test.
    """
```


---

## `src/v11/kinome_atlas_ksea.py`

### Block 1 — function docstring, lines 91-97 — guardrail: false

**Replacement one-liner:**

```python
    """Read the Ser/Thr kinase atlas into per-position log2 weights keyed {position}{residue}."""
```

**Archived original:**

```python
    """Read the Ser/Thr atlas into log2 weights.

    Column headers are ``{position}{residue}`` (e.g. ``-3R``); the
    position-normalised scaled matrix is positive everywhere, so ``log2`` of each
    cell is a clean per-position log-odds contribution that sums across the
    window.
    """
```

### Block 2 — function docstring, lines 163-168 — guardrail: false

**Replacement one-liner:**

```python
    """Summed log2 PSSM score per site x kinase over 13-mers; absent residues contribute zero."""
```

**Archived original:**

```python
    """Per-site, per-kinase summed log2 PSSM score over the atlas frame.

    ``motifs`` is a sequence of 13-mers centred on the phosphosite. Returns an
    ``(n_sites, n_kinases)`` array. Residues absent from a kinase's matrix (rare
    priming positions, terminal padding) contribute nothing for that position.
    """
```

### Block 3 — function docstring, lines 202-207 — guardrail: true

**Replacement one-liner:**

```python
    """Within-cohort per-kinase percentile of each site's score, standing in for Ochoa."""
```

**Archived original:**

```python
    """Within-cohort percentile of each site's score, computed per kinase.

    Column ``k`` is replaced by the percentile rank (0--100) of each site within
    the distribution of kinase ``k``'s scores across all sites -- the on-disk
    stand-in for Johnson's Ochoa reference distribution.
    """
```

### Block 4 — function docstring, lines 226-236 — guardrail: false

**Replacement one-liner:**

```python
    """Site x kinase mask: clears the kinase percentile floor and ranks in the site's top-k."""
```

**Archived original:**

```python
    """Boolean (n_sites x n_kinases) assignment mask.

    A site is assigned to a kinase when (a) its score clears the kinase's
    ``threshold`` percentile *and* (b) that kinase is among the site's
    ``k_per_site`` best-matching kinases. The per-site top-k constraint is what
    restores substrate specificity: because the percentile is computed *within
    cohort* (no external reference), a bare per-kinase floor would assign a fixed
    ``100 - threshold`` percent of sites to *every* kinase, collapsing the KSEA
    substrate sets to a near-constant size and an n-inflated z. Capping each site
    to its top-k kinases concentrates each kinase's net on its cognate motifs.
    """
```

### Block 5 — function docstring, lines 258-262 — guardrail: false

**Replacement one-liner:**

```python
    """Assign sites to top-ranked kinases in the kinase/substrate_gene/substrate_site schema."""
```

**Archived original:**

```python
    """Assign sites to their top-ranked kinases; emit the KSEA net schema.

    Emits exactly the ``kinase / substrate_gene / substrate_site`` schema that
    :func:`load_kinase_substrate_net` validates and :func:`ksea` consumes.
    """
```

### Block 6 — function docstring, lines 300-305 — guardrail: false

**Replacement one-liner:**

```python
    """One-sided Fisher enrichment of down sites among a kinase's substrate genes, one per gene."""
```

**Archived original:**

```python
    """One-sided Fisher: are a kinase's substrate *genes* enriched for down sites?

    Sites are collapsed to one representative per parent gene (lowest phospho
    p-value, then most-negative effect) so a single gene with many sites cannot
    dominate. ``down`` := representative effect < 0 and p < ``p_thresh``.
    """
```

### Block 7 — function docstring, lines 380-385 — guardrail: false

**Replacement one-liner:**

```python
    """Join verified phospho effects to atlas motifs, keeping Ser/Thr-centred single sites only."""
```

**Archived original:**

```python
    """Join the verified phospho effect table to atlas motifs; keep S/T sites.

    Returns ``gene_symbol, site_position, motif, phospho_effect, phospho_se,
    phospho_p_value`` for every quantified single phosphosite whose motif is
    centred on Ser or Thr (Tyr sites are dropped -- the atlas is Ser/Thr only).
    """
```


---

## `src/v11/observability_audit.py`

### Block 1 — function docstring, lines 30-35 — guardrail: false

**Replacement one-liner:**

```python
    """Peptide-weighted collapse of per-protein observability to one row per gene symbol."""
```

**Archived original:**

```python
    """Peptide-weighted per-gene collapse of per-protein observability.

    Inputs come from ``osd462_anchor/protein_effects_by_row.tsv`` which is
    one row per protein with ``n_channels_used`` (count of finite scaled
    S/N channels) and ``n_peptides``.  Output is one row per gene_symbol.
    """
```

### Block 2 — function docstring, lines 92-100 — guardrail: false

**Replacement one-liner:**

```python
    """Joint abundance x peptide x missingness stratum; falls back to 5x4 if missingness is flat."""
```

**Archived original:**

```python
    """Extended (abundance × peptide × missing-fraction) joint stratum label.

    Falls back to the standard 5×4 strata if every row has the same
    ``missing_fraction`` (most TMT 2-plex scaled-S/N pools).  In that case
    the extra dimension is uninformative and the audit's Module-2 vs
    Module-3 q-value comparison should show no difference — itself a
    reportable finding ("the proteome has effectively zero missingness;
    detectability bias cannot explain the mismatch").
    """
```

### Block 3 — function docstring, lines 129-137 — guardrail: false

**Replacement one-liner:**

```python
    """Fraction of genes also protein-quantified per RNA-effect-magnitude decile."""
```

**Archived original:**

```python
    """Per-RNA-effect-magnitude-decile fraction of genes also protein-quantified.

    ``rna_table`` is the RNA universe (every gene with a finite RNA effect);
    ``protein_pool`` is the subset additionally protein-quantified.

    A monotone decline across bins would suggest large-RNA-effect genes
    are systematically detection-limited at the protein level (a real
    confounder).  A flat profile rules that out at the RNA-effect level.
    """
```

### Block 4 — function docstring, lines 163-169 — guardrail: false

**Replacement one-liner:**

```python
    """Restrict to high-coverage rows (>= 3 peptides, <= 0.2 missing fraction by default)."""
```

**Archived original:**

```python
    """Restrict to high-confidence quantification.

    Defaults (``min_peptides=3``, ``max_missing_fraction=0.2``) mirror the
    informal "high-coverage subset" used by reviewer-prep audits in the
    proteomics literature; both can be loosened or tightened for a
    sensitivity sweep.
    """
```

### Block 5 — function docstring, lines 193-198 — guardrail: false

**Replacement one-liner:**

```python
    """Matched-null propagation per pathway drawing strata from an arbitrary ``pool`` column."""
```

**Archived original:**

```python
    """Per-pathway matched-null propagation test using ``pool[stratum_col]``.

    Mirrors :func:`src.v11.rna_protein_propagation.compute_propagation_per_pathway`
    but draws strata from an arbitrary column on ``pool`` (this is how
    Module 3 swaps in the observability-extended stratum).
    """
```

### Block 6 — function docstring, lines 263-274 — guardrail: true

**Replacement one-liner:**

```python
    """NCC/SPAK site observability audit; never upgrades position-only rows to site evidence."""
```

**Archived original:**

```python
    """Per-NCC/SPAK site: observability metrics + percentile vs the phosphoproteome.

    For each residue-indexed co-modified context feature (and the
    non-regulatory sentinels), report:
      - n_fl + n_gc (channels with a finite, positive scaled S/N),
      - missing_fraction percentile within the full phospho table,
      - intensity (effect magnitude) percentile,
      - an explicit role tag that does not imply isolated canonical occupancy.

    This is an observability audit only. It cannot upgrade position-only effect
    rows into isolated canonical-site evidence.
    """
```


---

## `src/v11/publication_figures.py`

### Block 1 — function docstring, lines 213-219 — guardrail: false

**Replacement one-liner:**

```python
    """Cross-cohort recurrence figure from canonical TSV artifacts, with explicit stat labels."""
```

**Archived original:**

```python
    """Clean replacement for the legacy cross-cohort figure.

    The older static Figure 1 mixed leave-one-pathway-out diagnostics with the
    headline OSD-513 cosine, which left a stale label in the rendered panel.
    This version uses the canonical TSV artifacts and labels each statistic
    explicitly.
    """
```


---

## `src/v11/recurrence_meta.py`

### Block 1 — function docstring, lines 53-58 — guardrail: true

**Replacement one-liner:**

```python
    """Map VST columns to flight/ground; vivarium, basal and OSD-253 GCrerun are excluded."""
```

**Archived original:**

```python
    """Map VST sample columns to ``flight`` / ``ground`` from GeneLab names.

    Flight := ``_FLT_``; ground := hardware ground control ``_GC_``. Vivarium,
    basal and the OSD-253 ``GCrerun`` batch are excluded (returned for neither
    arm) so each cohort contributes one clean hardware-ground contrast.
    """
```

### Block 2 — function docstring, lines 80-90 — guardrail: false

**Replacement one-liner:**

```python
    """Per-gene flight-minus-ground VST effect, Welch SE/df; ``moderate`` adds a variance floor."""
```

**Archived original:**

```python
    """Per-gene flight-vs-ground effect and SE on the VST scale (Welch).

    ``vst`` is genes x samples (index = gene id). ``design`` maps a subset of
    columns to ``flight`` / ``ground``. The effect is ``mean(flight) -
    mean(ground)`` (positive = up in flight; sign-faithful), with Welch SE
    ``sqrt(var_f/n_f + var_g/n_g)`` and a Welch-Satterthwaite df.

    ``moderate=True`` adds an empirical-Bayes variance floor (the 10th-percentile
    sampling variance across genes) to stabilise low-variance genes, a light
    stand-in for limma's variance shrinkage.
    """
```

### Block 3 — function docstring, lines 143-147 — guardrail: false

**Replacement one-liner:**

```python
    """DerSimonian-Laird random-effects pool for one gene: effect, SE, tau^2, Q, I^2 and Wald p."""
```

**Archived original:**

```python
    """DerSimonian-Laird random-effects pool for one gene.

    ``y`` = per-cohort effects, ``v`` = per-cohort variances (se^2). Returns the
    pooled effect, its SE, tau^2, Cochran's Q, I^2 and a two-sided Wald p.
    """
```

### Block 4 — function docstring, lines 190-196 — guardrail: false

**Replacement one-liner:**

```python
    """Per-gene DerSimonian-Laird meta across cohorts with BH-FDR and a Stouffer-Z cross-check."""
```

**Archived original:**

```python
    """Per-gene DerSimonian-Laird meta across cohorts, with BH-FDR.

    ``effects`` maps cohort -> per-gene table (output of :func:`per_cohort_effect`,
    index = gene id, columns include ``effect`` and ``se``). Genes present in at
    least ``min_cohorts`` cohorts are pooled. Adds a BH-FDR over the pooled
    p-values and a Stouffer's-Z cross-check column.
    """
```

### Block 5 — function docstring, lines 230-235 — guardrail: false

**Replacement one-liner:**

```python
    """Precision-weighted gene-set score over the meta table, with Stouffer Z, I^2 and FDR."""
```

**Archived original:**

```python
    """Aggregate the pooled per-gene statistic over a curated gene set.

    Maps set symbols (case-insensitive) to Ensembl ids, intersects with the meta
    table, and returns a precision-weighted set effect, a Stouffer combination of
    the per-gene meta z, median I^2 and median per-gene FDR.
    """
```

### Block 6 — function docstring, lines 279-283 — guardrail: false

**Replacement one-liner:**

```python
    """Re-pool each gene set dropping one cohort at a time, plus the ``__none__`` baseline."""
```

**Archived original:**

```python
    """Re-estimate each gene set's pooled score dropping one cohort at a time.

    Returns one row per (dropped_cohort, gene_set) plus the ``__none__`` full-set
    baseline, so stability of the set-level effect/FDR/I^2 is auditable.
    """
```


---

## `src/v13/compartment_adversarial_audit.py`

### Block 1 — module docstring, lines 1-13 — guardrail: true

**Replacement one-liner:**

```python
"""Post-hoc OSD-462 compartment audit of parent-protein annotation enrichment; no cell of origin."""
```

**Archived original:**

```python
"""Post-hoc adversarial closure audit for OSD-462 kidney compartments.

This module deliberately separates reference-only set construction from
effect-aware post-processing.  ``prepare`` may be run before the exact label
permutation because it reads only the external kidney atlas and the frozen
KEGG structural-control source.  ``postprocess`` consumes the emitted exact
run and writes the artifact, contributor, observability, Grey60, clinical-axis,
and decision summaries.

The output object is parent-protein annotation enrichment in whole kidney.
Nothing here identifies a phosphosite's cell of origin or repairs the perfect
condition-to-reporter-block alias in OSD-462.
"""
```


---

## `src/v13/continuous_phospho_inference.py`

### Block 1 — module docstring, lines 2-13 — guardrail: true

**Replacement one-liner:**

```python
"""Exact balanced-label inference for a continuous parent-gene phospho statistic; no p cutoff."""
```

**Archived original:**

```python
"""Continuous parent-gene phosphoproteomic inference.

This module implements the prospective analysis frozen in
``config/dct_asdn_phospho_reanalysis.yaml``.  Its primary endpoint is a
continuous, parent-gene-level competitive gene-set statistic.  It deliberately
does not use a nominal phosphosite p-value threshold to define membership.

The sharp no-flight null is generated by enumerating balanced FL/GC label
assignments independently within each isobaric plex.  Every assignment is
propagated through site contrasts, parent-gene aggregation, gene-specific null
calibration, and set-level inference.
"""
```

### Block 2 — function docstring, lines 363-368 — guardrail: true

**Replacement one-liner:**

```python
    """Parse single or composite localization scores as the minimum component score."""
```

**Archived original:**

```python
    """Parse a single or composite localization score conservatively.

    Composite workbook rows store one score per component as a semicolon-
    delimited string.  Requiring the minimum component score to pass prevents
    a well-localized component from masking a poorly localized one.
    """
```

### Block 3 — function docstring, lines 510-515 — guardrail: false

**Replacement one-liner:**

```python
    """Read the workbook's unscaled summed S/N channels, separate from the scaled-channel API."""
```

**Archived original:**

```python
    """Parse the workbook's unscaled summed signal-to-noise channels.

    The shared historical parser intentionally exposes only official scaled
    channels.  This local reader keeps the new normalization sensitivity
    isolated from that established API.
    """
```

### Block 4 — function docstring, lines 609-614 — guardrail: true

**Replacement one-liner:**

```python
    """Gene-level protein log2 values aligned to phosphosite samples; missing parents stay NA."""
```

**Archived original:**

```python
    """Return gene-level protein log2 values aligned to phosphosite samples.

    Multiple protein rows for a symbol are collapsed independently in each
    sample using quantified-peptide weights.  Missing parent proteins remain
    missing and therefore cannot enter the protein-subtracted sensitivity.
    """
```

### Block 5 — function docstring, lines 704-710 — guardrail: true

**Replacement one-liner:**

```python
    """Audit scaled vs summed S/N: a constant within-row/plex ratio makes them redundant."""
```

**Archived original:**

```python
    """Audit whether official scaled and summed S/N encode the same contrast.

    GeneLab's official scaled values are row/plex rescalings of summed S/N
    values.  A constant within-row/plex ratio cancels from an uncentered log2
    FL-GC contrast, making the two inputs algebraically redundant rather than
    independent robustness evidence.
    """
```

### Block 6 — comment run, lines 990-993 — guardrail: false

**Replacement one-liner:**

```python
    # Exact fast path for the common no-duplicate-key case; the de-duplication rule is unchanged.
```

**Archived original:**

```python
    # The production OSD-462 single-site universes currently contain no
    # duplicate canonical keys. Avoid thousands of one-row pandas reductions
    # in that common case; this is an exact fast path, not a change to the
    # de-duplication rule.
```

### Block 7 — function docstring, lines 1189-1195 — guardrail: false

**Replacement one-liner:**

```python
    """Enumerate or sample balanced within-plex label assignments; row zero is the observed one."""
```

**Archived original:**

```python
    """Enumerate or sample balanced label assignments within every plex.

    The observed assignment is always row zero.  In ``exact`` mode, every
    balanced assignment is included exactly once.  In ``sampled`` mode,
    ``n_null`` assignments are sampled without replacement from the remaining
    exact assignment space.
    """
```

### Block 8 — function docstring, lines 1270-1277 — guardrail: false

**Replacement one-liner:**

```python
    """Equal-weight within-plex FL-GC site contrasts (assignments x sites); both plexes valid."""
```

**Archived original:**

```python
    """Equal-weight within-plex FL-GC contrasts for each assignment.

    The result has shape ``assignments x sites``.  Both plex-specific contrasts
    must have at least ``min_per_group`` finite observations in each permuted
    group.  The two valid plex contrasts then receive equal weight.  With
    complete balanced data this equals the flight coefficient from
    ``value ~ flight + plex``; missing data do not silently reweight one plex.
    """
```

### Block 9 — function docstring, lines 1328-1333 — guardrail: false

**Replacement one-liner:**

```python
    """Collapse signed site effects per parent gene; positive = lower phosphorylation in flight."""
```

**Archived original:**

```python
    """Collapse signed site effects to one value per parent gene.

    Positive values mean lower phosphorylation in flight.  The primary median
    remains signed.  ``one_sided_maxmean`` is a sparse-hit sensitivity and
    truncates phosphorylation increases at zero.
    """
```

### Block 10 — function docstring, lines 1502-1508 — guardrail: false

**Replacement one-liner:**

```python
    """Load the frozen long-format reference gene sets, enforcing final_for_testing and status."""
```

**Archived original:**

```python
    """Load the frozen long-format reference-builder output.

    Expected columns are ``gene_symbol`` and ``gene_set``.  When present,
    ``final_for_testing`` and ``status`` are enforced.  The literature-frozen
    ASDN list in the production config is checked against, and if necessary
    supplies, the ASDN rows.
    """
```

### Block 11 — function docstring, lines 1581-1588 — guardrail: true

**Replacement one-liner:**

```python
    """Verify strict exclusion covers every Stage-0 canonical anchor; fail closed otherwise."""
```

**Archived original:**

```python
    """Verify the strict phosphoform exclusion against Stage-0 provenance.

    The strict analysis already removes the parent genes Slc12a3 and Stk39
    from both the tested set and its eligible background.  A feature-level
    removal is therefore exactly redundant when every Stage-0-qualified
    canonical anchor belongs to a strict-excluded parent gene.  We verify that
    condition explicitly and fail closed if it is ever false.
    """
```

### Block 12 — function docstring, lines 2094-2099 — guardrail: true

**Replacement one-liner:**

```python
    """Observability-matched annotation-label null; does not replace animal-label permutation."""
```

**Archived original:**

```python
    """Secondary annotation-label null matched on observability.

    This does not replace the animal-label permutation.  It asks whether a
    frozen set is unusual among genes with similar site count, intensity, and
    missingness, using the observed gene-specific Z values.
    """
```


---

## `src/v13/reporting.py`

### Block 1 — module docstring, lines 1-7 — guardrail: true

**Replacement one-liner:**

```python
"""Read-only v13 phospho reporting; recomputes no site, permutation, or claim-gate value."""
```

**Archived original:**

```python
"""Read-only reporting for v13 continuous phosphoproteomic inference outputs.

This module deliberately consumes the frozen inference artifacts instead of
recomputing any site, parent-gene, permutation, multiplicity, or claim-gate
quantity.  Its only derived values are descriptive joins, observable-member
counts for sets omitted as non-evaluable, labels, and display ordering.
"""
```


---

## `src/validation/continuous_target.py`

### Block 1 — module docstring, lines 2-9 — guardrail: true

**Replacement one-liner:**

```python
"""Fold-safe continuous kidney-injury validation; target marker genes excluded from predictors."""
```

**Archived original:**

```python
"""
Fold-safe continuous-target validation for kidney stress/injury scores.

Targets are built from independent markers such as Havcr1/KIM-1 and Lcn2/NGAL.
Those marker genes are excluded from predictor features to avoid circularity.
Network and expression baselines are evaluated inside cross-validation folds by
Pearson and Spearman correlation.
"""
```


---

## `src/validation/cross_validation.py`

### Block 1 — module docstring, lines 2-19 — guardrail: true

**Replacement one-liner:**

```python
"""Phase 8 leakage-safe CV: residualization, skeleton E and LIONESS are rebuilt inside each fold."""
```

**Archived original:**

```python
"""
Phase 8: Leakage-Safe Predictive Validation

Implements the validation framework described in Section 5 of the methodology:
  - Stratified K-fold CV (FLT vs GC, stratified within Age×Arm)
  - Fold-wise: residualization fit on train only, skeleton E built on train only,
    LIONESS computed for train and test relative to training pool
  - Sample-level features extracted via src.validation.sample_features
  - Classification with LogisticRegression (L2) and RandomForest
  - Reports per-fold accuracy, AUC, and confusion matrix

Usage:
    python -m src.validation.cross_validation \\
        --phase2_dir data/results/<run>/networks \\
        --meta data/results/<run>/phase1_residuals/meta_phase1.tsv.gz \\
        --rtech data/results/<run>/phase1_residuals/Rtech.tsv.gz \\
        --outdir data/results/<run>/phase8_validation
"""
```

### Block 2 — function docstring, lines 48-53 — guardrail: false

**Replacement one-liner:**

```python
    """Build skeleton E from training-fold data only, by Ledoit-Wolf shrinkage and top-k."""
```

**Archived original:**

```python
    """Build skeleton E using only training-fold data.

    Performs cell-standardization within Age×Arm×EnvGroup cells (training only),
    then computes Ledoit-Wolf shrinkage covariance, inverts, and keeps top-k
    neighbors per gene.  Returns edge_i, edge_j arrays.
    """
```

### Block 3 — function docstring, lines 106-116 — guardrail: false

**Replacement one-liner:**

```python
    """Raw LIONESS contributions relative to the training pool; test samples never join it."""
```

**Archived original:**

```python
    """Compute raw LIONESS correlation contributions for samples in sample_mask.

    The pooled network is computed from ALL samples indicated by sample_mask
    (the training set).  For each sample s, the leave-one-out network is
    computed by dropping s from the pool.

    For test samples, we compute their LIONESS relative to the training pool
    (the test sample influences only its own network, not the pool).

    Returns raw LIONESS contributions of shape (len(sample_mask), n_edges).
    """
```


---

## `src/validation/enhanced_cv.py`

### Block 1 — module docstring, lines 2-9 — guardrail: true

**Replacement one-liner:**

```python
"""Phase 8b leakage-safe CV: skeleton, weights, selection, scaling, PCA and fit inside each fold."""
```

**Archived original:**

```python
"""
Phase 8b: leakage-safe enhanced predictive validation.

For every fold, this module computes the network pool, skeleton, LIONESS or
alternative sample-specific weights, feature selection, scaling, PCA, and model
fit inside the fold. It preserves expression-only baselines and evaluates
multiple LIONESS pooling modes and feature sets.
"""
```


---

## `src/validation/external_replication.py`

### Block 1 — module docstring, lines 2-15 — guardrail: true

**Replacement one-liner:**

```python
"""Cohort roles: OSD-102 primary, OSD-513 sex checks, OSD-163/253 context, OSD-568 excluded."""
```

**Archived original:**

```python
"""
Protocol-guarded independent external cohort analysis.

OSD-102 is the primary LAR-Young-like replication partner. OSD-513 is secondary
and limited to sex-robustness/sex-stratification checks. OSD-163 and OSD-253 are
context-mapping cohorts for the biology-first remodeling panel; they are not
used as strict one-to-one replication cohorts for RRRM-2 gene claims. OSD-568 is
explicitly excluded from validation claims in this remediation pass.

This module does not require ComBat-seq. Each external cohort is analyzed
independently and compared against pre-registered direction, q-value, and
pathway criteria. Multi-study pooling is handled separately by
src.validation.multi_study_pool.
"""
```


---

## `src/validation/multi_study_pool.py`

### Block 1 — module docstring, lines 2-11 — guardrail: true

**Replacement one-liner:**

```python
"""Pre-registered OSD-102 + OSD-771 LAR-Young pooling; forbidden unless ComBat-seq PCA passes."""
```

**Archived original:**

```python
"""
Pre-registered OSD-102 + OSD-771 LAR-Young pooling.

This is separate from independent external replication. Pooling is restricted to
RRRM-2/OSD-771 LAR-Young and individual OSD-102 FLT/GC mice. ComBat-seq with
study as batch and treatment preserved is required before any pooled-network
claim. PCA checks must show that study separation is reduced while treatment
signal is not erased; otherwise the module recommends fixed-effects or
meta-analysis and forbids pooled network validation.
"""
```


---

## `src/validation/sample_features.py`

### Block 1 — function docstring, lines 35-47 — guardrail: false

**Replacement one-liner:**

```python
    """Node strength (sum of incident edge weights) per sample, returned as samples x genes."""
```

**Archived original:**

```python
    """
    Compute node strength (sum of incident edge weights) per sample.
    
    Args:
        lioness_z: LIONESS weights (samples x edges)
        edge_i: Source gene indices for each edge
        edge_j: Target gene indices for each edge
        n_genes: Total number of genes
        genes: Optional gene list for column names
        
    Returns:
        Node strength matrix (samples x genes)
    """
```

### Block 2 — function docstring, lines 67-80 — guardrail: false

**Replacement one-liner:**

```python
    """Aggregate within-pathway edge weights to one summary value per sample."""
```

**Archived original:**

```python
    """
    Compute summary statistics for edges within a pathway.
    
    Args:
        lioness_z: LIONESS weights (samples x edges)
        edge_i: Source gene indices
        edge_j: Target gene indices
        genes: Full gene list
        pathway_genes: Genes in the pathway of interest
        agg_func: Aggregation function ('mean', 'median', 'sum', 'std')
        
    Returns:
        Pathway summary per sample (n_samples,)
    """
```

### Block 3 — function docstring, lines 116-131 — guardrail: false

**Replacement one-liner:**

```python
    """Aggregate each silent shifter's top-k neighbor weights into samples x shifters features."""
```

**Archived original:**

```python
    """
    Compute shifter-centered connectivity scores.
    
    For each silent shifter, aggregate weights to its top neighbors.
    
    Args:
        lioness_z: LIONESS weights (samples x edges)
        edge_i: Source gene indices
        edge_j: Target gene indices
        genes: Full gene list
        shifter_genes: Silent shifter genes
        topk_neighbors: Number of top neighbors to consider
        
    Returns:
        Shifter connectivity features (samples x len(shifter_genes))
    """
```

### Block 4 — function docstring, lines 168-178 — guardrail: false

**Replacement one-liner:**

```python
    """PCA of LIONESS edge weights; returns (PC scores, explained variance ratios)."""
```

**Archived original:**

```python
    """
    Extract PC features from LIONESS edge weights.
    
    Args:
        lioness_z: LIONESS weights (samples x edges)
        n_components: Number of PCs to extract
        stable_edge_idx: Optional subset of stable edges to use
        
    Returns:
        Tuple of (PC scores, explained variance ratios)
    """
```

### Block 5 — function docstring, lines 204-219 — guardrail: false

**Replacement one-liner:**

```python
    """Assemble node-strength, pathway, shifter and PC features into ``SampleFeatures``."""
```

**Archived original:**

```python
    """
    Extract comprehensive sample-level features.
    
    Args:
        lioness_z: LIONESS weights (samples x edges)
        sample_ids: Sample identifiers
        genes: Gene list
        edge_i: Source gene indices
        edge_j: Target gene indices
        pathway_dict: Dictionary of pathway name -> gene list
        shifter_genes: Silent shifter genes
        n_pcs: Number of PCs to extract
        
    Returns:
        SampleFeatures container
    """
```


---

## `src/validation/wgcna_external_validation.py`

### Block 1 — function docstring, lines 41-44 — guardrail: false

**Replacement one-liner:**

```python
    """Mean z-scored module expression per sample, returned as modules x samples."""
```

**Archived original:**

```python
    """Compute module scores (mean z-scored expression) for each sample.

    Returns DataFrame: modules × samples.
    """
```


---

## `src/visualization/network_diagnostics.py`

### Block 1 — module docstring, lines 2-14 — guardrail: false

**Replacement one-liner:**

```python
"""Phase 2 skeleton diagnostics: edge-weight and degree distributions, hubs, graph, z stats."""
```

**Archived original:**

```python
"""
Phase 2 Skeleton Diagnostics and Visualization

Generates publication-quality figures for the network skeleton:
1. Partial correlation distribution (edge weight histogram)
2. Degree distribution (edges per gene)
3. Top hub genes table
4. Network graph visualization (top genes by degree)
5. LIONESS z-score statistics

Usage:
    python scripts/plot_skeleton_diagnostics.py
"""
```


---

## `scripts/align_metadata_to_counts.py`

### Block 1 - function docstring (`read_counts_csv`), lines 26-30 (5 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Read a GeneLab counts CSV into a DataFrame indexed by gene_id with sample columns."""
```

**Archived original:**

```python
    """
    Reads a GeneLab counts CSV shaped like:
      gene_id, sample1, sample2, ...
    Returns a DataFrame indexed by gene_id with sample columns.
    """
```

### Block 2 - function docstring (`align_metadata_to_samples`), lines 62-65 (4 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Reorder metadata rows to the canonical sample order; raise if any sample is missing."""
```

**Archived original:**

```python
    """
    Align metadata to the canonical sample ordering.
    Raises error if any samples are missing from metadata.
    """
```

### Block 3 - function docstring (`align_counts_to_samples`), lines 75-78 (4 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Reorder count-matrix columns to the canonical sample order; raise if any sample is missing."""
```

**Archived original:**

```python
    """
    Align count matrix columns to the canonical sample ordering.
    Raises error if any samples are missing from count matrix.
    """
```

### Block 4 - function docstring (`parse_design_from_sample_name`), lines 90-95 (6 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Parse EnvGroup, Arm, Age, and AnimalCode out of an RRRM2 sample name."""
```

**Archived original:**

```python
    """
    Optional: extract simple design factors from the RRRM2 sample name.
    Example: RRRM2_R-KDN_BSL_ISS-T_YNG_BY1
    
    Returns dict with keys: EnvGroup, Arm, Age, AnimalCode
    """
```


---

## `scripts/clinical_axes/run_compartment_context.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - function docstring (`add_disjoint_podocyte_variants`), lines 90-99 (10 lines) - guardrail: true - argparse_help: false

**Replacement one-liner:**

```python
    """Add barrier-depleted podocyte sensitivity sets into the same corrected max-|T| family, not as separate tests."""
```

**Preservation note:** Prohibition: the two sensitivities must stay inside the single 51-set max-|T| family. Treating either as a separately corrected test would inflate the claim. The module docstring of this file is 1 line and IS bound to argparse (description=__doc__ at line 272) - leave it alone.

**Archived original:**

```python
    """Add two conservative high-specificity podocyte sensitivity sets.

    The original atlas-defined family remains intact.  The two added sets use
    exactly the original high-specificity podocyte definition after removing
    (i) the six frozen barrier-core genes or (ii) those six plus the two
    expanded barrier markers.  Keeping the original family members and adding
    both sensitivities makes max-|T| inference conservative over the complete
    51-evaluable-set family rather than treating either sensitivity as a
    separately corrected test.
    """
```


---

## `scripts/clinical_axes/run_podocyte_gene_coherence.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - module docstring, lines 2-7 (6 lines) - guardrail: true - argparse_help: true

**Replacement one-liner:**

```python
"""Describe gene-level podocyte RNA coherence; descriptive follow-up only, not a confirmatory testing family."""
```

**Preservation note:** Prohibition: gene-level intervals/p-values are auditability artifacts and must never be read as a second confirmatory family.

**Archived original:**

```python
"""Describe gene-level coherence of the cross-mission podocyte RNA program.

This is a descriptive follow-up to the multiplicity-controlled compartment
family.  Gene-level intervals and p-values are retained for auditability but
must not be interpreted as a second confirmatory testing family.
"""
```


---

## `scripts/clinical_axes/run_podocyte_podxl_nid1_sensitivity.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - module docstring, lines 2-10 (9 lines) - guardrail: true - argparse_help: true

**Replacement one-liner:**

```python
"""Post-hoc sensitivity forcing Podxl and Nid1 into the podocyte program; reports both genes separately."""
```

**Preservation note:** Prohibition: this is post-hoc and adversarial, not prespecified. Per-gene reporting exists so a stable set-level result cannot conceal discordant individual genes.

**Archived original:**

```python
"""Force Podxl and Nid1 into the cross-mission podocyte RNA program.

This is a post-hoc, adversarial sensitivity analysis prompted by the PODXL/NID1
comparison with prior spaceflight literature.  The frozen high-specificity
podocyte set already contains Podxl; this script adds Nid1, verifies that both
genes pass the frozen CPM eligibility rule in every mission, and reruns the
animal-level blocked-label permutation.  It also reports Podxl and Nid1
separately so a stable set result cannot conceal discordant individual genes.
"""
```


---

## `scripts/clinical_axes/run_sensitivities.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - function docstring (`secondary_qc_covariate_sensitivity`), lines 190-196 (7 lines) - guardrail: true - argparse_help: false

**Replacement one-liner:**

```python
    """Recompute Hedges g on label-blind QC-residualized scores; a sensitivity, not the primary estimand."""
```

**Preservation note:** Prohibition: deliberately NOT the primary estimand, and the residualization model must stay label-blind (score ~ QC metric). This file's 1-line module docstring is argparse-bound (line 353).

**Archived original:**

```python
    """Residualize the flagged OSD-163 mapping metric, then recompute Hedges g.

    This is deliberately a sensitivity rather than the primary estimand.  The
    residualization model is label-blind (score ~ QC metric), and the resulting
    residuals enter the same standardized flight-control effect calculation as
    the other mission scores.
    """
```


---

## `scripts/clinical_axes/run_strict_podocyte_matching_audit.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - module docstring, lines 2-8 (7 lines) - guardrail: true - argparse_help: true

**Replacement one-liner:**

```python
"""Post-hoc adversarial audit: re-null the podocyte program against label-blind atlas-rule-matched candidates."""
```

**Preservation note:** Prohibition: post-hoc adversarial sensitivity; matching uses no flight labels and no flight-effect estimates.

**Archived original:**

```python
"""Strict atlas-selection-matched audit of the podocyte RNA program.

This is a post-hoc adversarial sensitivity analysis. It replaces the original
all-gene nearest-neighbour null with candidates selected by the *same frozen
high-specificity atlas rule* as the podocyte target. Matching uses no flight
labels or flight-effect estimates.
"""
```


---

## `scripts/map_top_hits.py`

### Block 1 - module docstring, lines 2-3 (2 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Map top Ensembl IDs from the Feb 19 regression results to symbols and report DCT/NCC-WNK status."""
```

**Archived original:**

```python
"""Map top Ensembl IDs from Feb 19 regression results to gene symbols,
and check DCT/NCC-WNK gene status in results."""
```


---

## `scripts/module_convergence_fisher.py`

### Block 1 - function docstring (`reframe_eigengene_contrasts`), lines 110-117 (8 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Derive per-module ISS-T vs LAR flight effects from module_trait fits (refs: Young, ISS-T, Control)."""
```

**Preservation note:** NOT a guardrail, but the archived body is the only written record of the contrast algebra (which coefficients sum to each arm x age flight effect). Relocate that algebra to a comment or docs rather than dropping it.

**Archived original:**

```python
    """Per-module ISS-T vs LAR flight-effect comparison from existing module_trait fits.

    Coding: Age reference = Young, Arm reference = ISS-T, FlightStatus reference = Control.
      ISS-T Young flight effect = Flight
      LAR   Young flight effect = Flight + Flight:ArmLAR
      ISS-T Old   flight effect = Flight + Flight:AgeOld
      LAR   Old   flight effect = Flight + Flight:AgeOld + Flight:ArmLAR + Flight:AgeOld:ArmLAR
    """
```


---

## `scripts/osd462/05_plot_dashboard.py`

### Block 1 - comment run, lines 202-204 (3 lines) - guardrail: true - argparse_help: false

**Replacement one-liner:**

```python
    # Context-only: the T53 and S383 rows are co-modified phosphoforms, not isolated canonical sites.
```

**Preservation note:** Prohibition: these features are context-only. The one-liner must keep "not isolated canonical" or the dashboard silently reads as canonical-site evidence.

**Archived original:**

```python
    # Context-only residue-indexed features.  The T53 row is a T53/Y65
    # phosphoform and the S383 row is an S382/S383 phosphoform; neither is an
    # isolated canonical-site measurement.
```


---

## `scripts/osd462/09_stage0_manuscript_reporting.py`

### Block 1 - module docstring, lines 2-7 (6 lines) - guardrail: true - argparse_help: false

**Replacement one-liner:**

```python
"""Build manuscript displays from frozen OSD-462 Stage 0 TSVs; reporting only, changes no qualification."""
```

**Preservation note:** Prohibition: this is a reporting layer and must not alter sample, feature, or assay qualification.

**Archived original:**

```python
"""Create manuscript-ready displays from frozen OSD-462 Stage 0 outputs.

This script is deliberately a reporting layer. It reads the Stage 0 TSV
artifacts, validates their internal contracts, and changes no sample, feature,
or assay qualification.
"""
```


---

## `scripts/osd462/_common.py`

### Block 1 - function docstring (`build_symbol_to_ensembl`), lines 104-110 (7 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Map lowercase symbol -> ENSMUSG, preferring the OSD-462 DE table over id_map; first hit wins on collisions."""
```

**Archived original:**

```python
    """Symbol (lowercase) -> ENSMUSG, OSD-462 DE table preferred, id_map fallback.

    The OSD-462 differential-expression table provides an authoritative
    SYMBOL<->ENSEMBL bridge for this organism/build; we prefer it and fall back
    to the project id_map for symbols it does not contain.  One-to-many
    collisions keep the first occurrence (logged by the caller).
    """
```


---

## `scripts/phase3_node_rewiring_from_deltas.py`

### Block 1 - module docstring, lines 2-12 (11 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Convert Phase 3.1 delta-z edge arrays into per-gene rewiring_abs, rewiring_signed, and degree."""
```

**Archived original:**

```python
"""
Phase 3.1: Convert *_delta_z.npy to per-gene rewiring scores

Simple, dependency-free node rewiring metrics:
  - rewiring_abs: sum of abs(delta) over incident edges
  - rewiring_signed: sum of delta over incident edges  
  - degree: number of incident edges

Usage:
    python scripts/phase3_node_rewiring_from_deltas.py
"""
```

### Block 2 - function docstring (`node_rewiring_from_delta_edges`), lines 30-35 (6 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Compute per-gene rewiring_abs, rewiring_signed, and degree over incident edges."""
```

**Archived original:**

```python
    """
    Simple, dependency-free node rewiring metrics:
      - rewiring_abs: sum of abs(delta) over incident edges
      - rewiring_signed: sum of delta over incident edges
      - degree: number of incident edges counted (each edge contributes to both endpoints)
    """
```


---

## `scripts/phase3_procrustes_rewiring.py`

### Block 1 - module docstring, lines 2-12 (11 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Phase 3.3: align per-seed embeddings by orthogonal Procrustes and score 1 - cosine rewiring."""
```

**Archived original:**

```python
"""
Phase 3.3: Procrustes alignment + cosine rewiring (multi-seed)

- Builds anchors from the 4 FLT–GC rewiring tables (lowest median rewiring_abs)
- Aligns embeddings per seed via orthogonal Procrustes
- Computes rewiring = 1 − cosine_similarity
- Reports seed mean/std + rank variance

Usage:
    python scripts/phase3_procrustes_rewiring.py
"""
```

### Block 2 - function docstring (`pick_phase2_anchors`), lines 34-39 (6 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Pick the K genes with lowest median rewiring_abs across the four FLT-GC tables as anchors."""
```

**Archived original:**

```python
    """
    Anchor recipe:
      - read the 4 FLT_minus_GC rewiring tables
      - compute per gene median rewiring_abs
      - take bottom K genes as anchors (lowest rewiring = most stable)
    """
```

### Block 3 - function docstring (`orthogonal_procrustes_align`), lines 74-78 (5 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Rotate B onto A by orthogonal Procrustes minimizing ||B R - A||_F; return (B_aligned, R)."""
```

**Archived original:**

```python
    """
    Rotate B to best match A using orthogonal Procrustes:
      find R minimizing ||B R - A||_F
      return (B_aligned, R)
    """
```


---

## `scripts/phase4_anchor_qc_report.py`

### Block 1 - module docstring, lines 2-18 (17 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Phase 4: report anchor cross-group stability and optional not-DE status into phase4_anchor_qc/."""
```

**Archived original:**

```python
"""
Phase 4 (prereg documentation): Anchor QC report.

Goal: produce a small report showing anchors are "stable" across groups
and optionally not-DE (if you provide a gene-level DE table).

Inputs:
  - anchors.txt from Phase 3.3
  - meta_phase1.tsv.gz
  - OPTIONAL: gene_DE.tsv with columns: gene, log2FC(or logFC), FDR(or adj.P.Val)

Outputs:
  data/results/phase4_anchor_qc/
    anchors_preregistered.txt
    anchor_qc.tsv
    anchor_qc_summary.json
"""
```


---

## `scripts/phase5_derive_interaction_persistence.py`

### Block 1 - module docstring, lines 2-14 (13 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Phase 5: derive arm interaction and ISS-minus-LAR persistence tables from Phase 3.3 rewiring aggregates."""
```

**Archived original:**

```python
"""
Phase 5 (derived): interaction + persistence/recovery metrics from Phase 3.3 outputs.

Inputs (from Phase 3.3):
  data/results/phase3_rewiring/*_rewiring_agg.tsv

Outputs:
  data/results/phase5_derived/
    ISS_T_interaction.tsv
    LAR_interaction.tsv
    ISS_minus_LAR_YNG_persistence.tsv
    ISS_minus_LAR_OLD_persistence.tsv
"""
```


---

## `scripts/phase7_grounding_fast.py`

### Block 1 - module docstring, lines 2-20 (19 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Phase 7: Fisher-test pre-registered gene sets against top-decile rewiring genes, plus optional cluster enrichment."""
```

**Archived original:**

```python
"""
Phase 7: Fast biological grounding (poster-friendly).

1) Test enrichment of pre-registered gene sets among high-Δ genes (top decile)
   using Fisher's exact test.

Optional:
2) Cluster a reference embedding (k-means) and test cluster enrichment of high-Δ genes.

Inputs:
  - rewiring agg table (has gene + rewiring_mean)
  - optional embedding npy for module clustering
  - optional gene mapping TSV (Ensembl -> Symbol) for readability

Outputs:
  data/results/phase7_grounding/
    gene_set_enrichment.tsv
    cluster_enrichment.tsv (if embedding provided)
"""
```


---

## `scripts/pick_anchors.py`

### Block 1 - module docstring, lines 2-11 (10 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Select minimally rewired anchor genes for Procrustes embedding alignment."""
```

**Archived original:**

```python
"""
Select anchor genes for Procrustes alignment.

Anchors are genes with minimal network rewiring between conditions,
making them suitable reference points for embedding alignment.

Usage:
    python scripts/pick_anchors.py
    python scripts/pick_anchors.py --groupA "YNG|ISS-T|FLT" --groupB "OLD|ISS-T|FLT" --k 150
"""
```

### Block 2 - function docstring (`pick_anchors`), lines 43-54 (12 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Select the k least-rewired genes between two group keys and write them to outname."""
```

**Archived original:**

```python
    """
    Select k anchor genes with minimal rewiring between two groups.
    
    Parameters
    ----------
    groupA, groupB : str
        Group keys (e.g., "YNG|ISS-T|FLT")
    k : int
        Number of anchors to select
    outname : str, optional
        Output filename (default: auto-generated)
    """
```


---

## `scripts/plot_skeleton_diagnostics.py`

### Block 1 - module docstring, lines 2-14 (13 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Plot Phase 2 skeleton diagnostics: edge-weight and degree distributions, hubs, graph, LIONESS z-stats."""
```

**Archived original:**

```python
"""
Phase 2 Skeleton Diagnostics and Visualization

Generates publication-quality figures for the network skeleton:
1. Partial correlation distribution (edge weight histogram)
2. Degree distribution (edges per gene)
3. Top hub genes table
4. Network graph visualization (top genes by degree)
5. LIONESS z-score statistics

Usage:
    python scripts/plot_skeleton_diagnostics.py
"""
```


---

## `scripts/regulator_activity/run_phenotype_anchor.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - comment run, lines 193-195 (3 lines) - guardrail: true - argparse_help: false

**Replacement one-liner:**

```python
        # Fail-closed legacy field: no isolated canonical feature qualifies, so emit NA, not the context score.
```

**Preservation note:** Prohibition: fail-closed. Legacy consumers must receive NA rather than the co-modified context score. Dropping "not the context score" reopens the exact miscall Stage 0 closed.

**Archived original:**

```python
        # Backward-compatible fail-closed field: no isolated canonical feature
        # qualifies, so legacy activity analyses receive NA rather than the
        # co-modified context score.
```


---

## `scripts/regulator_activity/run_regulator_activity.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - function docstring (`_load_rna_effects`), lines 115-119 (5 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Build a contrasts x genes matrix from per-cohort TSVs, auto-detecting id and effect columns."""
```

**Archived original:**

```python
    """Build a contrasts x genes matrix from per-cohort gene-effect tables.

    Each file must be a TSV with a gene-id column and a flight-effect column;
    the loader auto-detects common column names.
    """
```


---

## `scripts/run_full_pipeline.py`

### Block 1 - function docstring (`run_full_pipeline`), lines 29-41 (13 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Execute analysis phases 0-7, from deconvolution through leakage-safe validation."""
```

**Archived original:**

```python
    """
    Execute complete analysis pipeline.
    
    Pipeline Phases:
        0. Preprocessing & Deconvolution
        1. Global Residualization
        2. Shared Topology Construction
        3. LIONESS Sample-Specific Networks
        4. Edge-Wise Regression
        5. node2vec Embeddings & Alignment
        6. Rewiring Quantification & Statistics
        7. Leakage-Safe Validation
    """
```


---

## `scripts/run_phase1_networks.py`

### Block 1 - module docstring, lines 2-11 (10 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Run LIONESS on Phase 1 residualized expression; emit group networks and rewiring metrics."""
```

**Archived original:**

```python
"""
Phase 1 Network Analysis Pipeline

Runs LIONESS on Phase 1 residualized expression data and computes
group-level networks and rewiring metrics.

Usage:
    python scripts/run_phase1_networks.py
    python scripts/run_phase1_networks.py --max_genes 1000 --compare "EnvGroup:FLT-vs-GC"
"""
```


---

## `scripts/run_phase2_pipeline.py`

### Block 1 - module docstring, lines 2-13 (12 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Run Phase 2: cell-standardized skeleton, LIONESS edge weights, and edge-wise regression."""
```

**Preservation note:** The archived docstring is the only place recording that step 3 requires rpy2/limma. That dependency note is already surfaced by --skip_regression help text at line 87, so it is safe to drop here.

**Archived original:**

```python
"""
Phase 2 Unified Pipeline: Cell-Standardized Shared Skeleton Construction

Runs all Phase 2 steps in sequence:
  1. Build skeleton E (cell-standardized partial correlation)
  2. Compute raw/rank-normalized LIONESS weights on E
  3. Edge-wise regression + predicted networks (requires rpy2/limma)

Usage:
    python scripts/run_phase2_pipeline.py
    python scripts/run_phase2_pipeline.py --max_genes 2500 --topk 80 --skip_regression
"""
```


---

## `scripts/run_phase3_pipeline.py`

### Block 1 - module docstring, lines 2-13 (12 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Run Phase 3: node rewiring from delta-z, multi-seed node2vec embeddings, and Procrustes rewiring."""
```

**Preservation note:** Archived text records that node2vec runs multi-seed with topology held fixed - the only in-file statement of that design choice. Preserve "multi-seed" in the replacement.

**Archived original:**

```python
"""
Phase 3 Unified Pipeline: Node2vec Embedding

Runs all Phase 3 steps in sequence:
  1. Node rewiring from delta-z (fast)
  2. Node2vec embeddings (multi-seed, topology fixed)
  3. Procrustes alignment + cosine rewiring

Usage:
    python scripts/run_phase3_pipeline.py
    python scripts/run_phase3_pipeline.py --num_seeds 2 --num_walks 20  # quick test
"""
```


---

## `scripts/run_phase4_7_pipeline.py`

### Block 1 - module docstring, lines 2-17 (16 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Run Phases 4-7 in order: derived metrics, anchor QC, uncertainty, grounding, silent shifters."""
```

**Preservation note:** Archived text records the non-obvious execution order (Phase 5 runs before Phase 4, and Phase 5 silent-shifters runs last, after Phase 6). Keep "in order" in the replacement.

**Archived original:**

```python
"""
Master runner for Phases 4-7 of the kidney transcriptome network analysis.

Executes in order:
  1. Phase 5: Derived metrics (interaction + persistence)
  2. Phase 4: Anchor QC report
  3. Phase 6: Permutation + bootstrap uncertainty
  4. Phase 7: Biological grounding
  5. Phase 5: Silent shifters (with Phase 6 support)

Usage:
  python scripts/run_phase4_7_pipeline.py [--quick]

Options:
  --quick   Run Phase 6 with minimal iterations for testing (10 perm, 10 boot)
"""
```


---

## `scripts/sf_classifier_component_decomp.py`

### Block 1 - module docstring, lines 2-8 (7 lines) - guardrail: true - argparse_help: false

**Replacement one-liner:**

```python
"""Decompose the a-priori signature into its three axes with a permutation null; no supervised fitting."""
```

**Preservation note:** Prohibition: no supervised fitting, therefore no overfit claim. Each axis is scored in its prespecified flight-predicting direction against a within-cohort label-permutation null.

**Archived original:**

```python
"""Component decomposition of the cross-cohort signature score.

Which part of the a-priori signature actually generalizes across cohorts:
the remodeling(up) axis, the DCT/NCC-WNK transport(down) axis, or the
DCT2/aldosterone(down) axis? Each scored in its flight-predicting direction,
with a within-cohort label-permutation null. No supervised fitting -> no overfit.
"""
```


---

## `scripts/sf_classifier_negctrl_figure.py`

### Block 1 - module docstring, lines 2-9 (8 lines) - guardrail: true - argparse_help: false

**Replacement one-liner:**

```python
"""Random-gene-panel feature-specificity null and summary figure; not a non-kidney tissue control."""
```

**Preservation note:** Prohibition: this is a FEATURE-specificity control only. The non-kidney TISSUE control (OSD-105 muscle) is absent from disk and remains outstanding. Dropping that caveat overclaims specificity.

**Archived original:**

```python
"""Feature-specificity negative control (random-gene score null) + summary figure.

OSD-105 muscle (a non-kidney TISSUE control) is not on disk, so the download-free
negative control here is a FEATURE-specificity control: do random size-matched
up/down gene sets reach the a-priori signature's cross-cohort AUC? If not, the
biological signature is specifically informative, not a generic property of any
gene set. (A non-kidney tissue control remains a recommended one-time add.)
"""
```


---

## `scripts/spaceflight_kidney_classifier.py`

### Block 1 - module docstring, lines 2-23 (22 lines) - guardrail: true - argparse_help: false

**Replacement one-liner:**

```python
"""LOCO-validated flight-vs-ground kidney RNA classifier with within-cohort z-scoring and specificity nulls."""
```

**Preservation note:** Prohibition set (four honesty guardrails): within-cohort z-scoring so the model cannot win on batch/strain mean shifts; LOCO so held-out batch is never trained on; three specificity controls (a-priori score, random panels, permutation null); within-cohort CV reported alongside LOCO to show generalization loss. Relocate the enumerated list to docs - a one-liner cannot carry all four.

**Archived original:**

```python
"""
Cross-cohort spaceflight kidney distal-nephron suppression classifier.

Goal (STS / v11 extension): turn the interpretive multi-omics result into a
*portable, validated artifact*. Train a flight-vs-ground detector on bulk mouse
kidney RNA and test whether it generalizes across independent cohorts (different
strains, labs, missions) with LEAVE-ONE-COHORT-OUT (LOCO) validation.

Honesty guardrails baked in:
  1. Features are z-scored WITHIN each cohort, so the model cannot win on
     cohort-level (batch/strain) mean shifts.
  2. Validation is LOCO: the held-out cohort's batch is never seen in training.
  3. Three specificity controls:
       (a) an a-priori biological signature score (no training),
       (b) size-matched RANDOM gene panels (does any panel do this well?),
       (c) within-cohort label PERMUTATION null (is the AUC above chance?).
  4. Within-cohort 5-fold CV AUC is reported alongside LOCO to show, honestly,
     how much generalization is lost across cohorts.

Cohorts (kidney, on disk): OSD-102, OSD-163, OSD-253, OSD-462, OSD-513.
Outputs -> data/results/run_20260701_sf_classifier/
"""
```


---

## `scripts/stage0/axis_effect_size_anchor.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - module docstring, lines 2-32 (31 lines) - guardrail: true - argparse_help: true

**Replacement one-liner:**

```python
"""Stage 0: pooled Hedges g for the known-positive ECM axis as a design calibration, not a biological claim."""
```

**Preservation note:** Prohibition: design-calibration estimate only. Also archived: the stop rule (if pooled g for the known positive falls below the detectable threshold, the four-axis confirmatory study is not viable) and the note that signed Stouffer tests direction, not effect size. Relocate the stop rule to docs.

**Archived original:**

```python
"""Stage 0: observed effect size for the known-positive axis.

Why this exists
---------------
A power simulation needs a reference effect size. The ECM/remodeling axis is the
one program already reported as recurrent across this corpus (signed Stouffer
p = 7.0e-4 in prior repository analyses), so it is the natural anchor: if the
random-effects pooled Hedges g for the *known positive* falls below the
detectable threshold implied by the design, the four-axis confirmatory study is
not viable and no further Stage 0 work is warranted.

Signed Stouffer tests direction and is considerably more powerful than a
random-effects meta of magnitudes. It therefore does not tell you the effect
size, which is the quantity that governs power for the axis-ranking design.

Method
------
Per cohort: per-sample mean z-score across axis genes (z computed within cohort
on VST values), then Hedges g for flight versus ground control, then a
DerSimonian-Laird random-effects pool across cohorts with I-squared.

Boundary
--------
This is a design-calibration estimate, not a biological claim. Cohorts failing
the Stage 0B coverage-confounding gate are reported but flagged, and the pooled
estimate is given both with and without them.

Usage
-----
    python3 scripts/stage0/axis_effect_size_anchor.py
"""
```


---

## `scripts/stage0/coverage_confounding_gate.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - module docstring, lines 2-27 (26 lines) - guardrail: true - argparse_help: true

**Replacement one-liner:**

```python
"""Stage 0B gate: test whether 5'/3' coverage degradation differs between flight and control per cohort."""
```

**Preservation note:** Prohibition: the eligibility question is whether degradation DIFFERS between contrasted groups - degradation is universal - and no downstream modelling repairs a cohort that fails. Keep "differs between flight and control"; a generic "assess coverage quality" would invert the gate.

**Archived original:**

```python
"""Stage 0B: does transcript coverage differ between flight and control?

The fatal-confounding gate. Lai Polo et al. (iScience 2020) showed that
carcass preservation degrades 5' gene-body coverage while leaving RIN high, and
that this distortion can exceed the biological flight effect and completely
change the differential-expression landscape.

The question that decides cohort eligibility is therefore not "is coverage
degraded" -- it is degraded everywhere -- but **whether degradation differs
between the groups being contrasted**. If flight and control samples within a
cohort have systematically different 5'/3' coverage, that cohort's flight
estimate is confounded at the measurement layer and no downstream modelling
repairs it.

Metric: ratio of RSeQC mean gene-body coverage in the 5-20% bin to the 80-95%
bin, taken from GeneLab's published ``*_qc_metrics_*.csv``. A ratio near 1 is
uniform; below 1 indicates 5' loss.

Test: Welch t on the ratio between flight and ground-control samples, plus
Hedges g, per cohort. Reported alongside the same contrast for RIN so the two
can be compared -- RIN is expected to look fine.

Usage
-----
    python3 scripts/stage0/coverage_confounding_gate.py
"""
```


---

## `scripts/stage0/protocol_inventory.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - module docstring, lines 2-42 (41 lines) - guardrail: true - argparse_help: true

**Replacement one-liner:**

```python
"""Stage 0A/0C: inventory cohort protocol and integrity metrics; strata under three cohorts block the study."""
```

**Preservation note:** Prohibition (frozen before inspection): a cohort is ineligible if flight and control differ systematically in collection, preservation, or library prep; and if no preservation x library stratum holds at least three cohorts, the confirmatory four-axis study does not proceed. This stop rule must be relocated to docs/ or eligibility_summary.md, not deleted - the replacement only points at it.

**Archived original:**

```python
"""Stage 0A/0C: protocol and transcript-integrity inventory for the OSDR kidney corpus.

Purpose
-------
Before any biological axis is tested, establish for each cohort:

* how the animal was euthanised and where (on-orbit vs post-return);
* whether tissue was dissected immediately or from an intact frozen carcass;
* which library-preparation chemistry was used (polyA vs ribodepletion/total);
* GeneLab processing pipeline version;
* per-sample RNA integrity and, critically, gene-body coverage.

Motivation is Lai Polo et al., iScience 2020 (doi:10.1016/j.isci.2020.101733),
which showed that preservation method can exceed the flight effect in magnitude,
that the distortion is amplified by polyA selection, and that RIN does **not**
predict it -- 5' gene-body coverage does. Choi et al. 2016 (PLOS One
doi:10.1371/journal.pone.0167391) separately found kidney RNA *quality* to be
minimally affected by carcass freezing to 6-7 months, so the open question is
coverage-level distortion, not RIN.

GeneLab's consensus pipeline emits RSeQC gene-body coverage into the published
``*_qc_metrics_*.csv``, so the key metric is downloadable rather than requiring
realignment.

Outputs
-------
``cohort_protocol_inventory.tsv``   one row per cohort
``sample_integrity_metrics.tsv``    one row per sample where QC is available
``eligibility_summary.md``          strata and the frozen eligibility verdict

Eligibility rule (frozen before inspection)
-------------------------------------------
A cohort is *ineligible* if flight and control samples differ systematically in
collection, preservation, or library preparation. Cohorts are then grouped into
preservation x library strata. **If no stratum contains at least three cohorts,
the confirmatory four-axis study does not proceed.**

Usage
-----
    python3 scripts/stage0/protocol_inventory.py
"""
```


---

## `scripts/subtype_reference/02a_extract_gse228367_raw_pseudobulk.py`

### Block 1 - module docstring, lines 2-8 (7 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
"""Rebuild GSE228367 subtype pseudobulks from raw 10x matrices, verifying barcode coverage and integer counts."""
```

**Preservation note:** Not a claim boundary, but the verification contract (complete barcode coverage, integer counts) is a real precondition - keep it in the replacement.

**Archived original:**

```python
"""Reconstruct GSE228367 subtype pseudobulks from official raw 10x matrices.

The author-separated RDS objects supply cell membership only. This helper
matches those cells to the NK1--NK3 filtered 10x H5 matrices in the official
GEO raw archive, verifies complete barcode coverage and integer counts, then
emits gene-symbol-collapsed count and detection pseudobulks for edgeR.
"""
```


---

## `scripts/v11/02_run_v11_core_analysis.py`

### Block 1 - comment run, lines 720-722 (3 lines) - guardrail: true - argparse_help: false

**Replacement one-liner:**

```python
    # Stage 0 invalidated this mediation outcome: T53/S383 are co-modified phosphoforms. Fail closed.
```

**Preservation note:** Prohibition: fail closed. No isolated canonical NCC/SPAK feature qualifies, so the historical mediation cannot run. Removing this comment removes the only in-code trace of why the branch is disabled.

**Archived original:**

```python
    # Stage 0 invalidated the outcome used by the historical mediation:
    # position-indexed T53 and S383 rows are co-modified phosphoforms, and no
    # isolated canonical NCC/SPAK feature qualifies. Fail closed.
```


---

## `scripts/v13/intensity_confound_audit.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - module docstring, lines 2-30 (29 lines) - guardrail: true - argparse_help: true

**Replacement one-liner:**

```python
"""Re-test frozen v13 gene sets against an intensity-stratified null; reads emitted artefacts, refits no models."""
```

**Preservation note:** Prohibition: the v13 exact run controls gene-specific VARIANCE but not a systematic intensity-dependent SHIFT - that gap is why this audit exists. Also asserts the script never refits phosphosite models.

**Archived original:**

```python
"""Intensity-confound audit of the v13 continuous phosphosite enrichment.

Motivation
----------
The v13 exact run standardises every parent gene against its own balanced-label
null, which controls gene-specific *variance* but not a systematic
*intensity-dependent shift* in the observed effects. This script measures that
shift directly and re-tests every frozen gene set against an
intensity-stratified competitive null.

It consumes only emitted artefacts of
``run_20260729_v13_continuous_phospho_exact_final`` -- it does not refit the
phosphosite models -- so it is cheap, deterministic, and auditable.

Outputs
-------
``intensity_decile_gradient.tsv``
    Mean gene-level Z by decile of median phosphopeptide signal, per profile.
``intensity_stratified_set_enrichment.tsv``
    Per gene set and profile: raw competitive statistic, statistic after
    within-stratum centring, and an intensity-stratified gene-label
    permutation p-value.
``manifest.json``
    Inputs, digests, parameters, seed.

Usage
-----
    venv/bin/python scripts/v13/intensity_confound_audit.py
"""
```

### Block 2 - comment run, lines 226-228 (3 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
                # Vectorised stratified sampling: the statistic is linear in the sampled sum.
```

**Archived original:**

```python
                # Vectorised stratified sampling without replacement: the
                # competitive statistic is a linear function of the sampled
                # sum, so only the per-permutation sum is required.
```


---

## `scripts/v13/layer_block_shift_comparison.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - module docstring, lines 2-42 (41 lines) - guardrail: true - argparse_help: true

**Replacement one-liner:**

```python
"""Test whether the flight-block signal dip is phospho-specific or shared with protein; handling effects only."""
```

**Preservation note:** Prohibition: this discriminates shared upstream handling effects, NOT every technical explanation - the two layers were labelled in separate reactions (tc882-883 phospho, tc884-885 protein), so a labelling-batch effect could still differ. Also archived: the paired-by-parent-protein comparison is the one a "phosphorylation changes without abundance changing" claim actually requires. Relocate both to docs.

**Archived original:**

```python
"""Is the flight-block signal dip phospho-specific, or shared with protein?

Question
--------
The OSD-462 channel profile shows a block-shaped dip across the flight
reporter channels (129N-131) in both plexes, recovering immediately at 131C.
A within-block positional slope cannot see a step of that shape, so two
explanations remain:

  (a) a block-level technical effect (sample handling, loading, labelling
      batch) that happens to align with condition; or
  (b) a genuine reduction in phosphopeptide signal in flight animals.

Discriminator
-------------
The protein workbook measures the *same animals* in the *same channel layout*
but without Fe-NTA phosphopeptide enrichment. A handling/loading effect
upstream of the phospho/protein split should appear in both layers at similar
magnitude. Phospho-specific suppression should not.

Two comparisons are made:

1. **Marginal.** Flight-block minus ground-block mean centred log2 signal, per
   layer and plex.
2. **Paired by parent protein.** For every protein quantified in both layers,
   the phosphosite flight-minus-ground effect minus the same protein's own
   flight-minus-ground effect. This removes anything shared by the two layers
   for that protein and is the comparison a "phosphorylation changes without
   abundance changing" claim actually requires.

Boundary
--------
The two layers were labelled in separate reactions (tc882-883 phospho,
tc884-885 protein), so a labelling-batch effect could in principle differ
between them. This test discriminates shared upstream handling effects, not
every possible technical explanation.

Usage
-----
    python3 scripts/v13/layer_block_shift_comparison.py
"""
```


---

## `scripts/v13/reporter_position_diagnostic.py`

*File note: passes `description=__doc__` to argparse.*

### Block 1 - module docstring, lines 2-29 (28 lines) - guardrail: true - argparse_help: true

**Replacement one-liner:**

```python
"""Bound OSD-462 tag-position effects from within-block slopes; condition is aliased with reporter block."""
```

**Preservation note:** Prohibition: condition is perfectly aliased with reporter-tag block, so no between-block comparison can separate biology from tag position. Only within-block slopes (6 estimates: 3 blocks x 2 plexes) are confound-free. Dropping the aliasing clause would make between-block reads look legitimate.

**Archived original:**

```python
"""Within-block reporter-position diagnostic for OSD-462.

Why this is the clean version
-----------------------------
Condition is perfectly aliased with reporter-tag block (baseline 126-128C,
flight 129N-131N, ground 131C-133C, identical in both plexes), so a
between-block comparison cannot separate biology from tag position.

But *within* a block the five channels hold five biologically exchangeable
animals of the same condition. Any systematic trend across channel position
within a block is therefore a pure tag-position effect with no biological
confound. Six independent estimates are available: 3 blocks x 2 plexes.

If a within-block slope exists, extrapolating it across the block boundary
bounds how much of the flight-minus-ground contrast reporter position alone
could produce.

Outputs
-------
``within_block_position_slopes.tsv``  per block x plex slope and permutation p
``position_effect_bound.tsv``         implied flight-vs-ground positional shift
``channel_profile.tsv``               mean centred log2 signal per channel
``manifest.json``

Usage
-----
    python3 scripts/v13/reporter_position_diagnostic.py
"""
```


---

## `tests/test_continuous_phospho_inference.py`

### Block 1 - comment run, lines 91-93 (3 lines) - guardrail: true - argparse_help: false

**Replacement one-liner:**

```python
    # Frozen estimator is the equal mean of per-plex effects (+1), not missingness-weighted pooled OLS.
```

**Preservation note:** Prohibition: names the wrong estimator the test exists to exclude. A generic "checks the plex estimate" comment would let a regression to pooled OLS pass review.

**Archived original:**

```python
    # Plex 1 has effect +2 with 4/4 observations.  Plex 2 has effect 0
    # with 3/3 observations.  The frozen estimator is their equal mean (+1),
    # not the missingness-weighted pooled OLS coefficient.
```


---

## `tests/test_osd462_anchor.py`

### Block 1 - function docstring (`_synthetic_table`), lines 15-18 (4 lines) - guardrail: false - argparse_help: false

**Replacement one-liner:**

```python
    """Build a TmtTable where gene g has FL = GC + effects[g] in log2 in both plexes, with BL = GC."""
```

**Archived original:**

```python
    """Build a TmtTable where gene g has FL = GC + effects[g] (log2) in both plexes.

    Scaled values are 2**(base + condition shift); BL is set equal to GC.
    """
```

---

## Appendix A - compact guardrail comments: do not touch

These are 2-line `#` blocks, below the 3-line threshold, so they are **not** in scope for shortening. They are listed only so the applying wave does not sweep them up as incidental cleanup: each one states a claim boundary and is already as short as it can be.

- `scripts/celltype/run_celltype_decomposition.py` lines 117-118

```python
                          # Fail-closed compatibility field: this is not an
                          # isolated NCC regulatory-site score.
```

- `scripts/clinical_axes/run_barrier_specificity_adjustment.py` lines 82-83

```python
    # The proxy is disjoint by construction, so adjustment cannot mechanically
    # subtract the same genes used in the outcome.
```

- `scripts/clinical_axes/run_cross_mission.py` lines 274-275

```python
    # Prespecified median-score sensitivity. It receives the same four-axis
    # family and complete blocked permutation rather than an unadjusted shortcut.
```

- `scripts/grey60/01_run_internal.py` lines 159-160

```python
    # Rank-based bins avoid duplicate-edge failures and lock approximately
    # equal-sized observability strata.
```

- `scripts/grey60/01_run_internal.py` lines 226-227

```python
    # The flight-by-arm column is undefined without Flight and is removed for
    # the reduced nuisance fit when testing the ISS-terminal flight coefficient.
```

- `scripts/integration/plot_sixrec_integration.py` lines 57-58

```python
    # Residue-indexed but co-modified context features.  T53 is carried with
    # Y65 and S383 with S382; neither is isolated canonical-site evidence.
```

- `scripts/osd462/01_protein_concordance.py` lines 86-87

```python
        # within-study (same-cohort) concordance isolates transcript-protein
        # uncoupling from cross-strain divergence
```

- `scripts/osd462/02_phospho_axis.py` lines 74-75

```python
    # Target-axis site table.  Join residue/phosphoform provenance before any
    # canonical-site claim; a position-only match is not sufficient.
```

- `scripts/osd462/02_phospho_axis.py` lines 135-136

```python
    # Residue-indexed NCC display.  The strict isolated-canonical set is empty
    # for the current workbook because T53 is carried on a T53/Y65 phosphoform.
```

- `scripts/osd462/02_phospho_axis.py` lines 168-169

```python
    # Strict activity support requires isolated canonical NCC and upstream
    # features.  Broad gate-gene and co-modified rows are reported separately.
```

- `scripts/osd462/_common.py` lines 41-42

```python
# OSDR RNA differential-expression contrast: Log2fc_(A)v(B) = A - B, so
# Space Flight - Ground Control is the flight effect.
```

- `scripts/regulator_activity/run_regulator_activity.py` lines 78-79

```python
        # Fail closed. A position match to T53 or S383 is insufficient because
        # both OSD-462 rows are co-modified phosphoforms.
```

- `scripts/subtype_reference/03_atlas_pseudobulk.py` lines 96-97

```python
    # The atlas layer is restricted to 3,000 variable genes. Full raw.X is
    # required for broad-expression exclusion and comparator construction.
```

- `scripts/v11/03_h2_composition_aware_phospho.py` lines 461-462

```python
        # Absorb phosphosite baseline by within-site demeaning. This is the
        # scalable fixed-effect analogue of a random phosphosite intercept.
```

- `scripts/v11/06_observability_audit.py` lines 106-107

```python
    # Pool with observability columns added; keep Module 2's match_stratum and
    # rebuild a Module 3 (abundance × peptide × missing-fraction) stratum.
```

- `tests/test_continuous_phospho_inference.py` lines 132-133

```python
    # Observed assignment has only 2 FL observations in plex 1, while pooled
    # totals remain 5 FL and 6 GC.  Strict within-plex 3/3 validity must fail.
```

- `tests/test_continuous_phospho_inference.py` lines 249-250

```python
    # U is estimable under the observed labels (2/2 per plex) but not under
    # every balanced assignment, where its four observed samples can split 3/1.
```

- `tests/test_osd462_anchor.py` lines 64-65

```python
    # When up- and down-movers balance, per-channel medians are equal across
    # conditions so loading normalization must not change individual effects.
```

- `tests/test_osd462_stage0.py` lines 234-235

```python
    # Exact source-derived output contract. These hashes intentionally fail if
    # workbook parsing, row order, fields, or QC semantics drift.
```

- `tests/test_subtype_reference_builder.py` lines 209-210

```python
    # ExternalOnly is deliberately absent from discovery. Its membership must
    # therefore be independent of GSE228367 discovery membership.
```

---

## Appendix B — earlier cleanup pass

The block below is the pre-existing `COMMENT_ARCHIVE.md` from an earlier
comment-removal pass. It is reproduced here so this file is the single
archive. **Its line numbers refer to the tree as it stood at that time and
are not valid against the current source.** Some files it names have since
been renamed or removed.

## `scripts/align_metadata_to_counts.py`

### Module docstring, original lines 2-15

```python
"""
Metadata Alignment Script for GLDS-674 / OSD-771
Aligns metadata to count matrices (STAR, RSEM, VST) and creates aligned datasets.

This script:
1. Reads count files (STAR unnormalized, RSEM unnormalized, VST normalized)
2. Reads ISA metadata
3. Aligns metadata to samples across all count matrices
4. Optionally parses design factors from sample names
5. Saves aligned outputs to data/processed/aligned_outputs/

Usage:
    python scripts/align_metadata_to_counts.py
"""
```

## `scripts/audit_fix_status.py`

### Module docstring, original lines 2-14

```python
"""
Self-audit: report which of the 14 remediation-guide fixes are detectable in
the current source tree and configuration. This is a static check (no pipeline
run required) and is intended to be a quick pre-submission sanity check.

A fix is reported `present` if all of its detection probes succeed, `partial`
if some probes succeed, `absent` if none. The script does NOT verify behavior
correctness; it verifies that the implementing artifacts exist with the
expected semantic markers.

Usage:
    python scripts/audit_fix_status.py
"""
```

## `scripts/audit_stability_gates.py`

### Module docstring, original lines 2-19

```python
"""Convenience auditor for the Phase 3 stability gate (Guardrail A).

Per agents_instruction.md §6.1, this script lives alongside
``run_contrast_vector_framework.py`` and is meant for engineers / collaborators
to confirm that the stability artifacts exist and that the recorded decisions
match the pre-registered thresholds before downstream phases run.

It does *not* re-run the bootstrap. It only inspects the already-written
``agc_stability_report.tsv`` and ``agc_stability_decision.json`` and reports:

* whether the files exist for the requested run
* whether the recorded gate values match the configured thresholds
* per-(arm, resolution) pass/fail summary
* whether the global ``fallback_to_external_axis_only`` flag is set
* whether ``--bypass-stability`` would be required to proceed

Exit code is non-zero if anything is missing or inconsistent.
"""
```

## `scripts/build_hybrid_reference.R`

### Comment block, original lines 2-22

```r
# =============================================================================
# build_hybrid_reference.R
#
# Build a hybrid single-cell reference atlas for MuSiC deconvolution by
# combining:
#   - TMS (Tabula Muris Senis): background kidney cell types (PT, Podocyte,
#     Endothelial, Mesangial, Fibroblast, Immune, CD)
#   - Chen et al. 2021 (GSE150338): high-resolution distal nephron types
#     (DCT1, DCT2, CNT, CTAL)
#
# Cluster-to-celltype mapping for Chen (from marker expression analysis):
#   DCT1  = clusters 0, 1, 2, 8, 12  (Slc12a3+, Pvalb high)
#   DCT2  = cluster 3                 (Slc12a3+, Pvalb low, Calb1+)
#   CNT   = cluster 6                 (Slc12a3+, Calb1++, Scnn1g+, Pvalb-)
#   CTAL  = clusters 4, 5, 13         (Slc12a1+, Umod+)
#   DROP  = clusters 7, 9, 10, 11     (contaminants: stromal, CD-PC, PT, IC)
#
# Output: data/processed/deconvolution/hybrid_ref/
#   matrix.mtx, barcodes.tsv, features.tsv
# =============================================================================
```

### Comment block, original lines 133-135

```r
# ── 4. Remove DCT and TAL_LOH from TMS ──────────────────────────────────────
# We need to identify which TMS cells are DCT or TAL_LOH.
# Use the same segment_from_ct mapping logic from deconvolution.R
```

### Comment block, original lines 228-230

```r
# Write features.tsv in same format as TMS
# Need: Ensembl ID (row name), feature_name (symbol)
# Map Ensembl back to symbols for feature_name
```

## `scripts/celltype/run_celltype_decomposition.py`

### Module docstring, original lines 2-14

```python
"""Cell-type marker-panel decomposition (memo Recommendation 4).

Bounds the bulk-RNA "is this a DCT program change or a composition/dilution
change?" ambiguity without new single-cell data:

1. Cross-cohort marker-panel flight effects (DCT identity vs DCT transport vs
   PT / TAL / CNT-CD / endothelial / stromal / macrophage), reusing the
   per-cohort gene-effect tables built for regulator Layer B.
2. OSD-462 animal-matched test: per-sample compartment scores, the flight
   effect of each, and whether the DCT-transport / immune / stromal scores
   covary with the measured NCC/SPAK regulatory phospho score (dilution check).
3. Scenario decision per memo Section 8.3.
"""
```

## `scripts/compare_two_runs.py`

### Module docstring, original lines 2-20

```python
"""
Comprehensive comparison of two pipeline runs:
  run_20260312_203319_2500g  (Mar 12 — "NEW")
  run_20260226_132416_2500g  (Feb 26 — "OLD")

Sections covered:
  0. Run metadata & configuration
  1. Phase 0 – Deconvolution (cell-type proportions, pseudobulk recovery)
  2. Phase 1 – Residualization / gene counts
  3. Phase 2 – Network skeleton (edge & gene counts)
  4. Phase 3 – Node2Vec embeddings & seed stability
  5. Phase 3 – Procrustes rewiring (aggregated scores, top-50 genes)
  6. Phase 5 – Derived interaction / persistence metrics
  7. Phase 5 – Silent shifters (strict)
  8. Phase 6 – Edge regression (gene-level flight & interaction)
  9. Phase 6 – Permutation / bootstrap uncertainty
 10. Phase 7 – Gene-set enrichment / biological grounding
 11. Figure-level summary tables
"""
```

### Comment block, original lines 52-54

```python
# ═══════════════════════════════════════════════════════════════════════════════
#  0.  RUN METADATA
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 69-71

```python
# ═══════════════════════════════════════════════════════════════════════════════
#  1.  PHASE 0 — DECONVOLUTION
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 141-143

```python
# ═══════════════════════════════════════════════════════════════════════════════
#  2.  PHASE 1 — RESIDUALIZATION / GENE COUNTS
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 177-179

```python
# ═══════════════════════════════════════════════════════════════════════════════
#  3.  PHASE 2 — NETWORK SKELETON
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 248-250

```python
# ═══════════════════════════════════════════════════════════════════════════════
#  4.  PHASE 3 — NODE2VEC EMBEDDING STABILITY
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 282-284

```python
# ═══════════════════════════════════════════════════════════════════════════════
#  5.  PHASE 3 — PROCRUSTES REWIRING
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 362-364

```python
# ═══════════════════════════════════════════════════════════════════════════════
#  6.  PHASE 5 — DERIVED INTERACTION / PERSISTENCE
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 403-405

```python
# ═══════════════════════════════════════════════════════════════════════════════
#  7.  PHASE 5 — SILENT SHIFTERS (STRICT)
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 452-454

```python
# ═══════════════════════════════════════════════════════════════════════════════
#  8.  PHASE 6 — GENE-LEVEL EDGE REGRESSION
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 531-533

```python
# ═══════════════════════════════════════════════════════════════════════════════
#  9.  PHASE 6 — PERMUTATION / BOOTSTRAP UNCERTAINTY
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 567-569

```python
# ═══════════════════════════════════════════════════════════════════════════════
# 10.  PHASE 7 — GENE-SET ENRICHMENT
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 627-629

```python
# ═══════════════════════════════════════════════════════════════════════════════
# 11.  FIGURES-LEVEL DECONVOLUTION SUMMARIES
# ═══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 643-645

```python
# ═══════════════════════════════════════════════════════════════════════════════
# 12.  OVERALL STRUCTURAL DIFF — files unique to each run
# ═══════════════════════════════════════════════════════════════════════════════
```

## `scripts/explore_chen_atlas.R`

### Comment block, original lines 1-11

```r
# Marker-based cluster annotation for Chen atlas
# Key markers from Chen et al. 2021 and kidney biology:
#   - Slc12a1 (NKCC2): CTAL marker
#   - Umod (Uromodulin): TAL marker
#   - Slc12a3 (NCC): DCT marker
#   - Pvalb: DCT1 marker
#   - Trpm6: DCT2 marker (but also CNT)
#   - Calb1: CNT marker
#   - Aqp2: CD-PC marker
#   - Slc26a4 (Pendrin): CD-IC marker
#   - Slc34a1: PT marker
```

## `scripts/extract_ecm_edges.py`

### Comment block, original lines 92-96

```python
# What does a positive/negative interaction beta mean?
# flight=1(FLT), 0(GC); arm=1(ISS-T), 0(LAR).
# arm_flight interaction = (ISS-T FLT - ISS-T GC) - (LAR FLT - LAR GC)
# Positive: The increase in connectivity due to flight is GREATER on ISS-T than LAR.
# Negative: The increase in connectivity due to flight is WEAKER (or decreases) on ISS-T compared to LAR.
```

## `scripts/integration/plot_sixrec_integration.py`

### Module docstring, original lines 2-13

```python
"""Consolidated six-recommendation integration figure.

Ties together the memo's Recommendation analyses into one cross-omic
decoupling / evidence-ladder panel:

  A  RNA pathway-vector recurrence (RRRM-2 ISS-T vs OSD-462)        [Rec 3 / 1]
  B  Protein-abundance concordance null (matrix & DCT)              [Rec 3]
  C  Phospho layer: NCC regulatory vs non-regulatory + total NCC    [Rec 3]
  D  Regulator activity: KSEA WNK/SPAK + recurrent pathways/TFs     [Rec 2]
  E  Cell-type decomposition: DCT down, interstitial up (cohorts)   [Rec 4]
  F  Evidence ladder + honest grades (Rec 1-6)                      [Rec 5/6]
"""
```

## `scripts/make_consensus_anchors.py`

### Comment block, original lines 47-51

```python
    # score: lower is better
    # - median drives stability
    # - max punishes “good usually, bad once”
    # - iqr punishes inconsistency
    # - missing penalizes absent genes
```

## `scripts/map_top_hits.py`

### Comment block, original lines 18-20

```python
# ============================================================
# TASK 1: Map top regression hits
# ============================================================
```

### Comment block, original lines 67-69

```python
# ============================================================
# TASK 2: Find DCT/NCC-WNK genes in the map
# ============================================================
```

### Comment block, original lines 95-97

```python
# ============================================================
# TASK 3: Check DCT genes in regression results
# ============================================================
```

### Comment block, original lines 123-125

```python
# ============================================================
# TASK 4: Check DCT genes in rewiring results
# ============================================================
```

### Comment block, original lines 152-154

```python
# ============================================================
# TASK 5: Check permutation pvals for DCT genes
# ============================================================
```

### Comment block, original lines 177-179

```python
# ============================================================
# TASK 6: Summary of enrichment gene set overlap
# ============================================================
```

## `scripts/module_convergence_fisher.py`

### Module docstring, original lines 1-13

```python
"""
Module convergence: Fisher's-exact overlap between
  (1) 11 ISS-T-Young LIONESS edge-rewiring candidates (q<0.05) and each WGCNA module
  (2) ISS-T-Young strict silent shifters and each WGCNA module
Plus (3) reframe of existing eigengene FLT-vs-GC contrasts as cross-method consistency.

Universe = genes that were eligible for the relevant LIONESS/silent-shifter statistic
AND assigned to a WGCNA module. This avoids inflating Fisher enrichment by counting
WGCNA genes that could never have been selected by the candidate-generating analysis.
Dropped genes are logged in the output.

Outputs land in data/results/run_20260517_213205_2500g/module_convergence/.
"""
```

## `scripts/osd253_strain_module_projection.py`

### Module docstring, original lines 1-14

```python
"""Metadata-aware OSD-253 strain-stratified WGCNA module projection.

This is the decision analysis for the RRRM-2 grey60/TLR4 question:

1. Parse ISA sample and assay metadata, not filename regexes alone.
2. Score RRRM-2 WGCNA modules in OSD-253.
3. Test FLT-minus-control module shifts within Strain x Duration strata.
4. Test whether C3H/HeJ attenuates the C57BL/6J flight effect.
5. Run original-GC and white-light rerun-control sensitivity analyses separately.

The rerun-control analysis fixes the OSD-253 light-wavelength mismatch, but it
introduces sequencing/read-length/run differences. It is therefore a sensitivity
analysis, not a clean replacement for the original GC comparison.
"""
```

## `scripts/osd462/00_harmonize.py`

### Module docstring, original lines 2-13

```python
"""Layer 0 - Harmonization.

Builds the master ``osd462_flight_effects.tsv`` joining, per gene:
  * RRRM-2 ISS-T RNA flight effect (reference; from lar_reversal_gene_scatter)
  * OSD-462 RNA flight effect (Space Flight - Ground Control; OSDR DE table)
  * OSD-462 protein flight effect (FL - GC, plex-corrected; TMT workbook)
  * peptide count, plex coverage, protein abundance, matched-null strata bins

Usage::

    python scripts/osd462/00_harmonize.py [--run RUN_NAME]
"""
```

## `scripts/osd462/01_protein_concordance.py`

### Module docstring, original lines 2-19

```python
"""Layer 1 - Protein-level concordance with abundance/peptide-matched null.

For each targeted gene set, estimate three statistics on its protein-quantified
members and test each against an abundance x peptide-count-matched random
gene-set null (10,000 draws):

  1. signed mean protein flight effect (signed by the RRRM-2 RNA direction)
  2. Spearman concordance between RRRM-2 ISS-T RNA effect and OSD-462 protein effect
  3. RNA<->protein sign-agreement rate

``tubular_transport_broad`` is carried as a specificity control (broad transport
should be less concordant than the focused DCT/NCC/WNK axis).  The genome-wide
RNA<->protein correlation is reported as background context.

Usage::

    python scripts/osd462/01_protein_concordance.py --run RUN_NAME
"""
```

## `scripts/osd462/02_phospho_axis.py`

### Module docstring, original lines 2-14

```python
"""Layer 2 - Phosphoprotein activity of the WNK-SPAK/OSR1-NCC axis (conditional).

Hard checkpoint first: confirm the phosphoproteomics workbook quantifies the
NCC (Slc12a3) regulatory phosphosites and the SPAK/OSR1 (Stk39/Oxsr1) and WNK
sites.  If the key sites are absent, the activity claim is abandoned and we
report that honestly.  If present, estimate per-site FL - GC effects
(plex-corrected linear model, with CIs), both raw and normalized to total
protein abundance (phospho-occupancy), and integrate across layers.

Usage::

    python scripts/osd462/02_phospho_axis.py --run RUN_NAME
"""
```

### Comment block, original lines 139-141

```python
    # Gate-site support must include both the NCC regulatory cluster and the
    # upstream SPAK/OSR1 arm; one unrelated gate-site is not enough to claim
    # reduced NCC-pathway activity.
```

## `scripts/osd462/03_network_translation.py`

### Module docstring, original lines 2-17

```python
"""Layer 3 - Network-candidate translation check.

The only test in the anchor that speaks to the LIONESS/node2vec network layer.
Take the RRRM-2 top composite network candidates (gene_axis_priority.tsv) and
ask whether they are enriched among OSD-462 protein- and phospho-changing genes
relative to an abundance/detectability-matched random null.  Both outcomes are
reported honestly:

  * enriched     -> the network layer nominated genes with independent
                    proteomic / phosphoproteomic support
  * not enriched -> the network layer did not translate to protein biology

Usage::

    python scripts/osd462/03_network_translation.py --run RUN_NAME
"""
```

## `scripts/osd462/04_rna_recurrence.py`

### Module docstring, original lines 2-17

```python
"""Layer 4 - Same-modality RNA recurrence gate.

Before trusting any cross-modality (protein/phospho) result, confirm that the
OSD-462 *RNA* flight effect recurs the RRRM-2 ISS-T matrix-high / DCT-low
direction.  We build a pathway-effect vector for each cohort (mean flight
effect per mechanism gene set) and compute their cosine alignment, with a
sample-resampling bootstrap CI and a leave-one-pathway-out robustness check -
mirroring the OSD-513 cross-OSDR recurrence test.

This is a gate: if OSD-462 RNA does not recur the signal, a protein-concordance
result would be hard to interpret.

Usage::

    python scripts/osd462/04_rna_recurrence.py --run RUN_NAME
"""
```

## `scripts/osd462/05_plot_dashboard.py`

### Module docstring, original lines 2-14

```python
"""Layer 5 - Multi-omics anchor figures.

Produces two figures from the anchor outputs:
  1. ``fig_osd462_multiomics_dashboard`` (2x2): RNA->protein scatter; DCT-axis
     protein-abundance bars; WNK-SPAK/OSR1-NCC phosphosite effects with CIs;
     network-candidate translation enrichment vs matched null.
  2. ``fig_osd462_rna_recurrence`` (1x2): pathway-vector concordance (RNA gate)
     and the cross-layer DCT/NCC summary (RNA / protein abundance / phospho).

Usage::

    python scripts/osd462/05_plot_dashboard.py --run RUN_NAME
"""
```

## `scripts/osd462/06_compile_summary.py`

### Module docstring, original lines 2-13

```python
"""Layer 6 - Consolidate anchor results into a single summary + run manifest.

Reads every layer output and emits:
  * ``results_summary.json`` - headline numbers + the pre-registered hypothesis
    verdicts (H1-H4) and the decision-table row the data landed on.
  * ``manifest.json``        - run-level manifest referencing each layer manifest
    with input SHAs.

Usage::

    python scripts/osd462/06_compile_summary.py --run RUN_NAME
"""
```

## `scripts/osd462/07_build_compendium.py`

### Module docstring, original lines 2-10

```python
"""Layer 7 - Build the standalone LaTeX results compendium from anchor outputs.

Reads every layer artifact and emits ``latex_paper/osd462_multiomics_compendium.tex``
with all tables populated from data (so numbers cannot drift from the analysis).

Usage::

    python scripts/osd462/07_build_compendium.py --run RUN_NAME [--compile]
"""
```

## `scripts/osd462/_common.py`

### Module docstring, original lines 1-6

```python
"""Shared paths and helpers for the OSD-462 multi-omics anchor layer scripts.

All layer scripts (``00_harmonize`` ... ``04_rna_recurrence``) import from here
so that input locations, the ID bridge, gene-set loading, and manifest writing
are defined once.
"""
```

## `scripts/plot_manuscript_v4_figures.py`

### Module docstring, original lines 2-7

```python
"""Generate manuscript v4 figures from mechanism-axis outputs.

The visual style intentionally mirrors the original pipeline publication figures:
whitegrid background, high-contrast FLT/GC colors, individual mouse points,
median bars, compact multi-panel layouts, and PDF/PNG output.
"""
```

## `scripts/plot_output.py`

### Module docstring, original lines 2-9

```python
"""
plot_output.py - Generate publication-ready plots and tables for each pipeline phase.

Organizes outputs into a clean folder structure mirroring the export script's phase organization.

Usage:
    python scripts/plot_output.py --out_root plots --tag paper_figures
"""
```

## `scripts/regulator_activity/build_rna_effects.py`

### Module docstring, original lines 2-19

```python
"""Build per-cohort gene-level flight-effect tables for regulator Layer B.

decoupler PROGENy / CollecTRI priors are keyed by *gene symbol*, so each cohort
table is emitted as ``gene`` (mgi symbol) + ``stat`` (a signed flight-effect
statistic, Space Flight - Ground Control). The per-cohort effect definition is:

* RRRM-2 ISS-T / LAR (Young): Wald-style z = log2FC / lfcSE from the project
  limma gene-level DE tables (``data/processed/gene_level_DE``).
* OSD-462 RNA: the OSDR-provided moderated ``Stat_(Space Flight)v(Ground
  Control)`` (symbol bridge from the same DE table).
* OSD-513 / OSD-253: Welch t per gene (Space Flight vs Ground Control) computed
  directly from the GeneLab VST matrix (no DE table downloaded). OSD-253 pools
  strains/durations into a single cohort-level SF-vs-GC contrast and is treated
  as a context cohort.

All effects are *within cohort*; decoupler Layer B is run per cohort and only
the cross-cohort recurrence of activity sign is interpreted (never raw pooling).
"""
```

## `scripts/regulator_activity/run_phenotype_anchor.py`

### Module docstring, original lines 2-8

```python
"""Phenotype-anchoring layer -- animal-matched RNA state vs NCC activity.

Builds a per-animal NCC/SPAK regulatory-phosphorylation activity score and a
per-animal DCT/NCC-low RNA state score for OSD-462 / RR-10, matches them by
animal, and compares at group / all-sample / condition-adjusted levels with two
controls (non-regulatory phosphosites; RNA score with Slc12a3 removed).
"""
```

### Comment block, original lines 149-153

```python
    # RNA score = mean-z of DCT/NCC-WNK transport genes (higher = more DCT
    # transport program). Signed the same biological direction as the NCC
    # activity score (higher = more NCC activating phosphorylation), so a
    # concordant DCT-suppressed flight state lowers both and gives same-sign
    # group effects and a positive correlation.
```

## `scripts/regulator_activity/run_regulator_activity.py`

### Module docstring, original lines 2-21

```python
"""v10 regulator-activity prioritization -- orchestrator.

Layer A (kinase activity, KSEA) runs fully offline and is executed here.
Layer B (TF/pathway activity) uses decoupler with PROGENy + DoRothEA/CollecTRI
priors fetched from omnipath/zenodo; that fetch needs network access at run
time. If the priors cannot be fetched, Layer B is skipped with a clear message
and the run recipe is printed -- Layer A outputs are still produced.

Usage
-----
    python3 scripts/regulator_activity/run_regulator_activity.py \
        --outdir data/results/run_YYYYMMDD_regulator_activity

Inputs (defaults point at the verified v9 / OSD-462 anchor outputs):
    --phospho-sites   OSD-462 genome-wide phosphosite flight effects
    --ks-net          kinase-substrate table (curated core ships in-repo;
                      pass a PhosphoSitePlus-format table for the full panel)
    --rna-effects     optional JSON mapping {cohort_label: path} of gene-level
                      flight-effect tables for Layer B
"""
```

### Comment block, original lines 71-74

```python
# ---------------------------------------------------------------------------
# Layer A -- KSEA kinase activity
# ---------------------------------------------------------------------------
```

### Comment block, original lines 94-97

```python
# ---------------------------------------------------------------------------
# Layer B -- TF / pathway activity (decoupler, network-dependent priors)
# ---------------------------------------------------------------------------
```

## `scripts/regulator_activity/summarize_regulator_recurrence.py`

### Module docstring, original lines 2-17

```python
"""Summarize cross-cohort recurrence of decoupler Layer-B regulator activity.

Reads the PROGENy pathway and CollecTRI TF activity tables produced by
``run_regulator_activity.py`` and asks the memo's Recommendation-2 question: is
there a *non-obvious* upstream regulator that recurs across independent flight
cohorts, beyond generic injury/inflammation/stress?

Recurrence is evaluated over the *terminal-flight* cohorts only (ISS-T Young,
OSD-462, OSD-513, OSD-253). The RRRM-2 LAR (live-animal-return) arm is reported
separately as a recovery/reversal contrast, not as a recurrence cohort.

Outputs (into the regulator run dir):
* ``progeny_recurrence.tsv`` / ``collectri_recurrence.tsv`` -- per-source
  recurrence summary, ranked.
* ``regulator_recurrence_verdict.json`` -- the honest Rec-2 verdict.
"""
```

### Comment block, original lines 38-40

```python
# Curated set of regulators that are *expected* from generic kidney
# injury/inflammation/fibrosis/stress/proliferation -- recurrence here is
# literature-consistent, not a discovery. (lowercased match)
```

## `scripts/run_contrast_vector_framework.py`

### Module docstring, original lines 2-8

```python
"""Run the Cross-OSDR Network-Contrast Framework.

This script is intentionally separate from the legacy phase runner and is
invoked from ``src/run_all_phases.py`` only when the new opt-in flag is used.
It implements the RRRM-2 stability gate and within-cohort decomposition path,
with explicit external-axis and decision artifacts.
"""
```

## `scripts/run_deconvolution.R`

### Comment block, original lines 151-156

```r
# Handle potential transpose issues with mtx (Scanpy writes shape n_obs x n_vars, readMM reads as such)
# But SCE expects genes x cells (n_vars x n_obs). Scanpy usually writes cells x genes?
# Creating scipy.io.mmwrite(..., raw_adata.X) usually writes standard MM.
# If raw_adata.X is (cells, genes), then readMM returns (rows=cells, cols=genes).
# SCE needs (rows=genes, cols=cells). So we likely need to Transpose.
# Let's check dimensions carefully.
```

### Comment block, original lines 609-611

```r
# CRITICAL: donor_id diagnostics
# MuSiC uses cross-donor variance for cell-size (theta) estimation.
# If donor_id is NA, constant, or unique per cell, resolution fails.
```

### Comment block, original lines 655-658

```r
# Option to force single donor for testing
# Uncomment the line below to bypass donor issues:
# colData(sce2)$donor_id <- "D1"
```

### Comment block, original lines 671-674

```r
# In-silico pseudo-bulk test for MuSiC
# Requires: sce2 with counts, colData(sce2)$segment_use and $clusterType and $donor_id
# Also requires: clusters.type and group.markers as in your script
```

### Comment block, original lines 677-679

```r
# ---------
# Helpers
# ---------
```

### Comment block, original lines 697-699

```r
# ---------
# IMPROVED: Pseudobulk with exact proportions and all fixes
# ---------
```

### Comment block, original lines 769-772

```r
    # FIX 3: Independent DCT variation via two-level design
    # Sample tubule proportions with independent ranges

    # Step 1: Sample coarse group totals
```

### Comment block, original lines 913-915

```r
# ---------
# Build pseudo-bulk mixtures
# ---------
```

### Comment block, original lines 948-950

```r
  # Two ground-truth matrices:
  # P_cell: fraction of sampled cells from each segment
  # P_rna:  fraction of total UMIs contributed by each segment
```

### Comment block, original lines 1016-1018

```r
# ---------
# Experiment B: Per-donor pseudo-bulk (sample from single donor per mixture)
# ---------
```

### Comment block, original lines 1131-1133

```r
# ---------
# Run pseudo-bulk test (with normalize_flag parameter)
# ---------
```

### Comment block, original lines 1349-1351

```r
# align to pseudo-bulk sample names
# MuSiC returns rows = bulk samples; our pb_mat cols are samples
# So rownames(P_hat) should match colnames(pb_mat)
```

### Comment block, original lines 1486-1492

```r
# FIX APPLIED 2026-01-11: CPM Normalization for Validation
# Problem: prop_seg is library-size independent (fractions sum to 1), but raw
#          bulk counts scale with sequencing depth. This caused spurious low/negative
#          correlations in validation (e.g., DCT showed ρ=0.03 with raw, ρ=0.33 with CPM).
# Solution: Normalize bulk to CPM before computing marker scores.

# CPM normalize bulk matrix for validation
```

### Comment block, original lines 1708-1711

```r
# 6) Use segment-direct proportions (already at segment level)
# Since we ran MuSiC on segments directly, prop_seg is already at segment granularity
# No need to aggregate cell types - the output is already segment-level
```

## `scripts/run_full_pipeline.py`

### Module docstring, original lines 2-9

```python
"""
Master Pipeline Orchestration Script for RRRM-2 Kidney Network Analysis

Runs the complete analysis pipeline from raw data to validated rewiring signatures.

Usage:
    python scripts/run_full_pipeline.py --config config/hyperparameters.yaml
"""
```

## `scripts/run_lar_reversal_analysis.py`

### Module docstring, original lines 2-7

```python
"""Run LAR reversal and mechanism-switch analysis.

This secondary layer asks whether LAR flight effects are best described as
attenuation, true reversal of the ISS-T remodeling direction, or a mechanism
switch into other clock/DCT, S1P, immune, or preservation-context programs.
"""
```

## `scripts/run_mechanism_axis_prioritization.py`

### Module docstring, original lines 2-8

```python
"""Run recurrent remodeling-axis and network-neighborhood prioritization.

This phase sits downstream of the contrast-vector recurrence result.  It
projects RRRM-2, OSD-513, and OSD-253 samples onto the recurrent ISS-T
remodeling axis, tests candidate mechanism scores against that axis, and uses
LIONESS/node2vec artifacts only to prioritize exploratory local contributors.
"""
```

## `scripts/sanity_marker_panel.R`

### Comment block, original lines 1-14

```r
# Sanity Marker Panel Script
# Purpose: Validate marker expression across nephron segments to ensure correct
#          cell type labeling in the reference SCE before deconvolution.
#
# This script:
#   1) Maps marker panel genes into the gene ID space used by sce2
#   2) Computes mean logCPM and % cells expressing for each marker by segment_use
#   3) Identifies where each marker peaks to detect potential mislabeling
#   4) Optionally refines segment labels (DCT1, DCT2/CNT, CD_PC, CD_IC) using
#      module-score winner assignment for improved deconvolution
#
# Prerequisites: Run this AFTER run_deconvolution.R has created `sce2` so that
#                `sce2` is in your global environment with segment_use assigned.
```

### Comment block, original lines 148-166

```r
# How to interpret the output (the "skeptical" read):
#
# - If Aqp2 / Avpr2 / Scnn1g peak in "Other", you likely have true collecting
#   duct principal cells hiding in "Other".
# - If Atp6v1b1 / Foxi1 / Slc4a1 peak in "Other", you likely have intercalated
#   cells hiding in "Other".
# - If your current "CD" peaks for Calb1/Trpv5/Slc8a1/Klk1, that's CNT/DCT2,
#   NOT true collecting duct.
# - If Slc12a1/Umod don't peak in your TAL_LOH, your TAL label may be off
#   (or TAL is underrepresented).


# PART B: Rename/split distal segments (practical way)
# Option 1: Simple "module-score winner" relabel inside (DCT, CD, Other)
#
# This is the fastest approach, and it's surprisingly effective.
# After running this, use `segment_refined` instead of `segment_use` in MuSiC.

# Set to TRUE to enable the relabeling (run Part A first to inspect markers)
```

### Comment block, original lines 221-230

```r
    # What this accomplishes:
    #   - If true collecting duct is hiding in "Other", it becomes CD_PC / CD_IC.
    #   - Your current "CD" that's really CNT becomes DCT2CNT.
    #   - Your "DCT" that's really DCT1 stays DCT1.
    #   - You can keep TAL separate or ignore it.
    #
    # IMPORTANT: Tell MuSiC to use `segment_refined` instead of `segment_use`
    #   Wherever your deconvolution script currently uses `segment_use`,
    #   switch the cluster column to `segment_refined`.
```

## `scripts/v11/02_run_v11_core_analysis.py`

### Module docstring, original lines 2-8

```python
"""Execute the core v11 DCT1/phosphoproteome/mediation analyses.

This script intentionally keeps the external references as priors/scaffolds:
GSE228367 is a DCT subtype RNA reference, PXD001729 is a cultured DCT-lineage
phosphoproteomic reference, and OSD-462 remains the only spaceflight
phosphoproteomic anchor.
"""
```

## `scripts/v11/03_h2_composition_aware_phospho.py`

### Module docstring, original lines 2-13

```python
"""Composition-aware H2 robustness tests for OSD-462 phosphoproteomics.

This script tests whether the v11 DCT1-prior phosphosite signal survives a
conservative adjustment ladder using per-animal TMT phosphosite intensities.

Claim discipline:
  * This is not DCT-specific phosphoproteomics.
  * Composition scores are bulk RNA marker-panel estimates and may sit on the
    biological path from flight to phosphosite dilution.
  * The output should be described as composition-aware robustness, not
    deconvolution.
"""
```

## `scripts/v11/04_spatial_reference_projection.py`

### Module docstring, original lines 2-15

```python
"""External kidney injury/repair spatial reference projection for v11.

Primary spatial source:
  * GSE269622 Visium, whole-transcriptome Space Ranger archives.

Secondary spatial source:
  * GSE269719 processed Xenium AnnData, used only for annotation/neighborhood
    inventory because the panel is targeted and not a genome-wide projection
    source.

This script does not localize the RR-10 spaceflight lesion. It asks which IRI
timepoints or marker-enriched Visium spatial niches most resemble the bulk
spaceflight RNA remodeling vector.
"""
```

## `scripts/v11/06_observability_audit.py`

### Module docstring, original lines 2-28

```python
"""v11 Module 3 — proteome observability-bias audit.

Pre-empts the first-line reviewer objection to the v11 RNA→protein
discordance: "the mismatch is just proteome detectability."

Pipeline:

  1. Build per-gene observability features from
     ``protein_effects_by_row.tsv`` (peptide-weighted ``n_channels_used``,
     ``missing_fraction``).
  2. Detectability gradient — fraction of RNA-detected genes also
     protein-quantified, per RNA-effect-magnitude decile.
  3. Observability-matched re-test — extend the matched-null stratum with
     a missing-fraction bin; rerun the Module 2 per-pathway propagation
     tests on the extended strata and report q-values alongside Module
     2's.
  4. High-coverage subset re-test — restrict to ``n_peptides >= 3`` and
     ``missing_fraction <= 0.2``; rerun the per-pathway propagation
     tests.
  5. NCC/SPAK phosphosite observability check — confirm the
     pre-specified regulatory sites are NOT in the low-observability
     tail of the phosphoproteome.

Usage::

    python scripts/v11/06_observability_audit.py [--n-null 10000]
"""
```

## `src/data/build_id_map.py`

### Module docstring, original lines 2-24

```python
"""
Build Ensembl → Symbol ID Mapping
==================================

Queries the Ensembl REST API to map mouse Ensembl gene IDs to MGI symbols.
Generates cached TSVs for use by Phase 7 biological grounding.

Outputs:
  id_map.tsv          — universe genes only (safe for enrichment denominator)
  id_map_extras.tsv   — extra pathway symbols (annotation/reporting only)
  id_map.json         — combined dict for quick lookups

Uses two strategies:
  1) POST /lookup/id (batch, fast) — maps Ensembl → display_name
  2) GET /xrefs/symbol/mus_musculus/<symbol> (extra pathway symbols)
     followed by POST /lookup/id to get canonical display_name

Usage:
    python scripts/build_id_map.py \\
        --genes data/processed/networks/run_xxx/phase2_genes.txt \\
        --outdir data/processed/resources \\
        --extra_symbols Wnk1,Slc12a3,Kcnj10
"""
```

## `src/enrichment/__init__.py`

### Module docstring, original lines 1-8

```python
"""
RRRM-2 Enrichment Module

Phase 7: Biological Grounding

This module contains:
- biological_grounding: Gene set enrichment + pathway analysis
"""
```

## `src/markers/discover_dct.py`

### Module docstring, original lines 2-27

```python
"""
Phase 1.5: Dataset-Derived DCT Marker Discovery
=================================================

Identifies genes whose expression tracks the deconvolved DCT fraction,
after controlling for other nephron segments and experimental design.

KEY DESIGN DECISION:
  Uses VST-normalised expression (Y), NOT Rtech.
  Rtech has CLR segment proportions already regressed out (DESeq2.R line 224),
  so β_DCT ≈ 0 on Rtech. We need the pre-segment-residual expression.

Model (per gene g, sample i):
  Y_ig = α + β_DCT · CLR(DCT)_i + γ' · CLR(other)_i + θ_cell(i) + ε_ig

  where cell(i) = Age × Arm × EnvGroup (16 one-hot dummies, 15 df)

β_DCT > 0, q(BH) < 0.05, bootstrap_freq ≥ 0.70 → DCT marker gene

Usage:
    python scripts/discover_dct_markers.py \\
        --vst data/processed/vst_normalized/GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv \\
        --meta data/processed/phase1_residuals/meta_phase1.tsv.gz \\
        --clr data/processed/deconvolution/music_segment_direct_proportions_CLR.csv \\
        --outdir data/processed/dct_markers
"""
```

### Comment block, original lines 403-405

```python
    # If no single 'DCT' column, synthesize from DCT subtypes.
    # The hybrid reference splits DCT into DCT1, DCT2, CNT — we need a
    # combined DCT signal for marker discovery.
```

### Comment block, original lines 413-416

```python
            # CLR values are log-ratio transformed, so we need to:
            # 1) back-transform to proportions, 2) sum, 3) re-CLR
            # But since CLR = log(p/geom_mean), summing CLR is wrong.
            # Instead, exponentiate, sum, then log.
```

### Comment block, original lines 446-448

```python
    # ── 5. Filter genes (CPM ≥ 1 in ≥ 20% of samples) ──────────────────
    # VST values are already normalised, but we should filter low-expression
    # Use a variance floor: keep genes with var > 0
```

### Comment block, original lines 492-501

```python
    # ── Anti-confounding via DELTA approach ─────────────────────────────
    # Problem: CLR(DCT), CLR(TAL_LOH), CLR(CD) are highly correlated because
    # the two-stage Scale & Combine multiplies all distal subtypes by the same
    # distal_total. Full residualization destroys >99% of DCT variance.
    #
    # Fix: Instead of residualizing DCT against TAL/CD (which kills variance),
    # use a per-gene delta filter: gene passes if its correlation with DCT
    # exceeds its max correlation with confounders. This catches genes that
    # track DCT specifically rather than just "distal signal in general".
    
```

### Comment block, original lines 532-534

```python
        # Anti-confounded = gene tracks DCT MORE than it tracks any confounder
        # Use signed correlation: gene must correlate positively with DCT
        # AND that positive correlation must exceed the max confounder correlation
```

### Comment block, original lines 539-542

```python
        # Also store a "DCT_resid" correlation for reference (using partial correlation)
        # Partial r(Y, DCT | confounders) via formula:
        # r_partial = (r_xy - r_xz * r_yz) / sqrt((1 - r_xz²)(1 - r_yz²))
        # For multiple confounders, use the strongest one
```

### Comment block, original lines 624-627

```python
    # ── 10. Final panel ─────────────────────────────────────────────────
    # Core criteria: β_DCT > 0, q < threshold, bootstrap stable
    # NEW: Sign consistency AND Anti-confounding
    
```

### Comment block, original lines 657-659

```python
        # f3. Anti-confounding delta — use the STRICTER fallback_delta, not diff_threshold
        #     (diff_threshold=0.0 is fine for OLS which already controls for confounders;
        #      the marginal path has NO confounder control except this δ gate)
```

## `src/markers/discover_markers.py`

### Module docstring, original lines 2-29

```python
"""
Phase 1.5b: Generalized Segment Marker Discovery
==================================================

Discovers marker genes for ALL nephron segments in the CLR file, using the
same regression-based approach as discover_dct.py but parameterized by
target segment.

For each segment S in the CLR file:
  1. Fit OLS: Y_ig = α + β_S · CLR(S)_i + γ' · CLR(other)_i + θ_cell(i) + ε_ig
  2. Select genes where β_S > 0, BH q < 0.05 (or fallback to marginal r)
  3. Bootstrap stability within experimental cells
  4. Anti-confounding: gene tracks S more than any other segment

Outputs per segment:
  <outdir>/<segment>_marker_scores.tsv   — full regression results
  <outdir>/<segment>_marker_panel.txt    — final marker gene list

Usage:
    python -m src.markers.discover_markers \\
        --vst data/processed/vst_normalized/GLDS-674_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv \\
        --meta data/processed/phase1_residuals/meta_phase1.tsv.gz \\
        --clr data/processed/deconvolution/latest/music_segment_direct_proportions_CLR.csv \\
        --outdir data/processed/segment_markers

    # Single segment:
    python -m src.markers.discover_markers --segments DCT ...
"""
```

## `src/multiomics/celltype_panels.py`

### Module docstring, original lines 1-18

```python
"""Cell-type marker-panel decomposition (memo Recommendation 4).

Bulk kidney RNA cannot, by itself, distinguish a transcriptional *program*
change inside DCT cells from a change in nephron-segment *composition* or from
dilution by infiltrating/interstitial RNA. This module supplies curated kidney
cell-type marker panels and the scoring/decision helpers used to bound that
ambiguity without new single-cell data, following memo Section 8.2-8.3.

Key distinction the panels are designed to make:

* ``dct_identity`` -- DCT cell-type markers whose level tracks *how much DCT is
  present* (Pvalb, Trpm6, Calb1, ...).
* ``dct_transport`` -- the WNK-SPAK/OSR1-NCC functional transport program.

If ``dct_transport`` falls while ``dct_identity`` is preserved, the most likely
reading is transcriptional suppression of the transport program rather than
loss of DCT cells (memo: "strong support for DCT program interpretation").
"""
```

### Comment block, original lines 26-28

```python
# Curated mouse kidney marker panels. Sources: Chen et al. 2021 (GSE150338),
# Ransick et al. 2019 (kidney cell atlas), Park et al. 2018; standard segment
# markers. Symbols are mgi-style to match decoupler/id_map symbol space.
```

## `src/multiomics/osd462_anchor.py`

### Module docstring, original lines 1-28

```python
"""OSD-462 / RR-10 multi-omics anchor.

Dataset-agnostic estimation and inference helpers used by the
``scripts/osd462/*`` layer scripts to test whether the RRRM-2 RNA-level
matrix-high / DCT-low remodeling signal recurs at the protein and
phosphoprotein level in the independent OSD-462 spaceflight kidney cohort.

Design notes
------------
* Cross-study, cross-strain, cross-modality: no sample-level pooling.  The
  unit of comparison is the *direction of the flight effect* (Space Flight
  minus Ground Control), exactly the contrast-vector logic already used for
  OSD-513 in ``src/networks/cross_osdr_projection.py``.
* TMT 2-plex batch structure ("Samp1-5", "Samp6-10").  Plex is handled as a
  batch factor everywhere: flight effects are estimated *within* each plex and
  then averaged, which removes any per-plex (and, after channel centering,
  per-channel) loading constant.
* Primary inference is an abundance/peptide-matched random-gene-set null,
  which is robust to the generically weak RNA<->protein correlation and to
  gene-set size.

The TMT workbook layout (both proteomics and phosphoproteomics) is::

    row 1 : group banners ("... scaled ..." marks the scaled blocks)
    row 2 : per-column sample labels  (BL-01, FL-03, GC-05, ...)
    row 3 : machine column headers     (Samp1-5~rq_129n_sn scaled, ...)
    row 4+: one protein / phosphosite per row
"""
```

### Comment block, original lines 39-42

```python
# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 49-52

```python
# ─────────────────────────────────────────────────────────────────────────────
# TMT workbook parsing
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 201-204

```python
# ─────────────────────────────────────────────────────────────────────────────
# Flight-effect estimation (within-plex FL - GC, plex-averaged)
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 415-418

```python
# ─────────────────────────────────────────────────────────────────────────────
# Gene-set statistics
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 446-449

```python
# ─────────────────────────────────────────────────────────────────────────────
# Abundance / peptide-matched random-gene-set null
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 548-551

```python
# ─────────────────────────────────────────────────────────────────────────────
# Pathway-vector cosine recurrence (Layer 4)
# ─────────────────────────────────────────────────────────────────────────────
```

## `src/multiomics/phenotype_anchor.py`

### Module docstring, original lines 1-17

```python
"""Phenotype-anchoring layer: animal-matched RNA-state vs NCC-activity.

OSD-462 / RR-10 ran kidney RNA-seq and phosphoproteomics on the *same animals*
(channel ``FL-01`` <-> RNA sample ``FLT_F1`` etc.). This module builds a
per-animal NCC/SPAK regulatory-phosphorylation activity score and a per-animal
DCT/NCC-low RNA state score, matches them by animal, and compares them at three
levels of stringency:

1. group level (Space Flight vs Ground Control) -- the primary, robust claim;
2. all-sample correlation -- descriptive, inflated by the group means;
3. condition-adjusted correlation -- the real test of an animal-level link,
   on within-condition-centered residuals; underpowered at n~10/condition.

Controls: a non-regulatory NCC phosphosite score (should not track the RNA
state) and an RNA score with the Slc12a3 transcript removed (the link must not
be one gene covarying with its own phospho).
"""
```

## `src/multiomics/regulator_activity.py`

### Module docstring, original lines 1-17

```python
"""Regulator-activity prioritization (v10 trimmed plan).

Two activity-inference layers added on top of the verified v9/OSD-462 analysis:

* Layer A -- kinase-substrate enrichment (KSEA) on OSD-462 phosphoproteomics.
  Implemented as the classic Casado et al. (2013) z-score statistic, fully
  self-contained (no external database required for the curated core panel;
  a PhosphoSitePlus-format table can be supplied for the full panel).
* Layer B -- transcription-factor / pathway activity inference on RNA, using
  ``decoupler`` (ULM) with PROGENy and DoRothEA/CollecTRI priors. The decoupler
  resources are fetched from omnipath/zenodo and therefore require network
  access at run time; the computation itself is offline.

Both layers are *prioritization*, not causal mechanism discovery. The wording in
all outputs is deliberately constrained: "candidate upstream organizer",
"activity anchor", "negative boundary" -- never "mechanism" or "causal".
"""
```

### Comment block, original lines 32-35

```python
# ---------------------------------------------------------------------------
# Layer A -- KSEA (self-contained, no external dependency)
# ---------------------------------------------------------------------------
```

### Comment block, original lines 152-155

```python
# ---------------------------------------------------------------------------
# Layer B -- TF / pathway activity inference (decoupler ULM)
# ---------------------------------------------------------------------------
```

### Comment block, original lines 219-222

```python
# ---------------------------------------------------------------------------
# Integration -- evidence grading
# ---------------------------------------------------------------------------
```

## `src/networks/__init__.py`

### Module docstring, original lines 1-13

```python
"""
RRRM-2 Networks Module

Phase 2: Shared Topology Construction + LIONESS + Edge Regression
Phase 3: node2vec Embeddings + Procrustes Alignment

This module contains:
- shared_topology: Cell-standardized shared skeleton construction
- lioness: LIONESS sample-specific edge weights
- edge_regression: Edge-wise regression with full factorial design
- embeddings: PecanPy node2vec embeddings
- procrustes: Procrustes alignment + rewiring metrics
"""
```

## `src/networks/contrast_vectors.py`

### Module docstring, original lines 1-21

```python
"""Contrast-vector decomposition core (Cross-OSDR Network-Contrast Framework).

Implements the geometry of §1 of agents_instruction.md:
    A^a_GC  = N(Old, a, GC)  - N(Young, a, GC)
    A^a_FLT = N(Old, a, FLT) - N(Young, a, FLT)
    beta_a  = (A_FLT . A_GC) / (A_GC . A_GC)              # projection coefficient
    cos_a   = (A_FLT . A_GC) / (||A_FLT|| ||A_GC||)
    R_a     = A_FLT - beta_a * A_GC                       # redirected component
    rho_a   = ||R_a|| / ||A_FLT||                         # redirectedness

Optional precision (variance) weighting per Guardrail C (§2.3):
    beta_w = (A_FLT^T W A_GC) / (A_GC^T W A_GC)
    cos_w  = (A_FLT^T W A_GC) / sqrt(A_FLT^T W A_FLT * A_GC^T W A_GC)

Bootstrap (§2.4) and permutation (§4.4) helpers operate on
stratified row-indices over the sample table so callers control how subjects
are nested in cells.

This module is dataset-agnostic. Adapters that build N(condition) from
expression/network artifacts live in bootstrap_decomposition.py.
"""
```

### Comment block, original lines 33-36

```python
# ─────────────────────────────────────────────────────────────────────────────
# Core geometry
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 151-154

```python
# ─────────────────────────────────────────────────────────────────────────────
# Interpretation categories (§4.4)
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 234-237

```python
# ─────────────────────────────────────────────────────────────────────────────
# Stratified bootstrap and label permutation
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 276-279

```python
# ─────────────────────────────────────────────────────────────────────────────
# Vector builders + bootstrap orchestration
# ─────────────────────────────────────────────────────────────────────────────
```

## `src/networks/cross_osdr_projection.py`

### Module docstring, original lines 1-7

```python
"""Cross-OSDR contrast-vector alignment utilities.

This module implements the contrast-level pooling rule from
agents_instruction.md: external studies are reduced to within-study flight
vectors first, then compared by cosine/projection statistics. It contains no
raw cross-study expression pooling.
"""
```

## `src/networks/edge_regression.py`

### Comment block, original lines 122-124

```python
    # Covariates (numeric vs categorical). Residualized expression is already
    # adjusted for batch/deconvolution/SVs; adding those terms again would
    # double-adjust the network weights.
```

### Comment block, original lines 219-222

```python
    # Define contrasts
    
    # Helper to build coefficient names matching interaction(Age, Arm, EnvGroup) pattern
    # R's interaction() produces "Level1.Level2.Level3", then make.names() sanitizes
```

## `src/networks/external_aging_axis.py`

### Module docstring, original lines 1-7

```python
"""Guardrail E: project RRRM-2 flight vectors onto an external aging axis.

The implementation is deliberately feature-agnostic: callers provide a flight
effect vector and an independently estimated aging-axis vector indexed by the
same feature IDs. The CLI provides the gene-level RRRM-2 path used by the
contrast-vector orchestrator.
"""
```

## `src/networks/lar_reversal.py`

### Module docstring, original lines 1-13

```python
"""LAR reversal and mechanism-switch analysis utilities.

This module tests three competing interpretations of the RRRM-2 LAR arm:

* attenuation: LAR flight is near zero relative to ISS-T flight,
* reversal: LAR flight moves opposite to the ISS-T flight direction,
* mechanism switch: LAR does not simply negate ISS-T, but loads onto other
  biological axes such as clock/DCT, S1P, immune, or preservation-stress context.

The functions here are deliberately sample-level and contrast-level.  Gene and
network outputs are used as exploratory annotation, not confirmatory topology
claims.
"""
```

### Comment block, original lines 29-31

```python
# Independent state-space coordinates used for cosine/projection geometry.
# matrix_minus_dct is a scalar readout, not an independent coordinate, because
# it is exactly matrix_component - dct_transport_component.
```

## `src/networks/mechanism_axis.py`

### Module docstring, original lines 1-7

```python
"""Mechanism-axis scoring and exploratory network-neighborhood prioritization.

This module keeps the statistical pieces behind
``scripts/run_mechanism_axis_prioritization.py`` small enough to test.  The
analysis is intentionally pathway/module first; gene rankings produced here are
secondary prioritizations, not confirmatory gene-level rewiring calls.
"""
```

## `src/networks/stability_test.py`

### Module docstring, original lines 1-12

```python
"""Guardrail A — bootstrap stability of the control aging vector A^a_GC.

Per agents_instruction.md §2.1, the within-RRRM-2 projection layer is only
allowed when the GC aging direction is itself stable. This module quantifies
that stability and emits the pass/fail decision artifacts that downstream
phases gate on.

Key API:
    estimate_agc_stability(...)        -> StabilityReport
    apply_stability_gate(...)          -> StabilityDecision
    write_stability_artifacts(...)     -> persists report + decision JSON
"""
```

## `src/networks/tubulointerstitial_state.py`

### Module docstring, original lines 1-7

```python
"""Tubulointerstitial state-space analysis utilities.

This layer translates mechanism scores into a biological state-space:
matrix remodeling on the x-axis and native DCT/NCC-WNK transport on the
y-axis.  The DCT component is never sign-flipped; the scalar remodeling
summary is defined separately as matrix_minus_dct.
"""
```

## `src/networks/wgcna_analysis.R`

### Comment block, original lines 2-21

```r
# src/networks/wgcna_analysis.R
# ─────────────────────────────────────────────────────────────────────────
# WGCNA Module Discovery, Module-Trait Association, and Module Preservation
#
# Replaces the LIONESS/node2vec/Procrustes stack with:
#   Analysis A: Module discovery on FLT+GC samples (n≈40)
#   Analysis B: Module eigengene ~ Flight * Age * Arm association
#   Analysis C: Module preservation (GC ref → FLT test)
#   Analysis D: Pathway enrichment per module
#
# Usage:
#   Rscript src/networks/wgcna_analysis.R \
#     --rtech data/results/run_XYZ/phase1_residuals/Rtech.tsv.gz \
#     --meta  data/results/run_XYZ/phase1_residuals/meta_phase1.tsv.gz \
#     --id_map data/processed/resources/id_map.tsv \
#     --gene_sets config/gene_sets.yaml \
#     --outdir data/results/run_XYZ/wgcna \
#     --max_genes 5000 --n_pres_perms 200
# ─────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 147-150

```r
# ═══════════════════════════════════════════════════════════════════════════
# Analysis A: Soft Threshold Selection + Module Discovery
# ═══════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 262-265

```r
# ═══════════════════════════════════════════════════════════════════════════
# Analysis B: Module-Trait Association
# ═══════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 337-340

```r
# ═══════════════════════════════════════════════════════════════════════════
# Analysis C: Module Preservation (GC ref → FLT test)
# ═══════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 357-360

```r
# NOTE: Using cor (Pearson) for preservation instead of bicor.
# bicor + corOptions triggers a known WGCNA bug in modulePreservation
# where maxPOutliers is passed as a positional arg causing dimension errors.
# Modules are still built with bicor; only preservation testing uses Pearson.
```

### Comment block, original lines 430-433

```r
# ═══════════════════════════════════════════════════════════════════════════
# Analysis D: Module Pathway Enrichment
# ═══════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 518-521

```r
# ═══════════════════════════════════════════════════════════════════════════
# Save summary metadata
# ═══════════════════════════════════════════════════════════════════════════
```

## `src/networks/wgcna_followup.py`

### Comment block, original lines 117-119

```python
# ══════════════════════════════════════════════════════════════════════════
# 1. SIMPLE-EFFECT CONTRASTS
# ══════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 136-142

```python
# For ME ~ Flight * AgeOld * ArmLAR:
# Flight effect at (AgeOld=a, ArmLAR=r) = 
#   β_Flight + a*β_Flight:AgeOld + r*β_Flight:ArmLAR + a*r*β_Flight:AgeOld:ArmLAR
# Contrast vector c = [0, 1, 0, 0, a, r, 0, a*r] (matching intercept, Flight, AgeOld, ArmLAR, F:A, F:R, A:R, F:A:R)
```

### Comment block, original lines 198-200

```python
# ══════════════════════════════════════════════════════════════════════════
# 2. GREY60 HUB GENES + kME
# ══════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 237-239

```python
# ══════════════════════════════════════════════════════════════════════════
# 3. EIGENGENE PLOTS
# ══════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 308-310

```python
# ══════════════════════════════════════════════════════════════════════════
# 4. BROADER ENRICHMENT (GO/KEGG via gseapy)
# ══════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 365-367

```python
# ══════════════════════════════════════════════════════════════════════════
# 6. FINAL INTEGRATED TABLE
# ══════════════════════════════════════════════════════════════════════════
```

## `src/preprocessing/__init__.py`

### Module docstring, original lines 1-13

```python
"""
RRRM-2 Preprocessing Module

Phase 0: Pre-processing, Deconvolution, QC
Phase 1: Global Residualization (VST + SVA)

This module contains:
- data_alignment: Metadata-to-counts alignment
- export_counts: Raw count export utilities

Note: R scripts (residualization.R, deconvolution.R) are called via Rscript
and are stored alongside Python modules for organizational purposes.
"""
```

## `src/preprocessing/common_gene_universe.py`

### Module docstring, original lines 1-5

```python
"""Common gene-universe utilities for the Cross-OSDR framework.

The implementation only intersects already processed within-study matrices.
It does not pool raw expression across studies.
"""
```

## `src/preprocessing/data_alignment.py`

### Module docstring, original lines 2-15

```python
"""
Metadata Alignment Script for GLDS-674 / OSD-771
Aligns metadata to count matrices (STAR, RSEM, VST) and creates aligned datasets.

This script:
1. Reads count files (STAR unnormalized, RSEM unnormalized, VST normalized)
2. Reads ISA metadata
3. Aligns metadata to samples across all count matrices
4. Optionally parses design factors from sample names
5. Saves aligned outputs to data/processed/aligned_outputs/

Usage:
    python scripts/align_metadata_to_counts.py
"""
```

## `src/preprocessing/deconvolution.R`

### Comment block, original lines 136-138

```r
# 2) Load single-cell reference from H5AD
# 2) Load single-cell reference from Raw MTX
# Use hybrid reference (TMS + Chen) for improved DCT subtype resolution
```

### Comment block, original lines 158-163

```r
# Handle potential transpose issues with mtx (Scanpy writes shape n_obs x n_vars, readMM reads as such)
# But SCE expects genes x cells (n_vars x n_obs). Scanpy usually writes cells x genes?
# Creating scipy.io.mmwrite(..., raw_adata.X) usually writes standard MM.
# If raw_adata.X is (cells, genes), then readMM returns (rows=cells, cols=genes).
# SCE needs (rows=genes, cols=cells). So we likely need to Transpose.
# Let's check dimensions carefully.
```

### Comment block, original lines 309-314

```r
# ── DCT Sub-typing ──
# With hybrid reference (Chen + TMS), DCT subtypes are pre-annotated:
#   DCT1 (Slc12a3+/Pvalb+, thiazide-sensitive NCC)
#   DCT2 (Slc12a3+/Pvalb−, calcium transport)
#   CNT  (connecting tubule)
# No sub-clustering needed when using hybrid reference.
```

### Comment block, original lines 671-673

```r
# CRITICAL: donor_id diagnostics
# MuSiC uses cross-donor variance for cell-size (theta) estimation.
# If donor_id is NA, constant, or unique per cell, resolution fails.
```

### Comment block, original lines 717-720

```r
# Option to force single donor for testing
# Uncomment the line below to bypass donor issues:
# colData(sce2)$donor_id <- "D1"
```

### Comment block, original lines 733-736

```r
# In-silico pseudo-bulk test for MuSiC
# Requires: sce2 with counts, colData(sce2)$segment_use and $clusterType and $donor_id
# Also requires: clusters.type and group.markers as in your script
```

### Comment block, original lines 739-741

```r
# ---------
# Helpers
# ---------
```

### Comment block, original lines 759-761

```r
# ---------
# IMPROVED: Pseudobulk with exact proportions and all fixes
# ---------
```

### Comment block, original lines 832-835

```r
    # FIX 3: Independent DCT variation via two-level design
    # Sample tubule proportions with independent ranges

    # Step 1: Sample coarse group totals
```

### Comment block, original lines 980-982

```r
# ---------
# Build pseudo-bulk mixtures
# ---------
```

### Comment block, original lines 1015-1017

```r
  # Two ground-truth matrices:
  # P_cell: fraction of sampled cells from each segment
  # P_rna:  fraction of total UMIs contributed by each segment
```

### Comment block, original lines 1083-1085

```r
# ---------
# Experiment B: Per-donor pseudo-bulk (sample from single donor per mixture)
# ---------
```

### Comment block, original lines 1198-1200

```r
# ---------
# Run pseudo-bulk test (with normalize_flag parameter)
# ---------
```

### Comment block, original lines 1417-1419

```r
# align to pseudo-bulk sample names
# MuSiC returns rows = bulk samples; our pb_mat cols are samples
# So rownames(P_hat) should match colnames(pb_mat)
```

### Comment block, original lines 1468-1476

```r
# ══════════════════════════════════════════════════════════════════════════════
# TWO-STAGE DECONVOLUTION
# Stage 1: TMS-only MuSiC for coarse types (PT, Podo, Endo, Mes, Immune, Fibro, CD)
# Stage 2: Chen-only MuSiC for distal subtypes (TAL_LOH, DCT1, DCT2, CNT)
# Then scale Stage 2 relative proportions by distal_total from Stage 1.
# This avoids the block-diagonal donor problem where cross-donor variance
# between TMS segments and Chen segments cannot be modeled.
# ══════════════════════════════════════════════════════════════════════════════
```

### Comment block, original lines 1509-1515

```r
# MERGE DCT1 + DCT2 → DCT for Stage 2 deconvolution.
# Rationale: DCT1 and DCT2 share ~60% of canonical markers (Slc12a3, Calb1,
# Trpm6, Slc8a1, Trpv5, Fxyd2). The only DCT1-exclusive marker is Pvalb,
# which is insufficient for MuSiC to resolve a 4-type NNLS reliably.
# Stage 2 output: DCT1≈0, DCT2≈0.48 — a collapse artifact.
# Merging gives MuSiC 3 well-separated types (TAL_LOH, DCT, CNT) instead of 4
# near-collinear types.
```

### Comment block, original lines 1560-1568

```r
# PROBLEM: Stage 1 (TMS-only) has no distal tubule types (TAL, DCT1/2, CNT),
# so MuSiC absorbs all distal signal into PT and CD. The residual
# (1 - sum(Stage1)) is ~0 for most samples, zeroing out all distal types.
#
# FIX: Use a data-informed floor for distal_total. In mouse kidney, the distal
# nephron (TAL + DCT1 + DCT2 + CNT) represents ~20-40% of tubular epithelium.
# We set a floor so that Stage 2 relative proportions are always scaled by a
# meaningful amount, then re-normalize the final proportions.
```

### Comment block, original lines 1572-1575

```r
# Data-informed distal floor: use the median of Stage 2's own "confidence"
# Stage 2 gives relative proportions. If they're spread across 4 types,
# the reference thinks the bulk genuinely has distal signal. We combine this
# with a minimum floor of 10% (known kidney anatomy).
```

### Comment block, original lines 1696-1702

```r
# FIX APPLIED 2026-01-11: CPM Normalization for Validation
# Problem: prop_seg is library-size independent (fractions sum to 1), but raw
#          bulk counts scale with sequencing depth. This caused spurious low/negative
#          correlations in validation (e.g., DCT showed ρ=0.03 with raw, ρ=0.33 with CPM).
# Solution: Normalize bulk to CPM before computing marker scores.

# CPM normalize bulk matrix for validation
```

### Comment block, original lines 1797-1799

```r
# DCT-specific marker check (DCT1/DCT2/CNT/CD subtypes)
# Expanded canonical marker lists (15–25 per segment) to ensure ≥10 survive
# symbol→Ensembl mapping, making rho estimates reliable.
```

### Comment block, original lines 1912-1915

```r
# 2. Compute immune marker scores from bulk
# Expanded immune markers: pan-immune (Ptprc/Lyz2/Tyrobp/Ctss/Fcer1g),
# macrophage (Cd68/Adgre1/Csf1r), T-cell (Cd3d/Cd3e), B-cell (Ms4a1),
# monocyte (Lst1/Itgam), MHC-II (H2-Aa)
```

### Comment block, original lines 2006-2009

```r
# 6) Use segment-direct proportions (already at segment level)
# Since we ran MuSiC on segments directly, prop_seg is already at segment granularity
# No need to aggregate cell types - the output is already segment-level
```

## `src/preprocessing/deconvolution_sanity.py`

### Comment block, original lines 88-90

```python
    # Align samples
    # VST columns match CLR index?
    # VST cols = sample names. CLR index = sample names.
```

## `src/preprocessing/deconvolution_sensitivity.R`

### Comment block, original lines 3-10

```r
# MuSiC reference-atlas sensitivity analysis.
#
# This script compares primary MuSiC cell-type proportions against an alternate
# deconvolution/reference-free estimate (for example bMIND, CIBERSORTx-compatible
# output, or a TMS pseudo-bulk holdout). It also optionally compares key
# downstream Phase 3 top-decile rewiring and Phase 7 pathway summaries produced
# after rerunning residualization with alternate proportions.
```

## `src/preprocessing/export_phase1.R`

### Comment block, original lines 1-4

```r
# scripts/export_phase1_to_python.R
# Export Phase 1 residualized expression data to Python-friendly formats

# Parse --outdir argument
```

## `src/preprocessing/multi_study_harmonization.R`

### Comment block, original lines 2-30

```r
# Multi-study residualization + harmonization (Phase 1 of the Cross-OSDR framework).
#
# Per agents_instruction.md §5 Phase 1, this script:
#   1. Residualizes each MV study (OSD-771, OSD-513, OSD-253) INDEPENDENTLY
#      against estimated kidney cell composition + per-study technical covariates.
#      Raw expression is NEVER pooled across studies (§3 absolute pooling rule).
#   2. Restricts each study to the common protein-coding gene universe and
#      writes per-study residual matrices for downstream contrast-vector code.
#   3. Emits a per-study manifest that records sample counts, covariates, and
#      the formula used.
#
# Usage:
#   Rscript src/preprocessing/multi_study_harmonization.R \
#       --counts=data/processed/<study>/raw_counts.tsv.gz \
#       --meta=data/processed/<study>/metadata.tsv \
#       --clr=data/processed/<study>/music_composition_CLR.csv \
#       --gene_universe=data/processed/multi_study_gene_universe.tsv \
#       --study=OSD-513 \
#       --outdir=data/processed/OSD-513
#
# Downstream consumers:
#   - scripts/run_contrast_vector_framework.py  (Phase 3 stability gate + Phase 4 decomposition)
#   - src/networks/cross_osdr_projection.py     (Phase 6 alignment)
#
# This file is intentionally lightweight; the heavy preprocessing per study is
# already done elsewhere (e.g. nf-core/rnaseq, DESeq2 VST). The job here is to
# (a) regress out per-study cell composition + technical covariates and
# (b) restrict to the common gene universe.
```

## `src/run_all_phases.py`

### Module docstring, original lines 2-16

```python
"""
RRRM-2 Full Pipeline Runner
===========================
Runs all phases of the RRRM-2 network rewiring analysis pipeline.

Usage:
    python scripts/run_all_phases.py                    # Run all phases
    python scripts/run_all_phases.py --phases 2 3       # Run specific phases
    python scripts/run_all_phases.py --skip-r           # Skip R-dependent steps
    python scripts/run_all_phases.py --dry-run          # Show commands without running

Requires:
    - Python environment with dependencies from requirements.txt
    - R with DESeq2, limma, sva, MuSiC (for Phase 0-1 and edge regression)
"""
```

### Comment block, original lines 1334-1337

```python
    # ── WGCNA routing ─────────────────────────────────────────────────
    # When --network-method wgcna, replace LIONESS phases (2,3,5,6,7)
    # with a single WGCNA phase.  Phases 0, 1, 1.5, 8, 8.5, 10, 9 are
    # unchanged (preprocessing, validation, diagnostics, figures).
```

## `src/statistics/__init__.py`

### Module docstring, original lines 1-14

```python
"""
RRRM-2 Statistics Module

Phase 5: Rewiring Metrics + Silent Shifters
Phase 6: Uncertainty/Null Models (Permutation + Bootstrap)

This module contains:
- silent_shifters: Silent shifter identification (high rewiring, low DE)
- interaction_metrics: Interaction persistence analysis
- permutation_bootstrap: Permutation/bootstrap uncertainty estimation
- full_regression: Full edge regression (all 80 samples)

Note: R script (differential_expression.R) is called via Rscript.
"""
```

## `src/statistics/bootstrap_decomposition.py`

### Module docstring, original lines 1-6

```python
"""Bootstrap orchestration helpers for contrast-vector decomposition.

This module is the statistics-facing wrapper around
``src.networks.contrast_vectors``. It keeps file/summary conventions in one
place so Phase 3-7 runners emit consistent artifacts.
"""
```

## `src/statistics/differential_expression.R`

### Comment block, original lines 1-12

```r
# scripts/run_gene_level_DE.R
#
# Generate gene-level differential expression tables for FLT vs GC contrasts
# within each Age×Arm stratum. These outputs can be used with
# phase5_build_silent_shifters_strict.py --gene_de
#
# Outputs: data/processed/gene_level_DE/
#   - ISS_T_YNG_FLT_vs_GC_gene_DE.tsv
#   - ISS_T_OLD_FLT_vs_GC_gene_DE.tsv
#   - LAR_YNG_FLT_vs_GC_gene_DE.tsv
#   - LAR_OLD_FLT_vs_GC_gene_DE.tsv
```

## `src/statistics/direct_coexpression_test.py`

### Module docstring, original lines 2-32

```python
"""
Direct Differential Co-expression Diagnostic (Phase 10)
========================================================

Permutation-based aggregate diagnostic to determine whether LIONESS-derived
rewiring rankings reflect group-level covariance changes or are primarily
artifacts of the sample-specific network pipeline.

For each Age×Arm stratum:
  1. Compute per-edge Δz = atanh(r_FLT) − atanh(r_GC) on the shared skeleton.
     Individual edge p-values are NOT trusted (n=5 gives SE≈1).
  2. Aggregate to gene-level: S_g = mean(|Δz_e|) over incident edges.
  3. Build a label-permutation null (K shuffles of FLT/GC within stratum)
     recomputing the entire gene score each time.
  4. Compare direct-correlation gene ranking to LIONESS/node2vec rewiring
     ranking (Spearman ρ), with a permutation null for the ρ itself.
  5. Run pathway enrichment on the direct ranking and compare with
     LIONESS-derived pathway results.

The skeleton was built using all samples (including FLT/GC), so this is
labeled as a robustness diagnostic, not primary inference.

Usage:
    python -m src.statistics.direct_coexpression_test \\
        --rtech  data/results/run_XYZ/phase1_residuals/Rtech.tsv.gz \\
        --meta   data/results/run_XYZ/phase1_residuals/meta_phase1.tsv.gz \\
        --phase2_dir data/results/run_XYZ/networks \\
        --rewiring_dir data/results/run_XYZ/phase3_rewiring \\
        --outdir data/results/run_XYZ/phase10_direct_coexpr \\
        --n_perms 1000
"""
```

## `src/v11/aldosterone_axis.py`

### Module docstring, original lines 2-39

```python
"""Reach F -- explicit aldosterone / mineralocorticoid-axis test.

Motivation (manuscript glue). NCC is aldosterone-regulated and spaceflight's
cephalad fluid shift perturbs the renin-angiotensin-aldosterone system, yet the
manuscript never tests the aldosterone axis -- even though the regulator run
flagged a steroid-receptor ("Androgen") PROGENy pathway as up (MR and androgen
receptors share hormone-response elements). This module scores a curated
mineralocorticoid panel across the same five cohorts used by the recurrence
analysis, cross-checks it against the already-computed PROGENy activities, and
asks whether the aldosterone-axis direction is *coherent* with the predicted
distal-nephron (NCC/WNK-SPAK) suppression -- reported side by side with the
transport-anchor score.

Design notes.
  * The MR panel (``data/external/gene_sets/mineralocorticoid_axis_core.tsv``)
    is directional: ``up`` genes are aldosterone-induced (Sgk1, Tsc22d3, Per1,
    Klf15, Fxyd2, Scnn1a/b/g, Slc12a3); ``context`` genes (Hsd11b2, Nr3c2,
    Wnk1/4, Stk39, Oxsr1) gate or organize the axis and are reported but
    excluded from the directional score.
  * Scoring is decoupleR-free on purpose: the panel is small and curated, the
    sandbox has no decoupler, and a direction-corrected mean with a gene-label
    permutation null is fully reproducible. Cross-cohort pooling reuses the
    project's ``signed_stouffer_z`` and ``recurrence_class``.
  * Sign convention: a *positive* axis effect means the aldosterone-responsive
    program moves up; the manuscript's distal-nephron suppression prediction is
    a *negative* axis effect that tracks NCC/WNK-SPAK transport suppression.

Discipline (failure mode is informative). Bulk RNA cannot separate systemic
hormone exposure from intrinsic DCT response; a flat or incoherent result
argues against an endocrine-driven story and narrows the mechanism. Nothing
here is a mechanistic claim.

Inputs (on disk): the five cohort gene-stat tables under the regulator run's
``rna_effects/`` (gene, stat), the curated MR panel, and the existing
``rna_progeny_pathway_activity.tsv``.
Output: ``regulator_activity/aldosterone_axis_summary.tsv`` (+ verdict JSON).
Provenance key for the headline index: ``aldo_axis_meta_effect``.
"""
```

### Comment block, original lines 65-67

```python
# --------------------------------------------------------------------------- #
# Pure scoring (unit-tested; no I/O).                                          #
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 170-172

```python
# --------------------------------------------------------------------------- #
# Orchestration (I/O).                                                         #
# --------------------------------------------------------------------------- #
```

## `src/v11/channel_centering_sensitivity.py`

### Module docstring, original lines 2-14

```python
"""TMT channel-centering QC for the v11 phosphosite enrichment analysis.

This module answers two audit questions around the OSD-462 TMT phosphosite
layer:

1. Does the raw phosphosite channel pattern show lower flight-channel medians?
2. Does the DCT-subtype-prior enrichment persist when phosphosite effects are
   recomputed without the within-plex channel-centering step?

The primary analysis uses channel-centered phosphosite effects written by
``scripts/osd462/02_phospho_axis.py``.  The uncentered recomputation here is a
sensitivity check, not a replacement estimator.
"""
```

## `src/v11/cmap_screen.py`

### Module docstring, original lines 2-13

```python
"""Reach E -- hypothesis-generating LINCS/CMap appendix screen.

This is not a treatment-discovery analysis.  It builds a conservative human
L1000-compatible query from the Repair B cross-cohort mouse RNA meta-analysis
and, when the required LINCS metadata are present, computes an approximate local
connectivity score against the Level-5 matrix.

Critical guardrail: local signature scores are not interpretable without
``sig_info`` because the signature id alone does not tell us the perturbagen,
cell line, dose, or time.  If ``sig_info`` is missing, this module writes only
``query_genes.tsv`` plus a verdict marking the local screen as blocked.
"""
```

## `src/v11/core_analysis.py`

### Module docstring, original lines 2-10

```python
"""Execute the core v11 DCT-subtype-prior phosphoproteome analyses.

This script intentionally keeps the external references as priors/scaffolds:
GSE228367 is a DCT subtype RNA reference, PXD001729 is a cultured DCT-lineage
phosphoproteomic reference, and OSD-462 remains the only spaceflight
phosphoproteomic anchor. Historical output directories still use H2/DCT1 and
H3/mediation names for continuity, but manuscript-facing language treats them
as DCT-subtype-prior enrichment and exploratory covariance decomposition.
"""
```

## `src/v11/dct_continuous_gradient.py`

### Module docstring, original lines 2-39

```python
"""Repair C -- DCT1<->DCT2 *continuous* contrast for flight-suppressed phosphosites.

Motivation (audit response). The v11 H2 enrichment test bins the DCT-subtype
prior into top/bottom decile/quartile flags and asks whether flight-suppressed
phosphosites are over-represented in the DCT1-high bin. A reviewer can fairly
object that (a) hard binning discards most of the continuous DCT1<->DCT2
coordinate and (b) a bin-level signal could be carried by a handful of
NCC/SPAK/WNK anchor genes. This module answers both by regressing phospho
response on the *continuous* ``dct1_enrichment_score`` coordinate at the
parent-gene level, with:

  * a direct linear contrast (OLS slope, covariate-adjusted),
  * a natural-cubic-spline term to detect non-monotone structure the linear
    contrast would miss,
  * Spearman rho as a rank-based nonparametric backup, and
  * a logistic gradient on the binary suppression indicator,

each re-run with the NCC/SPAK/WNK anchor genes excluded and on
single-position sites only. It reuses ``build_parent_gene_table`` from
``core_analysis`` so the gene collapse, covariates, and coordinate definition
are identical to the binned H2 analysis -- only the *test* changes.

Discipline. ``dct1_enrichment_score`` is a GSE228367 DCT-subtype RNA prior, not
spaceflight evidence; this is an exploratory subtype-prior gradient test, not a
mechanistic claim. Input is the per-run artifact
``dct_prior/osd462_phosphosite_dct1_prior.tsv`` written by
``core_analysis.build_dct_prior_mapping``.

Sign convention (predicted, supportive direction):
  * phospho-effect outcomes (mean / most-negative): more DCT1 => more
    *negative* (suppressed) effect => supportive slope/rho is NEGATIVE.
  * suppression-indicator outcome: more DCT1 => higher P(suppressed) =>
    supportive log-odds is POSITIVE.

Output: ``h2_enrichment/h2_dct_continuous_gradient.tsv`` (+ a verdict JSON).
Provenance keys emitted for the headline index: ``dct_gradient_slope``,
``dct_gradient_p``, ``dct_gradient_spearman``.
"""
```

### Comment block, original lines 70-72

```python
# --------------------------------------------------------------------------- #
# Pure estimators (unit-tested; no I/O).                                       #
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 235-237

```python
# --------------------------------------------------------------------------- #
# Orchestration (I/O).                                                         #
# --------------------------------------------------------------------------- #
```

## `src/v11/h2_composition_aware_phospho.py`

### Module docstring, original lines 2-13

```python
"""Composition-aware H2 robustness tests for OSD-462 phosphoproteomics.

This script tests whether the v11 DCT-subtype-prior phosphosite signal survives
a conservative adjustment ladder using per-animal TMT phosphosite intensities.

Claim discipline:
  * This is not DCT-specific phosphoproteomics.
  * Composition scores are bulk RNA marker-panel estimates and may sit on the
    biological path from flight to phosphosite dilution.
  * The output should be described as composition-aware robustness, not
    deconvolution.
"""
```

## `src/v11/h2_occupancy_normalized_phospho.py`

### Module docstring, original lines 2-6

```python
"""Parent-protein-normalized OSD-462 phosphosite robustness analysis.

The historical output directory is named ``h2_occupancy`` for continuity, but
the manuscript-facing term is parent-protein-normalized phosphosite effect.
"""
```

## `src/v11/human_concordance.py`

### Module docstring, original lines 2-22

```python
"""Reach D -- human urine/fluid-axis concordance.

This module deliberately does *not* validate the mouse kidney omics signature in
human kidney tissue.  The human evidence is urine / clinical chemistry from the
NASA Twins Study, so the analysis is an analyte-level sign concordance check
against the mouse-derived fluid, mineralocorticoid, and distal-nephron transport
predictions.

Inputs on disk:
  * Twins Table S8: urinary electrolyte / biochemistry timepoints.
  * Twins SM PDF Fig. S4A: AQP2 / AGT / RENR urine-protein figure evidence.
  * Optional OSD-656 files: Inspiration4 urine inflammation/NULISAseq recovery
    context, never the main Reach D source.

Outputs:
  * ``human_concordance/twins_axis_concordance.tsv``
  * ``human_concordance/twins_axis_level_concordance.tsv`` (independent-axis sign test)
  * ``human_concordance/twins_aqp2_figure_evidence.tsv``
  * ``human_concordance/human_concordance_verdict.json``
  * optional OSD-656 catalog/summary files if an OSD-656 directory is present.
"""
```

### Comment block, original lines 55-62

```python
# ── Externalized, version-controlled prediction / scoring spec ───────────────
# The analyte predictions, scoring flags, and the digitized figure directions
# live in config YAML (auditable, citable) instead of as Python literals here.
#   * config/human_concordance_prereg.yaml  -- table predictions + Fig. S4A rows
#   * config/human_urine_marker_panel.yaml  -- OSD-656 context categories
# Observed directions for the table analytes are still COMPUTED from Table S8 at
# runtime; only the figure-row observed directions are stored, because Fig. S4A
# has no numeric deposit and was digitized from the SM PDF (see provenance).
```

## `src/v11/kinome_atlas_ksea.py`

### Module docstring, original lines 1-32

```python
"""Repair A -- kinome-wide KSEA from the Johnson 2023 Ser/Thr atlas.

Replaces the three-substrate curated WNK--SPAK/OSR1 net (which can only ever
*confirm* the gene that was hand-curated into it) with a motif-driven assignment
over the *entire* OSD-462 phosphoproteome, then runs the existing Casado-style
KSEA (:func:`src.multiomics.regulator_activity.ksea`). The question becomes
unbiased: of all 303 Ser/Thr kinases scoreable from substrate motifs, which come
back inferred-*down* in flight -- and do SPAK (atlas label ``STLK3``) and OSR1
fall out near the top without having been planted there?

Pipeline
--------
1. ``load_atlas_pssms`` -- read the 303 position-normalised, log2-scaled PSSMs
   (sheet ``ser_thr_all_norm_scaled_matrice``) into a (kinase x position x
   residue) tensor of log2 weights.
2. ``score_sites`` -- for each phosphosite's 13-mer motif, sum the log2 weights
   over the atlas -5..+4 frame (the central S/T at index 6 is unscored) to get a
   per-kinase log-odds score.
3. ``percentile_by_kinase`` + ``build_kinome_net`` -- convert each kinase's
   scores to a within-cohort percentile and assign a site to a kinase when its
   percentile clears ``assign_percentile`` (Johnson's standard >=90th-percentile
   call). This substitutes the OSD-462 phosphoproteome itself for Johnson's
   Ochoa reference distribution, which is not redistributed on disk -- a
   documented, honest adaptation.
4. ``ksea`` -- the existing, independently tested KSEA statistic, unchanged.
5. ``over_representation`` -- a parent-gene-aware one-sided Fisher cross-check:
   are a kinase's predicted substrate *genes* enriched for suppressed sites?

Positive control: STLK3 (SPAK) and OSR1 must return negative KSEA z (inferred
activity down), consistent with NCC-activating-cluster suppression in OSD-462.
The wording stays prioritisation, not mechanism.
"""
```

### Comment block, original lines 52-55

```python
# --------------------------------------------------------------------------- #
# On-disk locations
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 75-78

```python
# --------------------------------------------------------------------------- #
# Atlas / motif geometry
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 84-86

```python
#: Atlas labels for the manuscript's positive-control axis. SPAK is recorded in
#: the Johnson atlas as ``STLK3``; OSR1 keeps its name; the WNK panel is present
#: as WNK1/WNK3/WNK4 (WNK2 is absent from the Ser/Thr atlas).
```

### Comment block, original lines 97-100

```python
# --------------------------------------------------------------------------- #
# (1) atlas PSSMs
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 187-190

```python
# --------------------------------------------------------------------------- #
# (2) motif scoring
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 259-262

```python
# --------------------------------------------------------------------------- #
# (3) kinome substrate net
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 331-334

```python
# --------------------------------------------------------------------------- #
# (4) parent-gene-aware over-representation cross-check
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 393-396

```python
# --------------------------------------------------------------------------- #
# Inputs: join phospho effects to atlas motifs
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 452-455

```python
# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #
```

## `src/v11/observability_audit.py`

### Module docstring, original lines 1-41

```python
"""Module 3 — proteome observability-bias audit.

Pre-empts the first-line reviewer objection to the v11 cross-layer
discordance: *"the RNA→protein mismatch is just proteome detectability."*

Three audit deliverables:

  1. **Detectability gradient** — per-gene observability summary (n
     channels quantified out of total, missing fraction, peptide count,
     mean abundance) collapsed from ``protein_effects_by_row.tsv``,
     then binned by RNA-effect magnitude.  Shows whether large-RNA-
     effect genes are systematically missing from the proteome.

  2. **Observability-matched re-test** — extends the Module 2 matched-
     null stratum with a per-gene missing-fraction bin and re-runs the
     per-pathway propagation tests.  The Module 2 q-value (matched on
     abundance × peptide) and the Module 3 q-value (matched on
     abundance × peptide × missing-fraction) are reported side by side.
     If a pathway's effect drops below significance only after adding
     missingness to the null, the central claim about that pathway must
     be softened to "decoupling beyond what detectability explains is
     modest."  This is the calibrated reviewer-2 check.

  3. **High-coverage subset re-test** — restricts the pool to high-
     confidence genes (``n_peptides >= N_min`` and ``missing_fraction
     <= eps``) and re-runs the propagation tests on that subset.
     Confirms the headline layer assignments (ecm_inverted,
     dct_phospho_carry) persist when detectability uncertainty is at
     a minimum.

  4. **NCC/SPAK phosphosite observability check** — for each
     pre-specified NCC regulatory phosphosite (Slc12a3 S53/S58/S65/S68,
     Stk39 S382/S383), report observability metrics (n_fl, n_gc,
     missing-fraction percentile vs. the phosphoproteome).  Confirms
     the suppressed regulatory sites are NOT in the low-observability
     tail.

The matched-null engine is shared with Module 2 via
:mod:`src.v11.matched_null`; the per-gene per-pathway statistics are
shared with Module 2 via :mod:`src.v11.rna_protein_propagation`.
"""
```

### Comment block, original lines 59-62

```python
# ─────────────────────────────────────────────────────────────────────────────
# Per-gene observability features
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 162-165

```python
# ─────────────────────────────────────────────────────────────────────────────
# Detectability gradient (RNA effect bin → fraction quantified)
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 198-201

```python
# ─────────────────────────────────────────────────────────────────────────────
# High-coverage subset
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 226-229

```python
# ─────────────────────────────────────────────────────────────────────────────
# Per-pathway propagation re-test with a custom stratum column
# ─────────────────────────────────────────────────────────────────────────────
```

### Comment block, original lines 303-306

```python
# ─────────────────────────────────────────────────────────────────────────────
# NCC / SPAK phosphosite observability
# ─────────────────────────────────────────────────────────────────────────────
```

## `src/v11/perturbation_gse228367_lowk.py`

### Module docstring, original lines 2-15

```python
"""GSE228367 low-potassium DCT perturbation alignment.

This module asks whether the spaceflight kidney RNA vectors resemble, oppose,
or are unrelated to a native DCT-enriched low-potassium response.  The GEO raw
archive contains three normal-potassium (NK) and three potassium-depleted (KD)
filtered 10x matrices.  We compute sample-level log-normalized pseudobulk means
from those matrices, estimate KD - NK gene effects, and compare that effect
vector with spaceflight RNA effect vectors.

The low-K response is DCT-enriched pseudobulk, not a DCT1/DCT2 cell-isolated
contrast.  DCT1/DCT2 specificity is therefore evaluated by restricting genes to
the native DCT1-core and DCT2-core reference priors already built from
GSE228367 control Seurat objects.
"""
```

## `src/v11/publication_figures.py`

### Module docstring, original lines 2-10

```python
"""Publication-ready v11 figures.

The v11 figure labels use manuscript-facing terminology:
  * DCT-prior panels show enrichment among distal-nephron subtype-prior parent
    genes in whole-kidney OSD-462 phosphoproteomics.
  * Composition-aware panels separate robust top-decile enrichment from weak
    continuous DCT1-gradient models.
  * Spatial panels are external IRI reference contextualization.
"""
```

## `src/v11/recurrence_meta.py`

### Module docstring, original lines 1-28

```python
"""Repair B -- cross-cohort recurrence meta-analysis.

Turns the descriptive cosine recurrence (0.87 / 0.64 / -0.51 across three
cohorts) into a precision-weighted random-effects meta-analysis with per-gene
FDR and heterogeneity, across *five* on-disk mouse-kidney spaceflight cohorts
(OSD-102, OSD-163, OSD-253, OSD-462, OSD-513; two strains).

The contract the manuscript leans on is *sign-faithful, precision-weighted
pooling*: a gene coherently suppressed across cohorts must yield a negative
pooled effect with a small meta-p and (when the cohorts agree) a low I^2; the
per-cohort effects feeding the pool are honest flight-vs-ground contrasts on the
VST scale, computed in Python (Welch t; no R). The DCT/NCC-WNK transport program
and the ECM/matrix program are then re-scored with the *pooled* statistic in
place of cosine, and leave-one-cohort-out reports stability.

Design provenance
-----------------
Per-cohort flight/ground membership is read from the GeneLab VST column names:
``_FLT_`` -> flight, ``_GC_`` -> ground. Vivarium (``VIV``), basal (``BSL``) and
the OSD-253 ground re-run batch (``GCrerun``) are excluded so every cohort
contributes a clean hardware-ground-control contrast. OSD-462 uses the totRNA
library (rRNA-depleted total RNA), matching the rRNA-removed library type of the
other four cohorts.

Existing ``contrast_vectors/cross_osdr_recurrence/`` artifacts are PROGENy
*pathway*-level bootstraps, not per-gene effect+SE, so they cannot substitute for
the per-gene meta inputs; per-cohort effects are computed from VST here.
"""
```

### Comment block, original lines 45-48

```python
# --------------------------------------------------------------------------- #
# Configuration / on-disk locations
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 66-68

```python
#: OSD-163 is the BAL-TAL strain; the rest are C57-derived. Recorded for the
#: heterogeneity narrative (a strain difference *adds* robustness if the signal
#: survives it).
```

### Comment block, original lines 81-84

```python
# --------------------------------------------------------------------------- #
# Design resolution
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 104-107

```python
# --------------------------------------------------------------------------- #
# (1) per-cohort effect + SE  (pure)
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 175-178

```python
# --------------------------------------------------------------------------- #
# (2) random-effects meta-analysis  (pure)
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 255-258

```python
# --------------------------------------------------------------------------- #
# Gene-set scoring with the pooled statistic
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 309-312

```python
# --------------------------------------------------------------------------- #
# (3) leave-one-cohort-out  (pure)
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 343-346

```python
# --------------------------------------------------------------------------- #
# Orchestrator
# --------------------------------------------------------------------------- #
```

## `src/v11/spatial_reference_projection.py`

### Module docstring, original lines 2-15

```python
"""External kidney injury/repair spatial reference projection for v11.

Primary spatial source:
  * GSE269622 Visium, whole-transcriptome Space Ranger archives.

Secondary spatial source:
  * GSE269719 processed Xenium AnnData, used only for annotation/neighborhood
    inventory because the panel is targeted and not a genome-wide projection
    source.

This script does not localize the RR-10 spaceflight lesion. It asks which IRI
timepoints or marker-enriched Visium spatial niches most resemble the bulk
spaceflight RNA remodeling vector.
"""
```

## `src/validation/__init__.py`

### Module docstring, original lines 1-9

```python
"""
RRRM-2 Validation Module

Leakage-safe cross-validation and sample-level feature extraction.

This module contains:
- cross_validation: Leakage-safe CV framework
- sample_features: Mouse-level feature extraction for LIONESS networks
"""
```

## `src/validation/cross_validation.py`

### Comment block, original lines 40-43

```python
# ---------------------------------------------------------------------------
# Fold-wise skeleton + LIONESS (leakage-safe)
# ---------------------------------------------------------------------------
```

### Comment block, original lines 163-166

```python
# ---------------------------------------------------------------------------
# Main CV loop
# ---------------------------------------------------------------------------
```

### Comment block, original lines 275-277

```python
        # ── 3c) Compute LIONESS for test samples ─────────────
        # For test: compute relative to training pool
        # Each test sample is added to the pool individually
```

## `src/validation/external_replication.py`

### Comment block, original lines 273-275

```python
    # `same_direction` is object-dtype because non-directional rows hold pd.NA.
    # In pandas >= 2.2, .fillna on object dtype triggers a FutureWarning about
    # silent downcasting. Build the boolean column from raw values to avoid it.
```

## `src/validation/osd_external_validation.py`

### Module docstring, original lines 1-8

```python
"""
Protocol-bound external cohort validation and context mapping.

This module builds compact external feature tables from the downloaded GeneLab
processed VST matrices and immediately evaluates them through
src.validation.external_replication. It is intentionally separate from
multi-study pooling: each OSD cohort is analyzed independently.
"""
```

## `src/validation/sample_features.py`

### Module docstring, original lines 1-11

```python
"""
Sample-Level Feature Extraction for RRRM-2

Extract mouse-level features from LIONESS networks for supervised learning.

Per methodology Section 5:
    - Mean/median edge weight within pre-registered pathway subnetworks
    - Node strength (sum of incident weights) for candidate genes
    - Shifter-centered connectivity scores
    - Low-dimensional summaries (top PCs of edge weights)
"""
```

## `src/validation/wgcna_external_validation.py`

### Module docstring, original lines 1-12

```python
"""
WGCNA module-projection external validation.

Projects RRRM-2 WGCNA module gene sets into external OSD cohorts and tests
FLT vs GC module score shifts via permutation. Does NOT rebuild WGCNA in
external cohorts (they're too small).

For each external cohort × module:
  1. Compute module score = mean z-scored expression of module genes present
  2. Welch t-test (FLT vs GC) + label-permutation p-value
  3. BH correction across modules within each study
"""
```

## `src/visualization/__init__.py`

### Module docstring, original lines 1-9

```python
"""
RRRM-2 Visualization Module

Publication-ready plots and network diagnostics.

This module contains:
- publication_plots: Generate publication-ready plots for all phases
- network_diagnostics: Network skeleton visualization and diagnostics
"""
```

## `src/visualization/publication_figures.py`

### Module docstring, original lines 2-21

```python
"""
publication_figures.py — Generate all key publication-ready figures.

Figures produced (all saved under <results_dir>/figures/publication/):
  1. rewiring_vs_logfc.png        — 4-panel scatter: rewiring vs |log2FC|, annotated
  2. edge_sum_perm_rank.png       — Per-contrast edge-sum rank plots with candidates labeled
  3. pathway_dotplot.png          — OR dot plot across all contrasts × pathways
  4. pathway_barplot.png          — Horizontal bar chart of top pathway hits
  5. age_comparison.png           — ISS-T Young vs ISS-T Old rewiring scatter
  6. arm_comparison.png           — ISS-T vs LAR recovery (persistence) plot
  7. pipeline_schematic.png       — Text-based pipeline overview
  8. retinol_subnetwork.png       — Retinol gene rewiring across contrasts
  9. figures/v11/*.png/.pdf       — V11 DCT1/phosphoproteome/mediation figures, when present
 10. result_tables.tsv            — Table 1 (sig genes) + Table 2 (sig pathways)

Usage (standalone):
    python -m src.visualization.publication_figures --results_dir data/results/run_XYZ

Called automatically by Phase 9 of the pipeline.
"""
```

## `src/visualization/publication_plots.py`

### Module docstring, original lines 2-9

```python
"""
plot_output.py - Generate publication-ready plots and tables for each pipeline phase.

Organizes outputs into a clean folder structure mirroring the export script's phase organization.

Usage:
    python scripts/plot_output.py --out_root plots --tag paper_figures
"""
```

## `src/visualization/wgcna_publication_figures.py`

### Module docstring, original lines 2-15

```python
"""
WGCNA Publication Figures
=========================
Generates manuscript-ready figures for the WGCNA analysis when invoked via
--network-method wgcna in the pipeline.

Figures:
  1. Module-trait heatmap (Flight/Age/Arm correlations)
  2. Grey60 eigengene dot plot (by EnvGroup × Arm × Age)
  3. Module preservation bar chart
  4. Hub gene kME lollipop chart for grey60
  5. External validation cross-study comparison
  6. Simple-effect contrast forest plot
"""
```

## `tests/test_aldosterone_axis.py`

### Module docstring, original lines 1-7

```python
"""Unit tests for Reach F -- the aldosterone / mineralocorticoid-axis test.

The contract: the directional panel score is sign-faithful (a coherently
suppressed aldosterone program yields a negative axis effect with a small
permutation p), the meta layer reuses signed pooling and a sign-test honestly,
and the permutation null is competitive (random panels are not significant).
"""
```

## `tests/test_dct_continuous_gradient.py`

### Module docstring, original lines 1-8

```python
"""Unit tests for Repair C -- the continuous DCT1<->DCT2 phospho gradient.

These exercise the pure estimators on planted data. The contract the
manuscript leans on is sign-faithful slope recovery: a planted monotone
gradient must yield a positive slope/rho with a small p, and symmetric /
no-gradient data must yield a slope ~ 0 that is not significant. The spline
test must additionally catch curvature that the linear contrast misses.
"""
```

## `tests/test_human_concordance.py`

### Comment block, original lines 79-81

```python
    # AGT's predicted direction is read off the same Fig. S4A panel that supplies its
    # observation, so it cannot be discordant by construction -> it stays directionally
    # concordant but is no longer counted in the sign test.
```

### Comment block, original lines 122-124

```python
    # Sodium, 24 h volume, magnesium, and figure-level AQP2 are scored; AGT is now
    # report-only (circular) and potassium is context, so four scored analytes
    # collapse to three independent physiological axes.
```

## `tests/test_kinome_atlas_ksea.py`

### Module docstring, original lines 1-11

```python
"""Unit tests for Repair A -- kinome-wide atlas KSEA.

The contract: motif scoring is consensus-faithful (a kinase scores its own
consensus motif above a foreign one); percentile assignment + ``build_kinome_net``
hand the existing ``ksea`` a clean (kinase, gene, site) net; and a kinase whose
substrate sites are planted in the *down* set comes back with negative KSEA z and
is flagged by the parent-gene-aware over-representation Fisher -- without ever
being hand-curated into the net. SPAK (atlas label ``STLK3``) and OSR1 are the
real positive control; here we plant a synthetic basophilic kinase to prove the
machinery recovers a suppressed kinase unbiasedly.
"""
```

### Comment block, original lines 33-36

```python
# --------------------------------------------------------------------------- #
# fixtures: a two-kinase synthetic atlas + planted sites
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 83-86

```python
# --------------------------------------------------------------------------- #
# motif scoring
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 103-106

```python
# --------------------------------------------------------------------------- #
# assignment + net
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 122-125

```python
# --------------------------------------------------------------------------- #
# KSEA on the derived net
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 146-149

```python
# --------------------------------------------------------------------------- #
# over-representation cross-check
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 166-169

```python
# --------------------------------------------------------------------------- #
# real-atlas smoke test (skips if the 1.6 MB workbook isn't staged)
# --------------------------------------------------------------------------- #
```

## `tests/test_osd462_anchor.py`

### Module docstring, original lines 1-6

```python
"""Unit tests for the OSD-462 multi-omics anchor module.

These tests use small synthetic TMT tables so they do not depend on the large
workbooks; they validate the flight-effect arithmetic, plex correction,
gene collapse, the matched-null sampler, and the pathway-cosine geometry.
"""
```

## `tests/test_recurrence_meta.py`

### Module docstring, original lines 1-8

```python
"""Unit tests for Repair B -- the cross-cohort recurrence meta-analysis.

The contract: per-cohort effects are sign-faithful Welch contrasts on the VST
scale; the DerSimonian-Laird pool recovers a planted shared effect with I^2 ~ 0
when cohorts agree and high I^2 when they disagree; BH-FDR separates a truly
shared signal from noise; and gene-set re-scoring + leave-one-cohort-out are
sign-faithful and stable.
"""
```

### Comment block, original lines 21-24

```python
# --------------------------------------------------------------------------- #
# design resolution
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 41-44

```python
# --------------------------------------------------------------------------- #
# per-cohort effect + SE
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 78-81

```python
# --------------------------------------------------------------------------- #
# random-effects meta
# --------------------------------------------------------------------------- #
```

### Comment block, original lines 137-140

```python
# --------------------------------------------------------------------------- #
# gene-set scoring + leave-one-out
# --------------------------------------------------------------------------- #
```

## `tests/test_regulator_activity.py`

### Module docstring, original lines 1-7

```python
"""Unit tests for the v10 regulator-activity layer.

These tests exercise the offline machinery (KSEA statistic, recurrence
classification, evidence grading). The decoupler TF/pathway inference is not
unit-tested here because it requires network-fetched priors; its wrapper is
covered by an interface smoke check only.
"""
```

