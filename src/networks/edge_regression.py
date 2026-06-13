# src/networks/edge_regression.py
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
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd

from src.common import find_sample_col, normalize_labels


def validate_covariate_policy(expression_source: str, requested_covariates: list[str]) -> None:
    """Guard against double-adjusting residualized expression-derived edge weights."""
    if expression_source == "residualized" and requested_covariates:
        raise ValueError(
            "--add_covariates is not allowed with --expression-source=residualized. "
            "Use the residualized primary model edge_weight ~ Age * Arm * EnvGroup, "
            "or explicitly set --expression-source=non_residualized."
        )




def main():
    ap = argparse.ArgumentParser(description="Phase 2: edge-wise regression (limma) + predicted networks")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz")
    ap.add_argument("--phase2_dir", default="data/processed/networks/phase2")
    ap.add_argument("--edge-weights", "--z", dest="edge_weights", default="lioness_edges.npy",
                    help="LIONESS edge-weight matrix filename under phase2_dir")
    ap.add_argument("--outdir", default="data/processed/networks/phase2/regression")
    ap.add_argument("--expression-source", choices=["residualized", "non_residualized"],
                    default="residualized",
                    help="Use residualized for Rtech inputs after batch/deconv/SV removal. "
                         "Covariates are gated off in this mode to avoid double adjustment.")
    ap.add_argument("--add_covariates", default="",
                    help="Comma-separated nuisance covariates. Allowed only with --expression-source=non_residualized.")
    ap.add_argument("--pool-controls", action="store_true",
                    help="Pool GC+VIV+BSL as ground reference (increases control n from 5 to 15)")
    args = ap.parse_args()
    pool = args.pool_controls

    phase2 = Path(args.phase2_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 2 Step A2-A3: Edge-wise Regression")
    print("=" * 60)

    # Load LIONESS edge weights (N x E) and sample order
    weight_path = phase2 / args.edge_weights
    if not weight_path.exists() and args.edge_weights == "lioness_edges.npy":
        legacy = phase2 / "lioness_z_edges.npy"
        if legacy.exists():
            weight_path = legacy
            print("  WARNING: using legacy lioness_z_edges.npy. Prefer lioness_edges.npy for corrected defaults.")
    print(f"\nLoading LIONESS edge weights: {weight_path}")
    Z = np.load(weight_path)
    N, E = Z.shape
    print(f"  Shape: {N} samples × {E} edges")

    # CRITICAL: Load sample order from LIONESS script (safe alignment)
    samples_path = phase2 / "lioness_samples.txt"
    if not samples_path.exists():
        raise FileNotFoundError(
            f"Missing {samples_path}. Save sample order when computing LIONESS."
        )
    samples = [s.strip() for s in samples_path.read_text().splitlines() if s.strip()]
    if len(samples) != N:
        raise ValueError(f"lioness_samples.txt has {len(samples)} samples but Z has {N} rows")
    print(f"  Loaded sample order: {len(samples)} samples")

    # Load and align metadata
    print(f"\nLoading metadata: {args.meta}")
    meta = pd.read_csv(args.meta, sep="\t", compression="gzip")
    sample_col = find_sample_col(meta)
    meta = meta.set_index(sample_col, drop=False)

    # Align meta to LIONESS sample order
    missing = [s for s in samples if s not in meta.index]
    if missing:
        raise ValueError(f"{len(missing)} samples in lioness_samples.txt not in meta. First 5: {missing[:5]}")
    meta = meta.loc[samples].copy()
    print(f"  Aligned: {len(meta)} samples")

    # Normalize labels
    meta = normalize_labels(meta)

    # Factor encoding
    print("\nFactor encoding:")
    for col, allowed in [("Age", ["YNG", "OLD"]), 
                         ("Arm", ["ISS-T", "LAR"]), 
                         ("EnvGroup", ["BSL", "FLT", "GC", "VIV"])]:
        if col not in meta.columns:
            raise KeyError(f"Missing {col} in metadata")
        actual = set(meta[col].astype(str).unique())
        expected = set(allowed)
        if not actual.issubset(expected):
            print(f"  WARNING: {col} has values {actual - expected} not in expected {expected}")
        meta[col] = pd.Categorical(meta[col].astype(str), categories=allowed, ordered=False)
        print(f"  {col}: {list(meta[col].unique())}")

    # Covariates (numeric vs categorical). Residualized expression is already
    covs = [c.strip() for c in args.add_covariates.split(",") if c.strip()]
    validate_covariate_policy(args.expression_source, covs)
    covs = [c for c in covs if c in meta.columns]
    if args.expression_source == "non_residualized" and not covs:
        print("\n  WARNING: non_residualized mode selected without nuisance covariates.")
    print(f"\nExpression source: {args.expression_source}")
    print(f"Covariates found: {covs}")
    
    cov_terms = []
    for c in covs:
        if pd.api.types.is_numeric_dtype(meta[c]):
            cov_terms.append(f"scale(`{c}`)")
            print(f"  {c}: numeric (will scale)")
        else:
            meta[c] = meta[c].astype(str)
            n_levels = meta[c].nunique()
            if n_levels < 2:
                print(f"  {c}: categorical (SKIPPING - only {n_levels} level)")
                continue
            cov_terms.append(f"as.factor(`{c}`)")
            print(f"  {c}: categorical ({n_levels} levels)")
    cov_str = " + " + " + ".join(cov_terms) if cov_terms else ""

    # R setup via rpy2
    print("\nInitializing R/limma...")
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        from rpy2.robjects.conversion import localconverter
        from rpy2.robjects.packages import importr
        # Create converter for pandas <-> R conversions
        pandas_converter = ro.default_converter + pandas2ri.converter
    except Exception as e:
        raise ImportError(
            "Need rpy2 installed for limma regression. "
            "Install rpy2 + R + limma, or ask for the pure-Python fallback."
        ) from e

    limma = importr("limma")
    base = importr("base")
    stats = importr("stats")

    # Build design matrix with explicit cell factor (drops empty combinations)
    print("\nBuilding design matrix...")
    
    # Put meta into R (using converter context)
    with localconverter(pandas_converter):
        r_meta = pandas2ri.py2rpy(meta.reset_index(drop=True))

    ro.globalenv["meta_df"] = r_meta
    
    # Create a single "cell" factor with drop=TRUE to exclude empty combinations
    # This avoids rank-deficiency from missing Age×Arm×EnvGroup combinations
    ro.r("meta_df$cell <- interaction(meta_df$Age, meta_df$Arm, meta_df$EnvGroup, drop=TRUE)")
    
    # Show which cells are present
    cell_levels = list(ro.r("levels(meta_df$cell)"))
    print(f"  Cell factor levels ({len(cell_levels)}): {cell_levels}")
    
    # Formula: 0 + cell (cell means model) + gated covariates
    formula = f"~ 0 + cell{cov_str}"
    print(f"  Formula: {formula}")
    
    ro.r(f"formula_obj <- as.formula('{formula}')")
    ro.r("design <- model.matrix(formula_obj, data=meta_df)")
    
    # Sanitize column names for makeContrasts (R requires syntactically valid names)
    ro.r("colnames(design) <- make.names(colnames(design))")

    # Get column names (now sanitized)
    cn = list(ro.r("colnames(design)"))
    print(f"  Design columns ({len(cn)}): {cn[:8]}{'...' if len(cn) > 8 else ''}")
    
    # Save design column names for reproducibility
    (outdir / "design_colnames.txt").write_text("\n".join(cn) + "\n")

    # Transfer Z to R (safe column-major fill)
    print("\nTransferring LIONESS matrix to R...")
    # Safer approach: fill nrow=N, ncol=E (matches Z shape), then transpose in R
    # This guarantees Y[e, s] == Z[s, e] after transpose
    vec = ro.FloatVector(np.asarray(Z, dtype=np.float64).reshape(-1, order="F"))
    ro.globalenv["Y"] = ro.r.matrix(vec, nrow=N, ncol=E)
    ro.r("Y <- t(Y)")  # Now E x N
    
    # Add rownames for edge mapping (critical for topTable ordering)
    ro.r(f"rownames(Y) <- paste0('e', seq_len({E}))")
    print(f"  Y matrix: {E} edges × {N} samples (with rownames)")

    # Fit limma
    print("\nFitting limma model...")
    ro.r("fit <- lmFit(Y, design)")
    ro.r("fit <- eBayes(fit)")
    print("  lmFit + eBayes complete")

    # Define contrasts
    def term(age, arm, env):
        # interaction() pattern: Age.Arm.EnvGroup (dots as separators)
        raw_term = f"{age}.{arm}.{env}"
        # make.names() turns "-" into "."
        safe_term = raw_term.replace("-", ".")
        # Full coefficient name has "cell" prefix
        raw_coef = f"cell{raw_term}"
        safe_coef = f"cell{safe_term}"
        
        if raw_coef in cn:
            return raw_coef
        if safe_coef in cn:
            return safe_coef
        raise RuntimeError(f"Could not find coefficient for {raw_coef} (or {safe_coef}) in design. Available: {cn[:10]}...")

    # Backtick wrapper for bulletproof contrast expressions in R
    def bt(x):
        return f"`{x}`"

    # Validate all needed terms exist
    needed_combos = [
        ("YNG", "ISS-T", "FLT"), ("YNG", "ISS-T", "GC"),
        ("OLD", "ISS-T", "FLT"), ("OLD", "ISS-T", "GC"),
        ("YNG", "LAR", "FLT"),   ("YNG", "LAR", "GC"),
        ("OLD", "LAR", "FLT"),   ("OLD", "LAR", "GC"),
    ]
    for age, arm, env in needed_combos:
        term(age, arm, env)  # Will raise if not found

    # Helper: build the "ground" term — either GC alone or (GC+VIV+BSL)/3
    def gnd(age, arm):
        if pool:
            gc  = bt(term(age, arm, 'GC'))
            viv = bt(term(age, arm, 'VIV'))
            bsl = bt(term(age, arm, 'BSL'))
            return f"({gc} + {viv} + {bsl}) / 3"
        else:
            return bt(term(age, arm, 'GC'))

    # Suffix for contrast names
    ctrl_tag = "GND" if pool else "GC"
    if pool:
        # Validate VIV+BSL terms also exist
        for age, arm, env in [("YNG", "ISS-T", "VIV"), ("YNG", "ISS-T", "BSL"),
                              ("OLD", "ISS-T", "VIV"), ("OLD", "ISS-T", "BSL"),
                              ("YNG", "LAR", "VIV"),   ("YNG", "LAR", "BSL"),
                              ("OLD", "LAR", "VIV"),   ("OLD", "LAR", "BSL")]:
            term(age, arm, env)
        print(f"  Pool controls: FLT vs (GC+VIV+BSL)/3 — control n tripled")

    # Contrast effects (Δz for rewiring) - use backticks for safe R parsing
    contrasts = {
        f"ISS_T_YNG_FLT_minus_{ctrl_tag}": f"{bt(term('YNG','ISS-T','FLT'))} - {gnd('YNG','ISS-T')}",
        f"ISS_T_OLD_FLT_minus_{ctrl_tag}": f"{bt(term('OLD','ISS-T','FLT'))} - {gnd('OLD','ISS-T')}",
        f"LAR_YNG_FLT_minus_{ctrl_tag}":   f"{bt(term('YNG','LAR','FLT'))} - {gnd('YNG','LAR')}",
        f"LAR_OLD_FLT_minus_{ctrl_tag}":   f"{bt(term('OLD','LAR','FLT'))} - {gnd('OLD','LAR')}",
        f"ISS_T_AgeDep_Flight":    (f"({bt(term('OLD','ISS-T','FLT'))} - {gnd('OLD','ISS-T')}) - "
                                   f"({bt(term('YNG','ISS-T','FLT'))} - {gnd('YNG','ISS-T')})"),
        f"LAR_AgeDep_Flight":      (f"({bt(term('OLD','LAR','FLT'))} - {gnd('OLD','LAR')}) - "
                                   f"({bt(term('YNG','LAR','FLT'))} - {gnd('YNG','LAR')})"),
        f"ISS_minus_LAR_YNG_Flight": (f"({bt(term('YNG','ISS-T','FLT'))} - {gnd('YNG','ISS-T')}) - "
                                     f"({bt(term('YNG','LAR','FLT'))} - {gnd('YNG','LAR')})"),
        f"ISS_minus_LAR_OLD_Flight": (f"({bt(term('OLD','ISS-T','FLT'))} - {gnd('OLD','ISS-T')}) - "
                                     f"({bt(term('OLD','LAR','FLT'))} - {gnd('OLD','LAR')})"),
    }

    # Predicted condition-specific networks (z_hat per cell) - will extract from coefficients
    preds = {
        "Pred_YNG_ISS_T_FLT": term("YNG", "ISS-T", "FLT"),
        "Pred_YNG_ISS_T_GC":  term("YNG", "ISS-T", "GC"),
        "Pred_OLD_ISS_T_FLT": term("OLD", "ISS-T", "FLT"),
        "Pred_OLD_ISS_T_GC":  term("OLD", "ISS-T", "GC"),
        "Pred_YNG_LAR_FLT":   term("YNG", "LAR", "FLT"),
        "Pred_YNG_LAR_GC":    term("YNG", "LAR", "GC"),
        "Pred_OLD_LAR_FLT":   term("OLD", "LAR", "FLT"),
        "Pred_OLD_LAR_GC":    term("OLD", "LAR", "GC"),
    }

    # Save contrast definitions for reproducibility
    all_contrasts = {"effects": contrasts, "predictions": preds}
    with open(outdir / "contrasts.json", "w") as f:
        json.dump(all_contrasts, f, indent=2)

    # Run contrast effects (Δz for rewiring)
    print("\nRunning contrast effects...")
    
    for name, con in contrasts.items():
        ro.r(f"cm <- makeContrasts(`{name}` = {con}, levels=design)")
        ro.r("fit2 <- contrasts.fit(fit, cm)")
        ro.r("fit2 <- eBayes(fit2)")
        
        # Save edge-weight contrast from fit2$coefficients (legacy filename
        # keeps downstream compatibility).
        coef = np.array(ro.r("fit2$coefficients"))[:, 0].astype(np.float32)  # E-length
        np.save(outdir / f"{name}_delta_z.npy", coef)
        np.save(outdir / f"{name}_delta_edge_weight.npy", coef)
        
        # Save topTable for p-values etc
        ro.r("tt <- topTable(fit2, number=Inf, sort.by='none')")
        with localconverter(pandas_converter):
            tt = ro.conversion.get_conversion().rpy2py(ro.globalenv["tt"])
        
        # Guard: verify topTable row count matches expected edges
        if tt.shape[0] != E:
            raise RuntimeError(f"topTable returned {tt.shape[0]} rows, expected {E}.")
        
        tt = tt.reset_index(drop=True)
        edge_index = pd.read_csv(phase2 / "edge_index.tsv", sep="\t")
        if len(edge_index) == len(tt):
            provenance_cols = [
                c for c in ["gene_i", "gene_j", "edge_origin", "edge_source_detail", "is_fixed_prior"]
                if c in edge_index.columns
            ]
            tt = pd.concat([edge_index[provenance_cols].reset_index(drop=True), tt], axis=1)
            if "is_fixed_prior" in tt.columns:
                fixed = tt["is_fixed_prior"].fillna(False).astype(bool)
            else:
                fixed = pd.Series(False, index=tt.index)
            tt["edge_p_value_scope"] = np.where(
                fixed,
                "fixed_prior_edge_inference",
                "exploratory_post_selection_data_driven_skeleton",
            )
        tt.to_csv(outdir / f"{name}_limma.tsv", sep="\t", index=False)

        print(f"  {name}: delta_z saved ({coef.min():.3f} to {coef.max():.3f})")

    # Extract predicted condition networks directly from fit$coefficients
    print("\nExtracting predicted networks from fit$coefficients...")
    
    # Get coefficient matrix (E x ncol)
    coef_matrix = np.array(ro.r("fit$coefficients"))
    coef_names = list(ro.r("colnames(fit$coefficients)"))
    
    for name, coef_term in preds.items():
        if coef_term not in coef_names:
            print(f"  WARNING: {coef_term} not in coefficients, skipping {name}")
            continue
        
        coef_idx = coef_names.index(coef_term)
        z_hat = coef_matrix[:, coef_idx].astype(np.float32)
        np.save(outdir / f"{name}_z_hat.npy", z_hat)
        np.save(outdir / f"{name}_edge_weight_hat.npy", z_hat)
        print(f"  {name}: z_hat saved ({z_hat.min():.3f} to {z_hat.max():.3f})")

    # If pooling controls, also create averaged GND predicted networks
    if pool:
        print("\nCreating pooled ground (GND) predicted networks...")
        for age_lbl, arm_lbl, arm_tag in [
            ("YNG", "ISS-T", "ISS_T"), ("OLD", "ISS-T", "ISS_T"),
            ("YNG", "LAR", "LAR"),     ("OLD", "LAR", "LAR"),
        ]:
            gc_term  = term(age_lbl, arm_lbl, "GC")
            viv_term = term(age_lbl, arm_lbl, "VIV")
            bsl_term = term(age_lbl, arm_lbl, "BSL")
            
            gc_idx  = coef_names.index(gc_term)  if gc_term  in coef_names else None
            viv_idx = coef_names.index(viv_term) if viv_term in coef_names else None
            bsl_idx = coef_names.index(bsl_term) if bsl_term in coef_names else None
            
            if gc_idx is not None and viv_idx is not None and bsl_idx is not None:
                z_gnd = ((coef_matrix[:, gc_idx] + coef_matrix[:, viv_idx] + 
                          coef_matrix[:, bsl_idx]) / 3.0).astype(np.float32)
                name = f"Pred_{age_lbl}_{arm_tag}_GND"
                np.save(outdir / f"{name}_z_hat.npy", z_gnd)
                np.save(outdir / f"{name}_edge_weight_hat.npy", z_gnd)
                print(f"  {name}: z_hat saved ({z_gnd.min():.3f} to {z_gnd.max():.3f})")
            else:
                print(f"  WARNING: missing coefficient for {age_lbl}/{arm_lbl} GND average")

    # Copy edge index for downstream graph building
    edge_index = pd.read_csv(phase2 / "edge_index.tsv", sep="\t")
    edge_index.to_csv(outdir / "edge_index.tsv", sep="\t", index=False)
    model_metadata = {
        "edge_weight_file": str(weight_path),
        "expression_source": args.expression_source,
        "primary_model": "edge_weight ~ 0 + Age:Arm:EnvGroup cell means",
        "nuisance_covariates": covs,
        "double_adjustment_guard": args.expression_source == "residualized",
        "edge_p_values_on_data_driven_skeleton": "exploratory",
    }
    (outdir / "edge_regression_model_metadata.json").write_text(json.dumps(model_metadata, indent=2) + "\n")

    print(f"\n{'=' * 60}")
    print("Edge-wise regression complete")
    print(f"{'=' * 60}")
    print(f"\nOutputs in {outdir}:")
    print(f"  - design_colnames.txt")
    print(f"  - contrasts.json")
    print(f"  - *_limma.tsv (topTable per contrast effect)")
    print(f"  - *_delta_z.npy / *_delta_edge_weight.npy (effect sizes for rewiring)")
    print(f"  - Pred_*_z_hat.npy / Pred_*_edge_weight_hat.npy (predicted condition networks)")
    print(f"  - edge_index.tsv")
    print(f"  - edge_regression_model_metadata.json")


if __name__ == "__main__":
    main()
