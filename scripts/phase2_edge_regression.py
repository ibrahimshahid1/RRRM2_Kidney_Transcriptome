# scripts/phase2_edge_regression.py
"""
Phase 2 Step A2-A3: Edge-wise Regression + Predicted Networks

Fits limma edge-wise regression on LIONESS Fisher-z weights with the full
2×2×4 design (Age × Arm × EnvGroup) + technical covariates.

Outputs:
  - Contrast effects (Δz) for rewiring analysis
  - Predicted condition-specific networks (z_hat)
  - limma topTable results per contrast

Usage:
    python scripts/phase2_edge_regression.py
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd


def find_sample_col(meta: pd.DataFrame) -> str:
    """Find the sample ID column in metadata."""
    for col in ["Sample Name (raw_counts_colname)", "Sample Name", "sample"]:
        if col in meta.columns:
            return col
    return meta.columns[0]


def normalize_labels(meta: pd.DataFrame) -> pd.DataFrame:
    """Normalize factor labels to canonical forms."""
    meta = meta.copy()
    
    # Age normalization
    if "Age" in meta.columns:
        meta["Age"] = meta["Age"].astype(str).replace({
            "Young": "YNG", "Yng": "YNG", "young": "YNG",
            "Old": "OLD", "old": "OLD"
        })
    
    # Arm normalization
    if "Arm" in meta.columns:
        meta["Arm"] = meta["Arm"].astype(str).replace({
            "ISS": "ISS-T", "ISST": "ISS-T", "ISS_T": "ISS-T", "ISS T": "ISS-T",
            "LAR_T": "LAR", "LAR-T": "LAR", "LAR T": "LAR"
        })
    
    # EnvGroup normalization (HGC→GC, VGC→VIV)
    if "EnvGroup" in meta.columns:
        meta["EnvGroup"] = meta["EnvGroup"].astype(str).replace({
            "HGC": "GC", "VGC": "VIV",
            "HGC/GC": "GC", "VIV/VGC": "VIV"
        })
    
    return meta


def main():
    ap = argparse.ArgumentParser(description="Phase 2: edge-wise regression (limma) + predicted networks")
    ap.add_argument("--meta", default="data/processed/phase1_residuals/meta_phase1.tsv.gz")
    ap.add_argument("--phase2_dir", default="data/processed/networks/phase2")
    ap.add_argument("--z", default="lioness_z_edges.npy")
    ap.add_argument("--outdir", default="data/processed/networks/phase2/regression")
    ap.add_argument("--add_covariates", default="LibraryBatch,SeqInstr,ReadDepth,rRNA",
                    help="Comma-separated covariates to include if present")
    args = ap.parse_args()

    phase2 = Path(args.phase2_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Phase 2 Step A2-A3: Edge-wise Regression")
    print("=" * 60)

    # -------------------------------------------------------------------------
    # Load LIONESS z (N x E) and sample order
    # -------------------------------------------------------------------------
    print(f"\nLoading LIONESS z: {phase2 / args.z}")
    Z = np.load(phase2 / args.z)
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

    # -------------------------------------------------------------------------
    # Load and align metadata
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # Factor encoding
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # Covariates (numeric vs categorical)
    # -------------------------------------------------------------------------
    covs = [c.strip() for c in args.add_covariates.split(",") if c.strip()]
    covs = [c for c in covs if c in meta.columns]
    print(f"\nCovariates found: {covs}")
    
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

    # -------------------------------------------------------------------------
    # R setup via rpy2
    # -------------------------------------------------------------------------
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

    # -------------------------------------------------------------------------
    # Build design matrix with explicit cell factor (drops empty combinations)
    # -------------------------------------------------------------------------
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
    
    # Formula: 0 + cell (cell means model) + covariates
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

    # -------------------------------------------------------------------------
    # Transfer Z to R (safe column-major fill)
    # -------------------------------------------------------------------------
    print("\nTransferring LIONESS matrix to R...")
    # Safer approach: fill nrow=N, ncol=E (matches Z shape), then transpose in R
    # This guarantees Y[e, s] == Z[s, e] after transpose
    vec = ro.FloatVector(np.asarray(Z, dtype=np.float64).reshape(-1, order="F"))
    ro.globalenv["Y"] = ro.r.matrix(vec, nrow=N, ncol=E)
    ro.r("Y <- t(Y)")  # Now E x N
    
    # Add rownames for edge mapping (critical for topTable ordering)
    ro.r(f"rownames(Y) <- paste0('e', seq_len({E}))")
    print(f"  Y matrix: {E} edges × {N} samples (with rownames)")

    # -------------------------------------------------------------------------
    # Fit limma
    # -------------------------------------------------------------------------
    print("\nFitting limma model...")
    ro.r("fit <- lmFit(Y, design)")
    ro.r("fit <- eBayes(fit)")
    print("  lmFit + eBayes complete")

    # -------------------------------------------------------------------------
    # Define contrasts
    # -------------------------------------------------------------------------
    
    # Helper to build coefficient names matching interaction(Age, Arm, EnvGroup) pattern
    # R's interaction() produces "Level1.Level2.Level3", then make.names() sanitizes
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

    # Contrast effects (Δz for rewiring) - use backticks for safe R parsing
    contrasts = {
        "ISS_T_YNG_FLT_minus_GC": f"{bt(term('YNG','ISS-T','FLT'))} - {bt(term('YNG','ISS-T','GC'))}",
        "ISS_T_OLD_FLT_minus_GC": f"{bt(term('OLD','ISS-T','FLT'))} - {bt(term('OLD','ISS-T','GC'))}",
        "LAR_YNG_FLT_minus_GC":   f"{bt(term('YNG','LAR','FLT'))} - {bt(term('YNG','LAR','GC'))}",
        "LAR_OLD_FLT_minus_GC":   f"{bt(term('OLD','LAR','FLT'))} - {bt(term('OLD','LAR','GC'))}",
        "ISS_T_AgeDep_Flight":    (f"({bt(term('OLD','ISS-T','FLT'))} - {bt(term('OLD','ISS-T','GC'))}) - "
                                  f"({bt(term('YNG','ISS-T','FLT'))} - {bt(term('YNG','ISS-T','GC'))})"),
        "LAR_AgeDep_Flight":      (f"({bt(term('OLD','LAR','FLT'))} - {bt(term('OLD','LAR','GC'))}) - "
                                  f"({bt(term('YNG','LAR','FLT'))} - {bt(term('YNG','LAR','GC'))})"),
        "ISS_minus_LAR_YNG_Flight": (f"({bt(term('YNG','ISS-T','FLT'))} - {bt(term('YNG','ISS-T','GC'))}) - "
                                    f"({bt(term('YNG','LAR','FLT'))} - {bt(term('YNG','LAR','GC'))})"),
        "ISS_minus_LAR_OLD_Flight": (f"({bt(term('OLD','ISS-T','FLT'))} - {bt(term('OLD','ISS-T','GC'))}) - "
                                    f"({bt(term('OLD','LAR','FLT'))} - {bt(term('OLD','LAR','GC'))})"),
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

    # -------------------------------------------------------------------------
    # Run contrast effects (Δz for rewiring)
    # -------------------------------------------------------------------------
    print("\nRunning contrast effects...")
    
    for name, con in contrasts.items():
        ro.r(f"cm <- makeContrasts(`{name}` = {con}, levels=design)")
        ro.r("fit2 <- contrasts.fit(fit, cm)")
        ro.r("fit2 <- eBayes(fit2)")
        
        # Save Δz from fit2$coefficients (robust to topTable row reordering)
        coef = np.array(ro.r("fit2$coefficients"))[:, 0].astype(np.float32)  # E-length
        np.save(outdir / f"{name}_delta_z.npy", coef)
        
        # Save topTable for p-values etc
        ro.r("tt <- topTable(fit2, number=Inf, sort.by='none')")
        with localconverter(pandas_converter):
            tt = ro.conversion.get_conversion().rpy2py(ro.globalenv["tt"])
        
        # Guard: verify topTable row count matches expected edges
        if tt.shape[0] != E:
            raise RuntimeError(f"topTable returned {tt.shape[0]} rows, expected {E}.")
        
        tt = tt.reset_index(drop=True)
        tt.to_csv(outdir / f"{name}_limma.tsv", sep="\t", index=False)

        print(f"  {name}: delta_z saved ({coef.min():.3f} to {coef.max():.3f})")

    # -------------------------------------------------------------------------
    # Extract predicted condition networks directly from fit$coefficients
    # -------------------------------------------------------------------------
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
        print(f"  {name}: z_hat saved ({z_hat.min():.3f} to {z_hat.max():.3f})")

    # -------------------------------------------------------------------------
    # Copy edge index for downstream graph building
    # -------------------------------------------------------------------------
    edge_index = pd.read_csv(phase2 / "edge_index.tsv", sep="\t")
    edge_index.to_csv(outdir / "edge_index.tsv", sep="\t", index=False)

    print(f"\n{'=' * 60}")
    print("Edge-wise regression complete!")
    print(f"{'=' * 60}")
    print(f"\nOutputs in {outdir}:")
    print(f"  - design_colnames.txt")
    print(f"  - contrasts.json")
    print(f"  - *_limma.tsv (topTable per contrast effect)")
    print(f"  - *_delta_z.npy (effect sizes for rewiring)")
    print(f"  - Pred_*_z_hat.npy (predicted condition networks from coefficients)")
    print(f"  - edge_index.tsv")


if __name__ == "__main__":
    main()
