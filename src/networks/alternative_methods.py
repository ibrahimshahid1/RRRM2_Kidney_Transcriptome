# src/networks/alternative_methods.py
"""
Alternative sample-specific network methods for Phase 8 sensitivity.

Implemented:
  * SSN-style leave-one-sample perturbation of edge correlations.

Guarded wrappers:
  * CSN and scLink are exposed as explicit modes but require external
    implementations/dependencies; the benchmark records them as unavailable
    rather than silently substituting LIONESS.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from src.common import REPO_ROOT
from src.networks.lioness import CLIP_R, ZCAP


def _pearson_edges(X: np.ndarray, edge_i: np.ndarray, edge_j: np.ndarray) -> np.ndarray:
    """Pearson correlations for selected edges on genes × samples matrix."""
    Xi = X[edge_i, :]
    Xj = X[edge_j, :]
    Xi = Xi - Xi.mean(axis=1, keepdims=True)
    Xj = Xj - Xj.mean(axis=1, keepdims=True)
    denom = np.sqrt((Xi ** 2).sum(axis=1) * (Xj ** 2).sum(axis=1))
    r = np.divide((Xi * Xj).sum(axis=1), denom, out=np.zeros(len(edge_i)), where=denom > 1e-12)
    return np.clip(r, -CLIP_R, CLIP_R)


def ssn_delta_z(X: np.ndarray, edge_i: np.ndarray, edge_j: np.ndarray) -> np.ndarray:
    """SSN-style per-sample z-score perturbation against the reference network."""
    n_samples = X.shape[1]
    r_ref = _pearson_edges(X, edge_i, edge_j)
    z_ref = np.arctanh(r_ref)
    out = np.empty((n_samples, len(edge_i)), dtype=np.float32)
    for s in range(n_samples):
        keep = np.ones(n_samples, dtype=bool)
        keep[s] = False
        r_without = _pearson_edges(X[:, keep], edge_i, edge_j)
        z_without = np.arctanh(r_without)
        out[s, :] = np.clip(z_ref - z_without, -ZCAP, ZCAP).astype(np.float32)
    return out


def compute_alternative_network(
    method: str,
    X: np.ndarray,
    edge_i: np.ndarray,
    edge_j: np.ndarray,
) -> np.ndarray:
    method = method.lower()
    if method == "ssn":
        return ssn_delta_z(X, edge_i, edge_j)
    if method == "csn":
        raise ImportError("CSN mode requires an external CSN implementation; not vendored in this repository.")
    if method == "sclink":
        raise ImportError("scLink mode requires scLink/scRNA-specific dependencies; not vendored here.")
    raise ValueError(f"Unknown alternative network method: {method}")


def main() -> None:
    ap = argparse.ArgumentParser(description="Benchmark alternative sample-specific network methods")
    ap.add_argument("--rtech", default=str(REPO_ROOT / "data/processed/phase1_residuals/Rtech.tsv.gz"))
    ap.add_argument("--phase2_dir", default=str(REPO_ROOT / "data/processed/networks/phase2"))
    ap.add_argument("--method", choices=["ssn", "csn", "sclink"], default="ssn")
    ap.add_argument("--outdir", default=str(REPO_ROOT / "data/results/alternative_network_methods"))
    args = ap.parse_args()

    phase2 = Path(args.phase2_dir)
    genes = [g.strip() for g in (phase2 / "phase2_genes.txt").read_text().splitlines() if g.strip()]
    edge_i = np.load(phase2 / "edge_i.npy")
    edge_j = np.load(phase2 / "edge_j.npy")
    rtech = pd.read_csv(args.rtech, sep="\t", compression="gzip", index_col=0).loc[genes]
    X = rtech.values.astype(float)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        weights = compute_alternative_network(args.method, X, edge_i, edge_j)
        np.save(outdir / f"{args.method}_sample_specific_edges.npy", weights)
        status = "ok"
        message = ""
    except Exception as exc:
        status = "unavailable"
        message = str(exc)

    pd.DataFrame([{
        "method": args.method,
        "status": status,
        "message": message,
        "n_genes": len(genes),
        "n_edges": len(edge_i),
    }]).to_csv(outdir / "alternative_methods_benchmark.tsv", sep="\t", index=False)
    print(f"[OK] Alternative method benchmark status written to {outdir}")


if __name__ == "__main__":
    main()
