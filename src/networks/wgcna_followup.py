#!/usr/bin/env python3
"""
WGCNA Manuscript Follow-up: Simple contrasts, hub genes, eigengene plots,
GO/KEGG enrichment, and final integrated table.
"""
from __future__ import annotations
import json, os, warnings, gzip
from pathlib import Path

import numpy as np
from numpy.linalg import inv
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats as sp_stats
import yaml
import gseapy as gp

warnings.filterwarnings("ignore")
# Paths
ROOT  = Path(__file__).resolve().parent.parent.parent
_results_dir = os.environ.get("RRRM_RESULTS_DIR", "")

if _results_dir:
    # Running from pipeline — resolve paths from the active run
    WDIR  = Path(_results_dir) / "wgcna"
    # Find Phase 1 outputs: check current run first, then scan for latest
    _rtech_candidate = Path(_results_dir) / "phase1_residuals/Rtech.tsv.gz"
    if _rtech_candidate.exists():
        RTECH = _rtech_candidate
        META  = Path(_results_dir) / "phase1_residuals/meta_phase1.tsv.gz"
    else:
        # Phase 1 may have been run in a different run directory; scan
        _results_root = ROOT / "data/results"
        for rd in sorted(_results_root.glob("run_*"), reverse=True):
            if (rd / "phase1_residuals/Rtech.tsv.gz").exists():
                RTECH = rd / "phase1_residuals/Rtech.tsv.gz"
                META  = rd / "phase1_residuals/meta_phase1.tsv.gz"
                break
        else:
            raise FileNotFoundError("Cannot find Phase 1 Rtech.tsv.gz in any run directory")
else:
    # Standalone execution — use latest run
    _results_root = ROOT / "data/results"
    for rd in sorted(_results_root.glob("run_*"), reverse=True):
        if (rd / "wgcna/module_assignments.csv").exists():
            WDIR = rd / "wgcna"
            break
    else:
        raise FileNotFoundError("No WGCNA run found")
    for rd in sorted(_results_root.glob("run_*"), reverse=True):
        if (rd / "phase1_residuals/Rtech.tsv.gz").exists():
            RTECH = rd / "phase1_residuals/Rtech.tsv.gz"
            META  = rd / "phase1_residuals/meta_phase1.tsv.gz"
            break
    else:
        raise FileNotFoundError("No Phase 1 Rtech found")

IDMAP = ROOT / "data/processed/resources/id_map.tsv"
GSETS = ROOT / "config/gene_sets.yaml"
OUTDIR = WDIR / "manuscript_followup"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Load data
print("Loading data...")
rtech = pd.read_csv(RTECH, sep="\t", compression="gzip", index_col=0)
meta = pd.read_csv(META, sep="\t", compression="gzip")
scol = [c for c in meta.columns if "Sample Name" in c or c == "sample"][0]
meta = meta.set_index(scol, drop=False)
meta["EnvGroup"] = meta["EnvGroup"].replace({"HGC": "GC", "VGC": "VIV"})
meta["Age"] = meta["Age"].replace({"Young": "YNG", "young": "YNG", "Yng": "YNG", "old": "OLD"})
meta["Arm"] = meta["Arm"].replace({"ISS": "ISS-T", "ISST": "ISS-T", "ISS_T": "ISS-T"})

keep = meta["EnvGroup"].isin(["FLT", "GC"])
meta_fg = meta[keep].copy()
common = [s for s in rtech.columns if s in meta_fg.index]
rtech_fg = rtech[common]
meta_fg = meta_fg.loc[common]

# Module assignments + eigengenes
mod_assign = pd.read_csv(WDIR / "module_assignments.csv")
MEs = pd.read_csv(WDIR / "module_eigengenes.csv", index_col=0)

# ME → color map
me_map = {f"ME{k}": v for k, v in mod_assign.groupby("module_num")["module_color"].first().items()}
color_to_me = {v: k for k, v in me_map.items()}

# ID map
id_df = pd.read_csv(IDMAP, sep="\t", comment="#")
eid_col = [c for c in id_df.columns if "ensembl" in c.lower()][0]
sym_col = [c for c in id_df.columns if "symbol" in c.lower() or "mgi" in c.lower()][0]
sym_map = dict(zip(id_df[eid_col].astype(str), id_df[sym_col].astype(str)))

# Gene set membership for annotation

with open(GSETS) as f:
    gs_config = yaml.safe_load(f)

sym_to_ens = dict(zip(id_df[sym_col].str.lower(), id_df[eid_col].astype(str)))
pathway_genes = {}
for pw_name, pw in gs_config.items():
    if isinstance(pw, dict) and "genes" in pw:
        ens_set = set()
        for s in pw["genes"]:
            e = sym_to_ens.get(s.lower())
            if e:
                ens_set.add(e)
        pathway_genes[pw_name] = ens_set

print(f"  FLT+GC samples: {len(meta_fg)}, Genes: {len(rtech_fg)}")

# Significant modules to analyze
SIG_MODULES = ["grey60", "royalblue", "salmon", "pink", "green", "red"]

print("\n" + "="*60)
print("1. Simple-Effect Contrasts (FLT vs GC per Age×Arm stratum)")
print("="*60)

# Build design
meta_fg["Flight"] = (meta_fg["EnvGroup"] == "FLT").astype(float)
meta_fg["AgeOld"] = (meta_fg["Age"] == "OLD").astype(float)
meta_fg["ArmLAR"] = (meta_fg["Arm"] == "LAR").astype(float)

strata = {
    "ISS-T_YNG": {"AgeOld": 0, "ArmLAR": 0},
    "ISS-T_OLD": {"AgeOld": 1, "ArmLAR": 0},
    "LAR_YNG":   {"AgeOld": 0, "ArmLAR": 1},
    "LAR_OLD":   {"AgeOld": 1, "ArmLAR": 1},
}

# For ME ~ Flight * AgeOld * ArmLAR: Flight effect at (AgeOld=a, ArmLAR=r)

X = np.column_stack([
    np.ones(len(meta_fg)),
    meta_fg["Flight"].values,
    meta_fg["AgeOld"].values,
    meta_fg["ArmLAR"].values,
    meta_fg["Flight"].values * meta_fg["AgeOld"].values,
    meta_fg["Flight"].values * meta_fg["ArmLAR"].values,
    meta_fg["AgeOld"].values * meta_fg["ArmLAR"].values,
    meta_fg["Flight"].values * meta_fg["AgeOld"].values * meta_fg["ArmLAR"].values,
])
XtX_inv = inv(X.T @ X)
term_names = ["(Intercept)", "Flight", "AgeOld", "ArmLAR",
              "Flight:AgeOld", "Flight:ArmLAR", "AgeOld:ArmLAR",
              "Flight:AgeOld:ArmLAR"]

contrast_results = []
for me_name in sorted(me_map.keys(), key=lambda x: int(x[2:])):
    color = me_map[me_name]
    if color == "grey":
        continue
    y = MEs.loc[meta_fg.index, me_name].values
    beta = XtX_inv @ X.T @ y
    resid = y - X @ beta
    n, p = X.shape
    sigma2 = (resid @ resid) / (n - p)

    for stratum_name, levels in strata.items():
        a = levels["AgeOld"]
        r = levels["ArmLAR"]
        c = np.array([0, 1, 0, 0, a, r, 0, a*r])
        est = c @ beta
        se = np.sqrt(sigma2 * c @ XtX_inv @ c)
        t_val = est / se if se > 0 else 0
        p_val = 2 * sp_stats.t.sf(abs(t_val), n - p)
        contrast_results.append({
            "module": me_name, "color": color,
            "stratum": stratum_name,
            "estimate": round(est, 4), "SE": round(se, 4),
            "t": round(t_val, 3), "p": p_val,
        })

cdf = pd.DataFrame(contrast_results)
# BH within each stratum
for s in strata:
    mask = cdf["stratum"] == s
    cdf.loc[mask, "q"] = sp_stats.false_discovery_control(cdf.loc[mask, "p"].values)

cdf.to_csv(OUTDIR / "simple_effect_contrasts.csv", index=False)

# Print significant modules
print(f"\n{'Module':<6} {'Color':<14} {'Stratum':<12} {'Est':>7} {'SE':>6} {'t':>7} {'p':>9} {'q':>9}")
print("-"*80)
for _, r in cdf[cdf["q"] < 0.10].sort_values(["color", "stratum"]).iterrows():
    print(f"{r['module']:<6} {r['color']:<14} {r['stratum']:<12} {r['estimate']:>7.3f} {r['SE']:>6.3f} {r['t']:>7.2f} {r['p']:>9.4f} {r['q']:>9.4f}")

print("\n" + "="*60)
print("2. Grey60 Hub Genes (kME analysis)")
print("="*60)

grey60_genes = mod_assign[mod_assign["module_color"] == "grey60"]["gene"].tolist()
me17 = MEs.loc[meta_fg.index, "ME17"].values
expr_matrix = rtech_fg.loc[[g for g in grey60_genes if g in rtech_fg.index], meta_fg.index]

kme_vals = []
for gene in expr_matrix.index:
    r, _ = sp_stats.pearsonr(expr_matrix.loc[gene].values, me17)
    kme_vals.append({"gene": gene, "kME": round(r, 4), "abs_kME": round(abs(r), 4)})

kme_df = pd.DataFrame(kme_vals).sort_values("abs_kME", ascending=False)
kme_df["symbol"] = kme_df["gene"].map(sym_map).fillna("")
kme_df["rank"] = range(1, len(kme_df) + 1)

# Annotate pathway membership
def annotate_gene(ens_id):
    memberships = []
    for pw_name, gene_set in pathway_genes.items():
        if ens_id in gene_set:
            memberships.append(pw_name)
    return "; ".join(memberships) if memberships else "uncharacterized"

kme_df["pathway_annotation"] = kme_df["gene"].apply(annotate_gene)

kme_df.to_csv(OUTDIR / "grey60_hub_genes_kME.csv", index=False)

print(f"\nGrey60: {len(kme_df)} genes")
print(f"\nTop 20 hub genes by |kME|:")
print(f"{'Rank':>4} {'Symbol':<15} {'kME':>7} {'Ensembl':<20} {'Annotation'}")
print("-"*75)
for _, r in kme_df.head(20).iterrows():
    print(f"{r['rank']:>4} {r['symbol']:<15} {r['kME']:>7.3f} {r['gene']:<20} {r['pathway_annotation']}")

print("\n" + "="*60)
print("3. Eigengene Plots")
print("="*60)

plot_modules = ["grey60", "red", "royalblue", "salmon", "pink", "green"]
plot_dir = OUTDIR / "eigengene_plots"
plot_dir.mkdir(exist_ok=True)

meta_fg["group"] = meta_fg["EnvGroup"] + "\n" + meta_fg["Arm"] + "\n" + meta_fg["Age"]
meta_fg["group_short"] = meta_fg["EnvGroup"] + "_" + meta_fg["Arm"] + "_" + meta_fg["Age"]

group_order = [
    "FLT\nISS-T\nYNG", "GC\nISS-T\nYNG",
    "FLT\nISS-T\nOLD", "GC\nISS-T\nOLD",
    "FLT\nLAR\nYNG", "GC\nLAR\nYNG",
    "FLT\nLAR\nOLD", "GC\nLAR\nOLD",
]
colors_flt = "#E63946"
colors_gc = "#457B9D"

for mod_color in plot_modules:
    me_name = color_to_me.get(mod_color)
    if me_name is None or me_name not in MEs.columns:
        continue

    fig, ax = plt.subplots(figsize=(10, 5))
    me_vals = MEs.loc[meta_fg.index, me_name]
    plot_data = meta_fg.copy()
    plot_data["ME"] = me_vals.values

    # Strip plot with jitter
    positions = []
    for i, grp in enumerate(group_order):
        subset = plot_data[plot_data["group"] == grp]["ME"]
        if len(subset) == 0:
            continue
        jitter = np.random.default_rng(42).uniform(-0.15, 0.15, len(subset))
        c = colors_flt if "FLT" in grp else colors_gc
        ax.scatter([i] * len(subset) + jitter, subset.values, c=c, s=40, alpha=0.7,
                   edgecolors="white", linewidth=0.5, zorder=3)
        ax.plot([i - 0.25, i + 0.25], [subset.median(), subset.median()],
                c=c, linewidth=2, zorder=4)
        positions.append(i)

    ax.set_xticks(range(len(group_order)))
    ax.set_xticklabels([g.replace("\n", " / ") for g in group_order], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"{me_name} ({mod_color}) eigengene", fontsize=11)
    ax.set_title(f"Module eigengene: {mod_color} (n={len(mod_assign[mod_assign['module_color']==mod_color])} genes)", fontsize=13)
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")

    # Add vertical separators between Arm groups
    for sep in [1.5, 3.5, 5.5]:
        ax.axvline(sep, color="lightgrey", linewidth=0.5, linestyle=":")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


    ax.legend(handles=[Patch(facecolor=colors_flt, label="FLT"),
                       Patch(facecolor=colors_gc, label="GC")],
              loc="upper right", framealpha=0.8)

    fig.tight_layout()
    fig.savefig(plot_dir / f"eigengene_{mod_color}.pdf", dpi=150)
    fig.savefig(plot_dir / f"eigengene_{mod_color}.png", dpi=150)
    plt.close(fig)
    print(f"  Saved: eigengene_{mod_color}.pdf")

print("\n" + "="*60)
print("4. GO/KEGG Enrichment for Significant Modules")
print("="*60)



enrich_dir = OUTDIR / "enrichment"
enrich_dir.mkdir(exist_ok=True)

# Need gene symbols for enrichment
enrich_modules = ["grey60", "red", "royalblue", "salmon", "pink", "green",
                  "midnightblue", "blue"]

all_enrich = []
for mod_color in enrich_modules:
    mod_genes = mod_assign[mod_assign["module_color"] == mod_color]["gene"].tolist()
    symbols = [sym_map.get(g, "") for g in mod_genes]
    symbols = [s for s in symbols if s and s != "nan"]
    if len(symbols) < 5:
        print(f"  {mod_color}: too few mapped symbols ({len(symbols)}), skipping")
        continue

    # Background: all genes in WGCNA
    bg_symbols = [sym_map.get(g, "") for g in mod_assign["gene"]]
    bg_symbols = [s for s in bg_symbols if s and s != "nan"]

    for lib in ["GO_Biological_Process_2023", "KEGG_2021_Mouse", "Reactome_2022"]:
        try:
            enr = gp.enrichr(gene_list=symbols, gene_sets=lib,
                             background=bg_symbols,
                             organism="mouse", no_plot=True, cutoff=1.0)
            res = enr.results
            if len(res) > 0:
                res["module"] = mod_color
                res["library"] = lib
                all_enrich.append(res[res["Adjusted P-value"] < 0.10].head(10))
        except Exception as e:
            print(f"  {mod_color}/{lib}: {e}")

if all_enrich:
    enrich_all = pd.concat(all_enrich, ignore_index=True)
    enrich_all.to_csv(enrich_dir / "go_kegg_enrichment_all.csv", index=False)

    print(f"\n  Top enrichments per module (q < 0.10):")
    for mod in enrich_modules:
        mod_res = enrich_all[enrich_all["module"] == mod].head(5)
        if len(mod_res) > 0:
            print(f"\n  [{mod}]")
            for _, r in mod_res.iterrows():
                term = r["Term"][:60]
                print(f"    {term:<60} q={r['Adjusted P-value']:.4f} ({r['library']})")
else:
    print("  No enrichment results obtained")

print("\n" + "="*60)
print("6. Final Integrated Table")
print("="*60)

pres = pd.read_csv(WDIR / "module_preservation.csv")

final_rows = []
for mod_color in SIG_MODULES:
    me_name = color_to_me.get(mod_color)
    if me_name is None:
        continue

    size = len(mod_assign[mod_assign["module_color"] == mod_color])

    # Simple effects
    mod_contrasts = cdf[cdf["color"] == mod_color]
    se_strs = {}
    for _, r in mod_contrasts.iterrows():
        sig = "***" if r["q"] < 0.001 else "**" if r["q"] < 0.01 else "*" if r["q"] < 0.05 else "†" if r["q"] < 0.10 else ""
        se_strs[r["stratum"]] = f"{r['estimate']:+.3f}{sig}"

    # Preservation
    p_row = pres[pres["module_color"] == mod_color]
    zsummary = p_row["Zsummary"].values[0] if len(p_row) else np.nan
    pres_label = p_row["preservation"].values[0] if len(p_row) else "?"

    # Top pathway from curated
    pw_enrich = pd.read_csv(WDIR / "module_pathway_enrichment.csv")
    mod_pw = pw_enrich[(pw_enrich["module"] == mod_color) & (pw_enrich["q_value"] < 0.20)].sort_values("q_value")
    top_pw = mod_pw.iloc[0]["pathway"] if len(mod_pw) > 0 else "-"
    top_pw_q = f"{mod_pw.iloc[0]['q_value']:.4f}" if len(mod_pw) > 0 else "-"

    # Top hub genes
    if mod_color == "grey60":
        hubs = ", ".join(kme_df.head(5)["symbol"].tolist())
    else:
        # Quick kME for this module
        mg = mod_assign[mod_assign["module_color"] == mod_color]["gene"].tolist()
        me_v = MEs.loc[meta_fg.index, me_name].values
        kmes = []
        for g in mg:
            if g in rtech_fg.index:
                r, _ = sp_stats.pearsonr(rtech_fg.loc[g, meta_fg.index].values, me_v)
                kmes.append((g, abs(r)))
        kmes.sort(key=lambda x: -x[1])
        hubs = ", ".join([sym_map.get(g, g)[:10] for g, _ in kmes[:5]])

    # GO/KEGG top term
    if all_enrich:
        mod_go = enrich_all[enrich_all["module"] == mod_color].head(1)
        go_term = mod_go["Term"].values[0][:40] if len(mod_go) else "-"
    else:
        go_term = "-"

    final_rows.append({
        "module": me_name, "color": mod_color, "size": size,
        "curated_pathway": top_pw, "curated_q": top_pw_q,
        "GO_KEGG_top": go_term,
        "hub_genes": hubs,
        "ISS-T_YNG": se_strs.get("ISS-T_YNG", ""),
        "ISS-T_OLD": se_strs.get("ISS-T_OLD", ""),
        "LAR_YNG": se_strs.get("LAR_YNG", ""),
        "LAR_OLD": se_strs.get("LAR_OLD", ""),
        "Zsummary": round(zsummary, 1),
        "preservation": pres_label,
    })

final_df = pd.DataFrame(final_rows)
final_df.to_csv(OUTDIR / "final_integrated_table.csv", index=False)

print(f"\n{'Color':<14} {'n':>4} {'ISS-T YNG':>11} {'ISS-T OLD':>11} {'LAR YNG':>11} {'LAR OLD':>11} {'Zsum':>6} {'Curated':>18}")
print("-"*100)
for _, r in final_df.iterrows():
    print(f"{r['color']:<14} {r['size']:>4} {r['ISS-T_YNG']:>11} {r['ISS-T_OLD']:>11} {r['LAR_YNG']:>11} {r['LAR_OLD']:>11} {r['Zsummary']:>6} {r['curated_pathway']:>18}")

print(f"\nAll outputs saved to: {OUTDIR}")
