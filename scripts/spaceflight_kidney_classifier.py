#!/usr/bin/env python3
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
import os, json, warnings
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
warnings.filterwarnings("ignore")

RNG = np.random.default_rng(20260701)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OSDR = os.path.join(ROOT, "data", "external", "osdr")
OUT = os.path.join(ROOT, "data", "results", "run_20260701_sf_classifier")
os.makedirs(OUT, exist_ok=True)

VST = {
    "OSD-102": "OSD-102/GLDS-102_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
    "OSD-163": "OSD-163/GLDS-163_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
    "OSD-253": "OSD-253/GLDS-253_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
    "OSD-513": "OSD-513/GLDS-513_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
    "OSD-462": "OSD-462/GLDS-462_rna_seq_VST_Counts_mRNA_GLbulkRNAseq.csv",
}
DE_MAP = os.path.join(OSDR, "OSD-462/GLDS-462_rna_seq_differential_expression_mRNA_GLbulkRNAseq.csv")

# ---- a-priori signature (symbols), from the v11 manuscript gene sets ----
TRANSPORT_DOWN = ["Slc12a3","Stk39","Oxsr1","Wnk1","Wnk4","Wnk3","Kcnj10","Kcnj16",
                  "Clcnkb","Bsnd","Klhl3","Cul3","Calb1","Trpm6","Pvalb","Ppp1r1a"]
DCT2_ALDO_DOWN = ["Trpv5","Calb1","Scnn1a","Scnn1b","Scnn1g","Sgk1","Nr3c2","Hsd11b2",
                  "Klf15","Fxyd2","Tsc22d3","Per1"]
REMODEL_UP = ["Col1a1","Col1a2","Col3a1","Col4a1","Col5a2","Fn1","Sparc","Lox","Loxl1",
              "Timp1","Vim","Acta2","Tgfb1","Pecam1","Cdh5","Kdr","Emcn","Eng","Tek","Vwf"]

def parse_label(col):
    toks = col.replace('"','').split("_")
    if "FLT" in toks: return 1        # flight
    if "GC" in toks:  return 0        # ground control
    return np.nan                     # BSL / baseline / other -> excluded

def load_cohort(cid):
    df = pd.read_csv(os.path.join(OSDR, VST[cid]), index_col=0)
    df.index = [str(i).replace('"','') for i in df.index]
    labels = {c: parse_label(c) for c in df.columns}
    keep = [c for c in df.columns if not np.isnan(labels[c])]
    df = df[keep]
    y = np.array([labels[c] for c in keep], dtype=int)
    return df, y

# ---- load, intersect genes ----
cohorts, ys, common = {}, {}, None
for cid in VST:
    df, y = load_cohort(cid)
    cohorts[cid], ys[cid] = df, y
    common = set(df.index) if common is None else (common & set(df.index))
common = sorted(common)
print(f"[genes] common Ensembl IDs across 5 cohorts: {len(common)}")
for cid in VST:
    print(f"  {cid}: {cohorts[cid].shape[1]} samples  (flight={int(ys[cid].sum())}, ground={int((ys[cid]==0).sum())})")

# ---- Ensembl<->symbol map, resolve signature panels to Ensembl in common set ----
de = pd.read_csv(DE_MAP, usecols=["ENSEMBL","SYMBOL"])
sym2ens = {}
for e, s in zip(de["ENSEMBL"].astype(str), de["SYMBOL"].astype(str)):
    sym2ens.setdefault(s, e)
commonset = set(common)
def resolve(symbols):
    out = []
    for s in symbols:
        e = sym2ens.get(s)
        if e in commonset: out.append(e)
    return out
panel_transport = resolve(TRANSPORT_DOWN)
panel_dct2      = resolve(DCT2_ALDO_DOWN)
panel_remodel   = resolve(REMODEL_UP)
panel_signature = sorted(set(panel_transport) | set(panel_dct2) | set(panel_remodel))
print(f"[panels] transport_down={len(panel_transport)} dct2_aldo={len(panel_dct2)} "
      f"remodel_up={len(panel_remodel)} signature_union={len(panel_signature)}")

# ---- build within-cohort z-scored stacked matrix ----
def zscore(df):
    m = df.mean(axis=1); s = df.std(axis=1).replace(0, np.nan)
    z = df.sub(m, axis=0).div(s, axis=0).fillna(0.0)
    return z
X_parts, y_all, grp = [], [], []
for cid in VST:
    z = zscore(cohorts[cid].loc[common])   # genes x samples, z within cohort
    X_parts.append(z.T)                    # samples x genes
    y_all.append(ys[cid]); grp += [cid]*z.shape[1]
X = pd.concat(X_parts).loc[:, common]
y = np.concatenate(y_all); groups = np.array(grp)
print(f"[matrix] X = {X.shape[0]} samples x {X.shape[1]} genes")

DIRS = {g: (1.0 if g in set(panel_remodel) else -1.0) for g in panel_signature}  # up vs down

def signature_score(Xsub):
    # mean z of up-genes minus mean z of down-genes  -> higher = more flight-like
    up = [g for g in panel_remodel if g in Xsub.columns]
    dn = [g for g in (panel_transport+panel_dct2) if g in Xsub.columns]
    return Xsub[up].mean(axis=1) - Xsub[dn].mean(axis=1)

def auc_safe(yt, sc):
    return roc_auc_score(yt, sc) if len(set(yt)) == 2 else np.nan

# ---- (A) a-priori signature score, per held-out cohort ----
score_all = signature_score(X)
score_rows = []
for cid in VST:
    m = groups == cid
    score_rows.append({"cohort": cid, "auc": auc_safe(y[m], score_all[m].values),
                       "n_flt": int(y[m].sum()), "n_gc": int((y[m]==0).sum())})
score_df = pd.DataFrame(score_rows)
sig_score_mean = score_df["auc"].mean()
print(f"[A] a-priori signature-score LOCO-style per-cohort mean AUC = {sig_score_mean:.3f}")

# ---- (B) supervised LOCO logistic regression across feature panels ----
PANELS = {
    "full_transcriptome": common,
    "signature_union": panel_signature,
    "transport_down": panel_transport,
    "dct2_aldo_down": panel_dct2,
    "remodel_up": panel_remodel,
}
def loco_auc(cols, Xm=None, ym=None):
    Xm = X if Xm is None else Xm; ym = y if ym is None else ym
    rows = []
    for cid in VST:
        te = groups == cid; tr = ~te
        if len(set(ym[tr])) < 2 or len(set(ym[te])) < 2:
            rows.append((cid, np.nan)); continue
        clf = LogisticRegression(penalty="l2", C=0.1, class_weight="balanced", max_iter=2000)
        clf.fit(Xm.loc[tr, cols].values, ym[tr])
        p = clf.predict_proba(Xm.loc[te, cols].values)[:, 1]
        rows.append((cid, auc_safe(ym[te], p)))
    return rows

def cv_within_auc(cols):
    aucs = []
    for cid in VST:
        m = groups == cid
        Xi, yi = X.loc[m, cols].values, y[m]
        if yi.sum() < 2 or (yi==0).sum() < 2: continue
        skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=0)
        pr = np.zeros(len(yi))
        for tr, te in skf.split(Xi, yi):
            clf = LogisticRegression(penalty="l2", C=0.1, class_weight="balanced", max_iter=2000)
            clf.fit(Xi[tr], yi[tr]); pr[te] = clf.predict_proba(Xi[te])[:, 1]
        aucs.append(auc_safe(yi, pr))
    return np.nanmean(aucs)

loco_rows, summary = [], []
for name, cols in PANELS.items():
    rr = loco_auc(cols)
    for cid, a in rr:
        loco_rows.append({"panel": name, "held_out": cid, "auc": a})
    mean_loco = np.nanmean([a for _, a in rr])
    within = cv_within_auc(cols)
    summary.append({"panel": name, "n_features": len(cols),
                    "mean_loco_auc": round(mean_loco, 4),
                    "within_cohort_cv_auc": round(float(within), 4)})
    print(f"[B] {name:20s} n={len(cols):5d}  LOCO AUC={mean_loco:.3f}  within-CV AUC={within:.3f}")
loco_df = pd.DataFrame(loco_rows); summ_df = pd.DataFrame(summary)

# ---- (C) permutation null: shuffle labels within cohort ----
def perm_null(metric_fn, nperm=500):
    obs = metric_fn(y)
    null = np.empty(nperm)
    for i in range(nperm):
        yp = y.copy()
        for cid in VST:
            m = groups == cid; idx = np.where(m)[0]
            yp[idx] = RNG.permutation(y[idx])
        null[i] = metric_fn(yp)
    p = (1 + np.sum(null >= obs)) / (nperm + 1)
    return obs, null, p

def sig_score_metric(yy):
    return np.nanmean([auc_safe(yy[groups==cid], score_all[groups==cid].values) for cid in VST])
def sig_panel_metric(yy):
    return np.nanmean([a for _, a in loco_auc(panel_signature, ym=yy)])

obs_s, null_s, p_s = perm_null(sig_score_metric, nperm=1000)
obs_p, null_p, p_p = perm_null(sig_panel_metric, nperm=200)
print(f"[C] signature-score AUC={obs_s:.3f}  perm p={p_s:.4f}  (null mean {null_s.mean():.3f})")
print(f"[C] signature-panel LOCO AUC={obs_p:.3f}  perm p={p_p:.4f}  (null mean {null_p.mean():.3f})")

# ---- (D) random size-matched gene panels ----
k = len(panel_signature); ndraw = 200
rand_aucs = np.empty(ndraw)
for i in range(ndraw):
    cols = list(RNG.choice(common, size=k, replace=False))
    rand_aucs[i] = np.nanmean([a for _, a in loco_auc(cols)])
sig_loco = summ_df.loc[summ_df.panel=="signature_union","mean_loco_auc"].values[0]
rand_p = (1 + np.sum(rand_aucs >= sig_loco)) / (ndraw + 1)
print(f"[D] random {k}-gene panels LOCO AUC: mean={rand_aucs.mean():.3f} "
      f"95pct={np.quantile(rand_aucs,0.95):.3f}  signature={sig_loco:.3f}  p={rand_p:.4f}")

# ---- save ----
score_df.to_csv(os.path.join(OUT, "signature_score_auc_by_cohort.tsv"), sep="\t", index=False)
pd.DataFrame({"sample": X.index, "cohort": groups, "y_flight": y,
             "signature_score": score_all.values}).to_csv(
             os.path.join(OUT, "signature_score_per_sample.tsv"), sep="\t", index=False)
loco_df.to_csv(os.path.join(OUT, "loco_auc_by_panel.tsv"), sep="\t", index=False)
summ_df.to_csv(os.path.join(OUT, "loco_auc_summary.tsv"), sep="\t", index=False)
pd.DataFrame({"random_panel_loco_auc": rand_aucs}).to_csv(
    os.path.join(OUT, "random_panel_null.tsv"), sep="\t", index=False)
json.dump({
    "n_common_genes": len(common),
    "cohorts": {cid: {"flight": int(ys[cid].sum()), "ground": int((ys[cid]==0).sum())} for cid in VST},
    "panels": {"transport_down": len(panel_transport), "dct2_aldo": len(panel_dct2),
               "remodel_up": len(panel_remodel), "signature_union": len(panel_signature)},
    "signature_score_mean_auc": round(float(obs_s), 4), "signature_score_perm_p": round(float(p_s), 4),
    "signature_panel_loco_auc": round(float(obs_p), 4), "signature_panel_perm_p": round(float(p_p), 4),
    "random_panel_mean_auc": round(float(rand_aucs.mean()), 4),
    "random_panel_95pct": round(float(np.quantile(rand_aucs, 0.95)), 4),
    "signature_vs_random_p": round(float(rand_p), 4),
    "loco_summary": summary,
}, open(os.path.join(OUT, "summary.json"), "w"), indent=2)
print(f"\n[done] outputs -> {OUT}")
