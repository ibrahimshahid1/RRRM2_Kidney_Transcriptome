#!/usr/bin/env python3
"""Component decomposition of the cross-cohort signature score.

Which part of the a-priori signature actually generalizes across cohorts:
the remodeling(up) axis, the DCT/NCC-WNK transport(down) axis, or the
DCT2/aldosterone(down) axis? Each scored in its flight-predicting direction,
with a within-cohort label-permutation null. No supervised fitting -> no overfit.
"""
import os, json, warnings
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
warnings.filterwarnings("ignore")
RNG = np.random.default_rng(20260701)
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OSDR = os.path.join(ROOT, "data", "external", "osdr")
OUT = os.path.join(ROOT, "data", "results", "run_20260701_sf_classifier")

VST = {"OSD-102":"OSD-102/GLDS-102_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
       "OSD-163":"OSD-163/GLDS-163_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
       "OSD-253":"OSD-253/GLDS-253_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
       "OSD-513":"OSD-513/GLDS-513_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
       "OSD-462":"OSD-462/GLDS-462_rna_seq_VST_Counts_mRNA_GLbulkRNAseq.csv"}
DE = os.path.join(OSDR,"OSD-462/GLDS-462_rna_seq_differential_expression_mRNA_GLbulkRNAseq.csv")
TRANSPORT=["Slc12a3","Stk39","Oxsr1","Wnk1","Wnk4","Wnk3","Kcnj10","Kcnj16","Clcnkb","Bsnd","Klhl3","Cul3","Calb1","Trpm6","Pvalb","Ppp1r1a"]
DCT2=["Trpv5","Calb1","Scnn1a","Scnn1b","Scnn1g","Sgk1","Nr3c2","Hsd11b2","Klf15","Fxyd2","Tsc22d3","Per1"]
REMODEL=["Col1a1","Col1a2","Col3a1","Col4a1","Col5a2","Fn1","Sparc","Lox","Loxl1","Timp1","Vim","Acta2","Tgfb1","Pecam1","Cdh5","Kdr","Emcn","Eng","Tek","Vwf"]

def lab(c):
    t=c.replace('"','').split("_")
    return 1 if "FLT" in t else (0 if "GC" in t else np.nan)
coh={}; common=None
for cid,f in VST.items():
    df=pd.read_csv(os.path.join(OSDR,f),index_col=0); df.index=[str(i).replace('"','') for i in df.index]
    keep=[c for c in df.columns if not np.isnan(lab(c))]
    coh[cid]=(df[keep], np.array([int(lab(c)) for c in keep]))
    common=set(df.index) if common is None else common&set(df.index)
common=sorted(common); cs=set(common)
de=pd.read_csv(DE,usecols=["ENSEMBL","SYMBOL"]); s2e={}
for e,s in zip(de.ENSEMBL.astype(str),de.SYMBOL.astype(str)): s2e.setdefault(s,e)
res=lambda ss:[s2e[x] for x in ss if s2e.get(x) in cs]
pT,pD,pR=res(TRANSPORT),res(DCT2),res(REMODEL)

def z(df):
    m=df.mean(1); s=df.std(1).replace(0,np.nan); return df.sub(m,0).div(s,0).fillna(0.0)
Z={cid:z(df.loc[common]) for cid,(df,_) in coh.items()}
Y={cid:yy for cid,(_,yy) in coh.items()}

# component scores in flight-predicting direction (higher = flight)
def comp_scores(cid):
    zt=Z[cid]
    return {"remodel_up":  zt.loc[pR].mean(0).values,
            "transport_down": -zt.loc[pT].mean(0).values,
            "dct2_aldo_down": -zt.loc[pD].mean(0).values,
            "combined": zt.loc[pR].mean(0).values - zt.loc[pT+pD].mean(0).values}
COMPS=["remodel_up","transport_down","dct2_aldo_down","combined"]

def auc(y,s): return roc_auc_score(y,s) if len(set(y))==2 else np.nan
rows=[]
for comp in COMPS:
    per={}
    for cid in VST:
        sc=comp_scores(cid)[comp]; per[cid]=auc(Y[cid],sc)
    mean=np.nanmean(list(per.values()))
    # within-cohort permutation null on the mean-per-cohort AUC
    nperm=2000; null=np.empty(nperm)
    for i in range(nperm):
        a=[]
        for cid in VST:
            yp=RNG.permutation(Y[cid]); a.append(auc(yp,comp_scores(cid)[comp]))
        null[i]=np.nanmean(a)
    p=(1+np.sum(null>=mean))/(nperm+1)
    rows.append({"component":comp,"mean_auc":round(float(mean),3),"perm_p":round(float(p),4),
                 **{f"auc_{cid}":round(float(per[cid]),3) for cid in VST}})
    print(f"{comp:16s} mean AUC={mean:.3f}  perm p={p:.4f}  | " +
          "  ".join(f"{cid.split('-')[1]}:{per[cid]:.2f}" for cid in VST))
pd.DataFrame(rows).to_csv(os.path.join(OUT,"signature_component_auc.tsv"),sep="\t",index=False)
json.dump(rows, open(os.path.join(OUT,"signature_component_auc.json"),"w"), indent=2)
print("saved ->", os.path.join(OUT,"signature_component_auc.tsv"))
