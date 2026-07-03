#!/usr/bin/env python3
"""Feature-specificity negative control (random-gene score null) + summary figure.

OSD-105 muscle (a non-kidney TISSUE control) is not on disk, so the download-free
negative control here is a FEATURE-specificity control: do random size-matched
up/down gene sets reach the a-priori signature's cross-cohort AUC? If not, the
biological signature is specifically informative, not a generic property of any
gene set. (A non-kidney tissue control remains a recommended one-time add.)
"""
import os, json, warnings
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
RNG=np.random.default_rng(7)
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OSDR=os.path.join(ROOT,"data","external","osdr")
OUT=os.path.join(ROOT,"data","results","run_20260701_sf_classifier")
VST={"OSD-102":"OSD-102/GLDS-102_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
     "OSD-163":"OSD-163/GLDS-163_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
     "OSD-253":"OSD-253/GLDS-253_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
     "OSD-513":"OSD-513/GLDS-513_rna_seq_VST_Counts_rRNArm_GLbulkRNAseq.csv",
     "OSD-462":"OSD-462/GLDS-462_rna_seq_VST_Counts_mRNA_GLbulkRNAseq.csv"}
DE=os.path.join(OSDR,"OSD-462/GLDS-462_rna_seq_differential_expression_mRNA_GLbulkRNAseq.csv")
TRANSPORT=["Slc12a3","Stk39","Oxsr1","Wnk1","Wnk4","Wnk3","Kcnj10","Kcnj16","Clcnkb","Bsnd","Klhl3","Cul3","Calb1","Trpm6","Pvalb","Ppp1r1a"]
DCT2=["Trpv5","Calb1","Scnn1a","Scnn1b","Scnn1g","Sgk1","Nr3c2","Hsd11b2","Klf15","Fxyd2","Tsc22d3","Per1"]
REMODEL=["Col1a1","Col1a2","Col3a1","Col4a1","Col5a2","Fn1","Sparc","Lox","Loxl1","Timp1","Vim","Acta2","Tgfb1","Pecam1","Cdh5","Kdr","Emcn","Eng","Tek","Vwf"]
lab=lambda c:1 if "FLT" in c.replace('"','').split("_") else (0 if "GC" in c.replace('"','').split("_") else np.nan)
coh={}; common=None
for cid,f in VST.items():
    df=pd.read_csv(os.path.join(OSDR,f),index_col=0); df.index=[str(i).replace('"','') for i in df.index]
    keep=[c for c in df.columns if not np.isnan(lab(c))]
    coh[cid]=(df[keep],np.array([int(lab(c)) for c in keep])); common=set(df.index) if common is None else common&set(df.index)
common=sorted(common); cs=set(common)
de=pd.read_csv(DE,usecols=["ENSEMBL","SYMBOL"]); s2e={}
for e,s in zip(de.ENSEMBL.astype(str),de.SYMBOL.astype(str)): s2e.setdefault(s,e)
res=lambda ss:[s2e[x] for x in ss if s2e.get(x) in cs]
pT,pD,pR=res(TRANSPORT),res(DCT2),res(REMODEL); nUP,nDN=len(pR),len(pT)+len(pD)
def z(df):
    m=df.mean(1); s=df.std(1).replace(0,np.nan); return df.sub(m,0).div(s,0).fillna(0.0)
Z={cid:z(df.loc[common]) for cid,(df,_) in coh.items()}; Y={cid:yy for cid,(_,yy) in coh.items()}
auc=lambda y,s: roc_auc_score(y,s) if len(set(y))==2 else np.nan
def mean_auc(up,dn):
    a=[]
    for cid in VST:
        sc=Z[cid].loc[up].mean(0).values - Z[cid].loc[dn].mean(0).values
        a.append(auc(Y[cid],sc))
    return np.nanmean(a)
obs=mean_auc(pR,pT+pD)
# random-gene score null: random up/down sets size-matched
ndraw=2000; null=np.empty(ndraw)
for i in range(ndraw):
    pick=RNG.choice(common,size=nUP+nDN,replace=False)
    null[i]=mean_auc(list(pick[:nUP]),list(pick[nUP:]))
p_rand=(1+np.sum(null>=obs))/(ndraw+1)
print(f"[neg-ctrl] a-priori combined score AUC={obs:.3f}  random-gene-score null mean={null.mean():.3f} "
      f"95pct={np.quantile(null,0.95):.3f}  p={p_rand:.4f}")
json.dump({"observed_combined_auc":round(float(obs),4),"random_score_null_mean":round(float(null.mean()),4),
           "random_score_null_95pct":round(float(np.quantile(null,0.95)),4),"p_signature_vs_random_score":round(float(p_rand),4),
           "note":"OSD-105 non-kidney TISSUE control not on disk; recommended one-time add."},
          open(os.path.join(OUT,"negative_control.json"),"w"),indent=2)

# ---- figure ----
comp=pd.read_csv(os.path.join(OUT,"signature_component_auc.tsv"),sep="\t")
order=["remodel_up","transport_down","dct2_aldo_down","combined"]
labels=["Remodeling/ECM\n(up)","DCT/NCC-WNK\ntransport (down)","DCT2/aldosterone\n(down)","Combined\nsignature"]
comp=comp.set_index("component").loc[order]
cids=["OSD-102","OSD-163","OSD-253","OSD-513","OSD-462"]
fig,ax=plt.subplots(figsize=(8.4,5.0))
xs=np.arange(len(order)); bars=comp["mean_auc"].values
cols=["#1f7a4d","#b5651d","#8a8a8a","#204a78"]
ax.bar(xs,bars,color=cols,alpha=.85,width=.62,zorder=2)
for i,cid in enumerate(cids):
    ax.scatter(xs+RNG.uniform(-.16,.16,len(order)),[comp.loc[o,f"auc_{cid}"] for o in order],
               s=34,color="k",alpha=.55,zorder=3,label="per-cohort AUC" if i==0 else None)
ax.axhline(.5,ls="--",c="#444",lw=1,zorder=1); ax.axhline(np.quantile(null,0.95),ls=":",c="#a11",lw=1.2,zorder=1)
ax.text(len(order)-.5,np.quantile(null,0.95)+.006,"random-gene-score 95th pct",ha="right",va="bottom",color="#a11",fontsize=8)
for i,o in enumerate(order):
    ax.text(i,bars[i]+.02,f"AUC {bars[i]:.2f}\np={comp.loc[o,'perm_p']:.3f}",ha="center",va="bottom",fontsize=8.5,fontweight="bold")
ax.set_xticks(xs); ax.set_xticklabels(labels,fontsize=9); ax.set_ylim(0.15,1.02)
ax.set_ylabel("Flight-vs-ground AUC (mean across 5 cohorts)",fontsize=10)
ax.set_title("Cross-cohort portability of the spaceflight kidney signature\n(within-cohort z; leave-one-cohort-out in spirit; label-permutation p)",fontsize=10.5)
ax.legend(loc="lower left",fontsize=8,framealpha=.9)
plt.tight_layout()
for ext in ("png","pdf"):
    fig.savefig(os.path.join(OUT,f"signature_portability.{ext}"),dpi=150,bbox_inches="tight")
print("saved figure ->",os.path.join(OUT,"signature_portability.png"))
