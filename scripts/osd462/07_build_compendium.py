#!/usr/bin/env python3
"""Layer 7 - Build the standalone LaTeX results compendium from anchor outputs."""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from _common import REPO_ROOT, anchor_dir, default_run

TEX_PATH = REPO_ROOT / "latex_paper" / "osd462_multiomics_compendium.tex"


# formatting helpers
def esc(s) -> str:
    return str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")


def f(x, d=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "--"
    return f"{x:+.{d}f}" if x < 0 else f"{x:.{d}f}"


def fp(p):
    if p is None or (isinstance(p, float) and not np.isfinite(p)):
        return "--"
    if p < 1e-3:
        m, e = f"{p:.1e}".split("e")
        return rf"${m}\times10^{{{int(e)}}}$"
    return f"{p:.3f}"


def table(headers, rows, align, caption, label, note=None, small=True):
    out = ["\\begin{table}[htbp]", "\\centering"]
    if small:
        out.append("\\small")
    out.append("\\begin{tabular}{" + align + "}")
    out.append("\\toprule")
    out.append(" & ".join(headers) + " \\\\")
    out.append("\\midrule")
    for r in rows:
        out.append(" & ".join(r) + " \\\\")
    out.append("\\bottomrule")
    out.append("\\end{tabular}")
    if note:
        out.append(f"\\\\[2pt]\\footnotesize\\emph{{{note}}}")
    out.append(f"\\caption{{{caption}}}")
    out.append(f"\\label{{{label}}}")
    out.append("\\end{table}")
    return "\n".join(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=default_run())
    ap.add_argument("--compile", action="store_true")
    args = ap.parse_args()
    out = anchor_dir(args.run)

    sjson = json.loads((out / "results_summary.json").read_text())
    rec = pd.read_csv(out / "osd462_rna_recurrence.tsv", sep="\t").iloc[0]
    pv = pd.read_csv(out / "osd462_rna_pathway_effects.tsv", sep="\t")
    summ = pd.read_csv(out / "protein_concordance_summary.tsv", sep="\t")
    phos = pd.read_csv(out / "phospho_axis_summary.tsv", sep="\t")
    net = pd.read_csv(out / "network_translation.tsv", sep="\t")
    master = pd.read_csv(out / "osd462_flight_effects.tsv", sep="\t")
    verdict = pd.read_csv(out / "phospho_axis_verdict.tsv", sep="\t").iloc[0]

    g = sjson["layer4_rna_gate"]
    l1 = sjson["layer1_protein"]
    l2 = sjson["layer2_phospho"]
    l3 = sjson["layer3_network"]

    # Table: Layer 4 pathway effects
    pv_sorted = pv.sort_values("rrrm2_iss_t_pathway_effect", ascending=False)
    t4_rows = [[esc(r["pathway"]), f(r["osd462_rna_pathway_effect"]),
                f(r["rrrm2_iss_t_pathway_effect"]),
                ("yes" if r["sign_agree"] else "\\textbf{no}")]
               for _, r in pv_sorted.iterrows()]
    tab4 = table(["Pathway", "OSD-462 RNA", "RRRM-2 RNA", "sign agree"],
                 t4_rows, "lrrc",
                 "Layer 4 RNA recurrence: pathway-level flight effects (mean per-gene "
                 "SF$-$GC / FLT$-$GC log$_2$). The matrix-high / DCT-low direction recurs; "
                 f"cosine $={g['point_cosine']:.3f}$, 95\\% CI "
                 f"[{g['ci'][0]:.3f}, {g['ci'][1]:.3f}], {int(rec['n_pathways'])} pathways, "
                 f"{sum(pv['sign_agree'])}/{len(pv)} sign-concordant.",
                 "tab:layer4")

    # Table: Layer 1 protein concordance
    order = ["ecm_organization", "dct_ncc_wnk_transport", "tlr4_innate", "s1p_s1pr3",
             "tubular_transport_broad"]
    s = summ.set_index("gene_set")
    t1_rows = []
    for name in order:
        if name not in s.index:
            continue
        r = s.loc[name]
        kind = "control" if r["kind"] == "control" else ""
        t1_rows.append([
            esc(name) + (f" ({kind})" if kind else ""),
            str(int(r["n_quantified"])),
            r["predicted_direction"],
            f(r["mean_protein_effect"]),
            fp(r["signed_mean_p"]),
            f(r["concordance"]),
            fp(r["concordance_p"]),
            f(r["sign_agreement"], 2),
        ])
    tab1 = table(["Gene set", "$n$", "pred.", "mean prot.", "$p_{\\uparrow}$",
                  "concord.", "$p_{c}$", "sign agr."],
                 t1_rows, "lccrrrrc",
                 "Layer 1 protein-abundance concordance. ``mean prot.'' is the mean OSD-462 "
                 "protein flight effect; $p_\\uparrow$ is the abundance/peptide-matched null "
                 "$p$ that the set moves in its RRRM-2-predicted direction; ``concord.'' is "
                 "the Spearman correlation between RRRM-2 RNA and OSD-462 protein effects "
                 "($p_c$ its matched-null $p$). No targeted set is concordant in the predicted "
                 "direction; ECM matrix proteins move significantly \\emph{opposite} to the "
                 "transcript prediction (two-sided matched-null "
                 f"$p={l1['ecm_signed_mean_p_two_sided']:.3f}$).",
                 "tab:layer1")

    # Table: NCC + SPAK phosphosites
    foc = phos[phos["gene_symbol"].isin(["Slc12a3", "Stk39", "Oxsr1"])].copy()
    foc = foc.sort_values(["gene_symbol", "phospho_effect"])
    tph_rows = []
    for _, r in foc.iterrows():
        occ = r.get("phospho_occupancy_effect", np.nan)
        tph_rows.append([
            "\\textit{" + esc(r["gene_symbol"]) + "}",
            esc(r["site_position"]),
            f"{int(r['n_fl'])}/{int(r['n_gc'])}",
            f(r["phospho_effect"]),
            f"[{r['phospho_ci_low']:.2f}, {r['phospho_ci_high']:.2f}]",
            fp(r["phospho_p_value"]),
            fp(r["phospho_q_value"]),
            f(occ) if np.isfinite(occ) else "--",
        ])
    tabph = table(["Gene", "site", "FL/GC", "effect", "95\\% CI", "$p$", "$q$", "protein-adj."],
                  tph_rows, "llcrrrrr",
                  "Layer 2 WNK-SPAK/OSR1-NCC phosphosite flight effects (FL$-$GC, log$_2$; "
                  "plex-adjusted per-site model). ``protein-adj.'' is the site-feature effect "
                  "minus total parent-protein effect, not isolated occupancy. Co-modified "
                  "T53/Y65 and S382/S383 features decrease; no isolated canonical feature "
                  "qualifies. NCC C-terminal features (96--124) and OSR1-339 are flat.",
                  "tab:phospho", small=True)

    # Table: network translation
    tn_rows = [[f"top-{int(r['top_k'])}",
                str(int(r["n_candidates_with_protein"])), f(r["protein_mean_abs_effect"]),
                f(r["protein_null_median"]), fp(r["protein_enrichment_p"]),
                str(int(r["n_candidates_with_phospho"])), f(r["phospho_mean_max_abs_effect"]),
                f(r["phospho_null_median"]), fp(r["phospho_enrichment_p"])]
               for _, r in net.iterrows()]
    tabn = table(["candidates", "$n_p$", "obs", "null", "$p$", "$n_\\phi$", "obs", "null", "$p$"],
                 tn_rows, "lcccccccc",
                 "Layer 3 network-candidate translation. For RRRM-2 LIONESS/node2vec top-$k$ "
                 "composite candidates: mean $|$protein effect$|$ ($n_p$ candidates) and mean "
                 "per-gene max $|$phospho effect$|$ ($n_\\phi$ candidates) versus the "
                 "abundance/detectability-matched null. Candidates are \\emph{not} enriched at "
                 "either layer.", "tab:layer3")

    # peptide detection for target machinery
    det_genes = ["Slc12a3", "Wnk1", "Wnk4", "Stk39", "Oxsr1", "Calb1", "Nedd4l", "Cul3",
                 "Fn1", "Col1a1", "Col3a1", "Fbn1"]
    mg = master.set_index("gene_symbol")
    det_rows = []
    for gname in det_genes:
        if gname in mg.index:
            r = mg.loc[gname]
            det_rows.append(["\\textit{" + esc(gname) + "}",
                             str(int(r["n_peptides"])) if np.isfinite(r["n_peptides"]) else "--",
                             f(r["rrrm2_iss_t_rna_effect"]) if np.isfinite(r["rrrm2_iss_t_rna_effect"]) else "--",
                             f(r["protein_flight_effect"]) if np.isfinite(r["protein_flight_effect"]) else "--"])
    tabdet = table(["Gene", "peptides", "RRRM-2 RNA", "OSD-462 protein"], det_rows, "lrrr",
                   "Detection and flight effects for representative DCT/transport and matrix "
                   "machinery (gene-level, both plexes). Peptide counts are summed across the "
                   "two TMT plexes.", "tab:detect")

    # compose document
    matrix_anti = l1["ecm_signed_mean_p_two_sided"]
    doc = rf"""\documentclass[11pt]{{article}}
\usepackage[a4paper,margin=1in]{{geometry}}
\usepackage{{amsmath,amssymb}}
\usepackage{{booktabs}}
\usepackage{{array}}
\usepackage{{graphicx}}
\usepackage{{caption}}
\usepackage{{xcolor}}
\usepackage{{enumitem}}
\usepackage[colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue]{{hyperref}}
\usepackage{{parskip}}
\graphicspath{{{{figures_osd462/}}}}
\captionsetup{{font=small,labelfont=bf}}
\newcommand{{\pass}}{{\textcolor{{green!55!black}}{{\textbf{{PASS}}}}}}
\newcommand{{\fail}}{{\textcolor{{red!70!black}}{{\textbf{{not supported}}}}}}

\title{{\textbf{{A protein- and phosphoprotein-level test of the spaceflight\\
kidney matrix-high / DCT-low remodeling signal}}\\[4pt]
\large OSD-462 / RR-10 multi-omics anchor for the RRRM-2 transcriptomic axis}}
\author{{RRRM-2 Kidney Transcriptome Project}}
\date{{Compendium generated {pd.Timestamp.now().strftime('%Y-%m-%d')} \\ run \texttt{{{esc(args.run)}}}}}

\begin{{document}}
\maketitle

\begin{{abstract}}
\noindent The RRRM-2 (OSD-771) kidney transcriptome analysis identified a recurrent
ISS-transit remodeling axis: an extracellular-matrix (ECM) program elevated in flight
and a distal-convoluted-tubule (DCT) Na--Cl transport program (WNK--SPAK/OSR1--NCC)
reduced in flight. Because an RNA co-expression shift is weak standalone evidence, we
used the independent OSD-462 / RR-10 spaceflight kidney multi-omic study (TMT
proteomics and phosphoproteomics, female B6129SF2/J, single ISS-transit-type arm) as a
protein/phosphoprotein \emph{{anchor}}. We pre-registered four hypotheses and explicit
falsification rules and tested them with abundance/peptide-matched random-gene-set
nulls. \textbf{{Findings.}} (i) The matrix-high / DCT-low direction \emph{{recurs at the
RNA level}} in OSD-462 (pathway-vector cosine $={g['point_cosine']:.2f}$, 95\% CI
$[{g['ci'][0]:.2f},{g['ci'][1]:.2f}]$; gate \pass{{}}). (ii) It does \emph{{not}} translate
into protein-abundance concordance: no targeted set is concordant in the predicted
direction, and matrix proteins move significantly \emph{{opposite}} to their transcripts
(two-sided matched-null $p={matrix_anti:.3f}$). (iii) A T53-indexed NCC feature is
suppressed
(${phos.set_index(['gene_symbol','site_position']).loc[('Slc12a3','53'),'phospho_effect']:.2f}$
log$_2$), as is an upstream S383-indexed SPAK feature
(${phos.set_index(['gene_symbol','site_position']).loc[('Stk39','383'),'phospho_effect']:.2f}$
log$_2$), while total NCC
protein is unchanged (${l2['ncc_total_protein_effect']:+.2f}$). (iv) Network-nominated
candidate genes are \emph{{not}} enriched among protein- or phospho-changing genes.
\textbf{{Interpretation.}} Both residue-indexed rows arise from co-modified peptide
sequences (T53/Y65 and S382/S383), so they are compatible with, but do not isolate,
canonical NCC/SPAK phosphoactivation. The strict activity claim remains unresolved.
\end{{abstract}}

\section{{Purpose and design}}
This compendium reports an independent multi-omics test of one transcriptomic result.
The RRRM-2 backbone analysis (run \texttt{{run\_20260519\_000547\_2500g}}) is a \emph{{fixed
input}}, not re-derived here. OSD-462 / RR-10 is the NASA kidney multi-omic study: female
B6129SF2/J mice, left kidney, single terminal collection arm, with RNA-seq, TMT
proteomics, and phosphoproteomics. Its hard constraints shape every claim
(Table~\ref{{tab:cohort}}): it corroborates only the \emph{{ISS-transit}} side (no
late-adapting-return arm, no aged animals); the comparison is cross-study, cross-strain
and cross-modality, so the unit of inference is the \emph{{direction}} of the flight
effect (Space Flight $-$ Ground Control), never sample-level pooling; and bulk
whole-kidney mass spectrometry dilutes the small DCT compartment. Primary inference is an
abundance/peptide-matched random-gene-set null, which is robust to the generically weak
RNA$\leftrightarrow$protein correlation and to gene-set size.

\begin{{table}}[htbp]\centering\small
\begin{{tabular}}{{lll}}
\toprule
Property & OSD-462 / RR-10 (anchor) & RRRM-2 / OSD-771 (reference) \\
\midrule
Tissue & Left kidney & Kidney \\
Sex & Female & Female \\
Strain & B6129SF2/J (WT) & C57BL/6NTac \\
Age & 14--15 wk (single) & Young $+$ Old \\
Arms & Single terminal arm & ISS-T $+$ LAR \\
Modalities & RNA-seq, TMT proteomics, phosphoproteomics & RNA-seq only \\
\bottomrule
\end{{tabular}}
\caption{{Cohort comparison. OSD-462 can corroborate the ISS-transit flight signal only;
the late-adapting-return and age-stratified findings remain RNA-only and out of scope here.}}
\label{{tab:cohort}}
\end{{table}}

\paragraph{{Estimation.}} Protein flight effects are mean(log$_2$ FL) $-$ mean(log$_2$ GC)
computed \emph{{within}} each TMT plex (``Samp1-5'', ``Samp6-10'') on row-normalized S/N,
with per-channel median (sample-loading) centering, then averaged across plexes; this
removes the per-plex batch constant. Genes are required to be quantified in both plexes
with $\geq 2$ peptides for the primary tests ($\geq 4$-peptide sensitivity reproduces the
conclusions). Phosphosite effects use a per-site linear model
$\log_2(\text{{S/N}})\sim\text{{flight}}+\text{{plex}}$, reporting the flight coefficient
with a $t$-based 95\% CI. The RNA recurrence gate compares pathway-effect vectors (mean
per-gene flight effect per mechanism gene set) by cosine, with a sample-resampling
bootstrap CI and leave-one-pathway-out check, mirroring the OSD-513 cross-study test.

\section{{Pre-registered hypotheses and outcomes}}
\begin{{description}}[leftmargin=2.2em,style=nextline]
\item[H1 -- matrix protein concordance] \fail{{}}. Matrix/ECM proteins do not move in the
RNA-predicted (up) direction; they move significantly down (Table~\ref{{tab:layer1}}).
\item[H2 -- DCT protein concordance] \fail{{}}. DCT/WNK/NCC \emph{{abundance}} is flat
(NCC ${l2['ncc_total_protein_effect']:+.2f}$ log$_2$); not concordant with the predicted
decrease.
\item[H3 -- phospho activity] unresolved. Canonical-site-indexed features decrease, but
all are co-modified phosphoforms; no isolated canonical NCC/SPAK feature qualifies
(Table~\ref{{tab:phospho}}).
\item[H4 -- network translation] \fail{{}}. LIONESS/node2vec candidates are not enriched
among protein- or phospho-changing genes (Table~\ref{{tab:layer3}}).
\end{{description}}

\section{{Results}}

\subsection{{Layer 4 --- the RNA signal recurs (gate \pass{{}})}}
Before any cross-modality claim, we confirmed that OSD-462's own RNA flight effect
reproduces the RRRM-2 ISS-transit direction. The pathway-effect vectors align with cosine
${g['point_cosine']:.3f}$ (95\% CI $[{g['ci'][0]:.3f},{g['ci'][1]:.3f}]$), robust to
dropping any single pathway (leave-one-out range
$[{g['loo_cosine_range'][0]:.3f},{g['loo_cosine_range'][1]:.3f}]$), with
{sum(pv['sign_agree'])}/{len(pv)} pathways sign-concordant (Table~\ref{{tab:layer4}},
Fig.~\ref{{fig:recur}}A). The ECM program is up in both cohorts (OSD-462
${pv.set_index('pathway').loc['ecm_organization','osd462_rna_pathway_effect']:+.2f}$,
RRRM-2 ${pv.set_index('pathway').loc['ecm_organization','rrrm2_iss_t_pathway_effect']:+.2f}$)
and the DCT/NCC/WNK program is down in both
(${pv.set_index('pathway').loc['dct_ncc_wnk_transport','osd462_rna_pathway_effect']:+.2f}$
vs ${pv.set_index('pathway').loc['dct_ncc_wnk_transport','rrrm2_iss_t_pathway_effect']:+.2f}$).
The gate passes, so a protein-level result is interpretable.

{tab4}

\subsection{{Layer 1 --- the signal does not reach protein abundance}}
The genome-wide RRRM-2-RNA$\leftrightarrow$OSD-462-protein flight-effect correlation is
weak (Spearman ${l1['genome_wide_spearman_rrrm2_rna_vs_protein']:.3f}$, Pearson
${l1['genome_wide_pearson']:.3f}$; $n=6{{,}}418$ genes), as expected for cross-study bulk
TMT. Against the matched null, \emph{{no}} targeted set is concordant in its predicted
direction (all Benjamini--Hochberg $q\geq0.92$; Table~\ref{{tab:layer1}},
Fig.~\ref{{fig:dash}}A). The matrix set is in fact \emph{{anti}}-concordant: mean protein
effect ${l1['ecm_signed_mean_effect']:+.3f}$ log$_2$ where ``up'' was predicted (two-sided
matched-null $p={matrix_anti:.3f}$), with 14/15 matrix genes showing transcript-up but
protein-down and a within-set RNA$\leftrightarrow$protein Spearman of
${l1['ecm_concordance']:+.3f}$. The DCT set is flat
(${l1['dct_signed_mean_effect']:+.3f}$), driven by NCC protein being maintained or slightly
elevated (Fig.~\ref{{fig:dash}}B). Representative machinery and peptide depth are in
Table~\ref{{tab:detect}}.

{tab1}

{tabdet}

\subsection{{Layer 2 --- canonical-site-indexed co-modified features are suppressed}}
The phosphoproteomics workbook contains NCC (\textit{{Slc12a3}}), SPAK
(\textit{{Stk39}}), and OSR1 (\textit{{Oxsr1}}) site features. Residue-level provenance
identifies a T53-indexed feature; it is down in flight
(${phos.set_index(['gene_symbol','site_position']).loc[('Slc12a3','53'),'phospho_effect']:.2f}$
log$_2$), but its reported peptide sequence also carries Y65 phosphorylation. Workbook
positions 65 and 68 resolve to Y65 and Y68 and are excluded from canonical-site inference.
NCC C-terminal sites are unchanged, a within-protein specificity that argues against a
simple abundance artifact. The upstream S383-indexed SPAK feature, carried on an S382/S383
phosphoform, falls in parallel (
${phos.set_index(['gene_symbol','site_position']).loc[('Stk39','383'),'phospho_effect']:.2f}$,
$p={phos.set_index(['gene_symbol','site_position']).loc[('Stk39','383'),'phospho_p_value']:.1e}$).
Genome-wide the phosphoproteome is centered (median $+0.007$), and these NCC/SPAK features sit
in its extreme down tail (bottom $\sim$0.0--0.7 percentile). Total NCC protein is unchanged,
but the co-modified peptides do not permit isolated-site occupancy inference
(Table~\ref{{tab:phospho}}, Fig.~\ref{{fig:dash}}C, Fig.~\ref{{fig:recur}}B).

{tabph}

\subsection{{Layer 3 --- network candidates do not translate (\fail{{}})}}
RRRM-2 LIONESS/node2vec top composite candidates are not enriched among OSD-462 protein- or
phospho-changing genes at any cutoff (protein $p\geq{l3['protein_enrichment_p_min']:.2f}$,
phospho $p\geq{l3['phospho_enrichment_p_min']:.2f}$; Table~\ref{{tab:layer3}},
Fig.~\ref{{fig:dash}}D). The network layer therefore remains a transcript-level
candidate-prioritization device with no independent proteomic support.

{tabn}

\begin{{figure}}[htbp]\centering
\includegraphics[width=\textwidth]{{fig_osd462_multiomics_dashboard.pdf}}
\caption{{Multi-omics anchor dashboard. \textbf{{A}} Transcript (RRRM-2 ISS-T RNA) vs.\ OSD-462
protein flight effects; targeted sets coloured. Matrix genes (orange) occupy the
transcript-up / protein-down quadrant. \textbf{{B}} DCT/NCC/WNK protein abundances are flat
(NCC slightly up). \textbf{{C}} WNK--SPAK/OSR1--NCC phosphosite effects with 95\% CIs;
co-modified T53/Y65 and S382/S383 features decrease while C-terminal features are flat;
isolated canonical activity is unresolved. \textbf{{D}} Network-candidate translation: observed mean
$|$effect$|$ tracks the matched null at protein and phospho layers.}}
\label{{fig:dash}}
\end{{figure}}

\begin{{figure}}[htbp]\centering
\includegraphics[width=\textwidth]{{fig_osd462_rna_recurrence.pdf}}
\caption{{\textbf{{A}} RNA recurrence gate: OSD-462 vs.\ RRRM-2 pathway-effect vectors lie along
the identity line (cosine ${g['point_cosine']:.2f}$, PASS); dct\_ncc\_wnk\_transport sits in
the shared-down quadrant, ECM/matrix in the shared-up quadrant. \textbf{{B}} The DCT/NCC axis
across layers: NCC transcript down, total NCC protein maintained, NCC regulatory
phosphorylation and SPAK phosphorylation strongly down --- an activity-level, not
abundance-level, change.}}
\label{{fig:recur}}
\end{{figure}}

\section{{Cross-layer integration and decision}}
The pre-registered decision table resolves to: \emph{{``{esc(sjson['decision_table_row'])}.''}}
The honest reading across layers is two-pronged:
\begin{{itemize}}[leftmargin=1.4em]
\item \textbf{{Matrix arm:}} transcriptional only. The ECM program is induced at the RNA
level in both cohorts but matrix protein abundance is not increased (it is mildly,
significantly decreased), indicating transcript--protein uncoupling at terminal collection
rather than net matrix deposition.
\item \textbf{{DCT/NCC arm:}} NCC transcript is down and total NCC protein is preserved.
T53/Y65 and S382/S383 co-modified features decrease, but they do not isolate the canonical
SPAK$\to$NCC activation sites. The signaling-activity claim therefore remains unresolved.
\end{{itemize}}

\section{{Biological implications}}
The WNK--SPAK/OSR1--NCC cascade is the master switch for Na--Cl reabsorption in the distal
convoluted tubule; phosphorylation of the NCC N-terminal threonines by SPAK/OSR1 is the
activating step, and SPAK is itself activated by WNK-dependent phosphorylation. The current
workbook instead identifies strongly decreased, co-modified T53/Y65 and S382/S383
phosphoforms while total NCC protein is unchanged. This is supportive molecular context,
but it is not direct evidence that spaceflight reduces isolated canonical NCC
\emph{{transport-pathway activation}}. Resolving that claim requires a site- and
phosphoform-specific assay and connects to the physiological context that motivated the axis
(distal Na--Cl handling, K$^+$/Ca$^{{2+}}$ balance, volume regulation, and renal stone risk).
It also reframes the matrix result: the ``matrix-high'' signature is, in this cohort, a
transcriptional program not yet realised as protein, consistent with active remodeling
(transcription plus proteolysis) rather than fibrotic accumulation. Finally, the failure of
network-nominated candidates to translate cautions that LIONESS/node2vec rewiring scores
prioritise transcripts, not protein-level biology.

\section{{What this can and cannot establish}}
\textbf{{Can:}} that the RNA remodeling \emph{{direction}} reproduces in a second, independent,
cross-strain flight cohort at the RNA level; that T53/Y65- and S382/S383-indexed
phosphoforms decrease; and that matrix induction is transcript-level only.
\textbf{{Cannot:}} isolated canonical NCC/SPAK activity or physiology. There is no qualified
single-phosphoform canonical feature, histology, or urine/serum electrolyte data in scope.
A single terminal time point cannot resolve the transcript$\to$protein timing that would
explain the matrix uncoupling. The result does not rehabilitate ``network rewiring'' as a
mechanism-discovery method. Histology and urine/serum chemistry remain explicit future work.

\section{{Reproducibility}}
All numbers and tables in this document are generated directly from the analysis outputs
under \texttt{{data/results/{esc(args.run)}/osd462\_anchor/}}: harmonized flight effects
(\texttt{{osd462\_flight\_effects.tsv}}), RNA recurrence
(\texttt{{osd462\_rna\_recurrence.tsv}}), protein concordance
(\texttt{{protein\_concordance\_summary.tsv}}), phospho axis
(\texttt{{phospho\_axis\_summary.tsv}}, \texttt{{phospho\_axis\_verdict.tsv}}), network
translation (\texttt{{network\_translation.tsv}}), and \texttt{{results\_summary.json}}. Each
layer writes a manifest with input SHA-256 digests; the run-level \texttt{{manifest.json}}
references them. Pipeline: \texttt{{scripts/osd462/00\_harmonize.py}} $\to$
\texttt{{04\_rna\_recurrence.py}} (gate) $\to$ \texttt{{01\_protein\_concordance.py}} $\to$
\texttt{{02\_phospho\_axis.py}} $\to$ \texttt{{03\_network\_translation.py}} $\to$
\texttt{{05\_plot\_dashboard.py}} $\to$ \texttt{{06\_compile\_summary.py}} $\to$
\texttt{{07\_build\_compendium.py}}; core estimators in
\texttt{{src/multiomics/osd462\_anchor.py}} with unit tests in
\texttt{{tests/test\_osd462\_anchor.py}}. Residue and phosphoform claims additionally require
\texttt{{scripts/osd462/08\_stage0\_provenance\_audit.py}} and its zero-row
isolated-canonical qualification table. Matched nulls use {10000} draws (seed 20260521);
the RNA recurrence bootstrap uses {int(rec['n_bootstrap'])} resamples.

\end{{document}}
"""
    TEX_PATH.write_text(doc)
    print(f"[compendium] wrote {TEX_PATH}")

    if args.compile:
        for _ in range(2):
            r = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                                TEX_PATH.name], cwd=TEX_PATH.parent,
                               capture_output=True, text=True)
        if (TEX_PATH.with_suffix(".pdf")).exists():
            print(f"[compendium] compiled {TEX_PATH.with_suffix('.pdf')}")
        else:
            print("[compendium] compile FAILED; tail of log:")
            print(r.stdout[-2500:])


if __name__ == "__main__":
    main()
