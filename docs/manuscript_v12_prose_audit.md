# Final-draft prose audit — `manuscript_v12` (version / pipeline "artifact" language)

**Manuscript:** `latex_paper/manuscript_v12.tex` (875 lines; = the Google Doc `manuscript_v12`, verified verbatim on the flagged passages, and = `manuscript_v12.docx`).
**Date:** 2026-07-23
**Scope:** A *different* axis from the June 2026 audit (`manuscript_v11_prose_audit.md`), which removed editorializing/affect ("unexpectedly," "provocative") — all of those fixes carried into v12 cleanly. This pass targets the register problem you flagged now: places where the paper **refers to itself as a software version/draft** or **describes its own outputs as engineering artifacts** ("emitted with the v11 run artifacts," "v11 baseline artifacts," "the project pipeline," "sheets," "Draft v12"), and a few **clumsy, notes-like sentences** (the RRRM-2 "design is the richest" example). These read like an internal build log, not a final manuscript. **No numbers, statistics, citations, or scientific claims are changed** in any rewrite below.

The guiding principle, same as before: a published paper describes *what was done and where the data live*, never *which internal run produced the file*. The reference for "how to write this" is the paper's own strong Methods register (e.g. §2.5, §2.6) and the Siew 2024 / npj Microgravity house style — neither ever names a run folder or a `.py` module inside a sentence.

---

## 1. Bottom line

The science and hedging are unchanged and strong. The problem is a thin but pervasive layer of **provenance/tooling vocabulary** that leaked from the repo into the prose. It clusters in three predictable places — the title page, the Methods reproducibility sentences, and the two supplementary-pointer sentences you quoted — plus one genuinely clumsy cohort-description sentence.

**26 items** are flagged, grouped by type. The single highest-impact fix is the smallest: the **title page literally says "Draft v12"** (§3.1). After that, the recurring tells are the word **"emitted"** (3×), the token **"v11"** in visible prose (5 sentences + 1 table caption), two leaked **`run_20260709`** IDs, inline **`src/v11/*.py`** paths mid-sentence, and the **"the project pipeline / repository preprocessing pipeline / run artifacts"** family. Fixing the ~15 High/Medium items removes essentially all of the "prototype" feel.

A ready-to-run find/replace table is in §7, and a note on the Google-Doc-vs-source workflow is in §8.

---

## 2. Severity key

- **High** — visibly marks the paper as a draft/build artifact to any reader (title page, "emitted," "v11 analysis," leaked run IDs). Fix strongly advised.
- **Medium** — tooling register a reviewer will notice as informal; not a draft-marker per se. Fix recommended.
- **Low** — consistency/polish; optional.

---

## 3. HIGH — self-references to a draft version or internal run

### 3.1 — line 56, title page — *the most visible one*
> `\date{Draft v12, July 2026}`

**Issue.** The title page announces the document as a numbered internal draft. This is the first thing a reader (or reviewer, or the Google-Doc viewer) sees, and it undercuts "final."
**Rewrite.** `\date{July 2026}` — or, if you want to signal preprint status, `\date{Preprint, July 2026}`.

### 3.2 — line 672, Table S1 caption
> "Multiple-testing families used in the **main v11 analysis**."

**Issue.** "v11 analysis" names an internal manuscript version as if it were the study.
**Rewrite.** "Multiple-testing families used in the **main analysis**."

### 3.3 — line 198, Methods §2.3 — *your example*
> "…the frozen member lists, pathway-family labels, leave-one-family recurrence check, and paired pathway bootstrap **are emitted with the v11 run artifacts**."

**Issue.** Triple tell: "emitted" (a program writes files), "v11 run artifacts" (internal version + build vocabulary), and "frozen member lists" (engineering jargon for a pinned list).
**Rewrite.** "…the member lists, pathway-family labels, leave-one-family recurrence check, and paired pathway bootstrap **are provided in the companion results compendium and project repository**."
*(Point it at whatever actually holds them — Supplementary Information if you add a table, otherwise the compendium/repo. "frozen" → dropped; "emitted" → "provided.")*

### 3.4 — line 240, Methods §2.4 — *your example*
> "**Raw-workbook** TMT QC summaries of channel medians and missingness **are emitted** for protein, single-position phosphosite, and composite phosphosite **sheets with the v11 baseline artifacts**."

**Issue.** "Raw-workbook," "emitted," "sheets," and "v11 baseline artifacts" together read as a spreadsheet-export note.
**Rewrite.** "**TMT quality-control summaries of channel medians and missingness for the protein, single-position phosphosite, and composite phosphosite tables are reported in the Supplementary Information.**"

### 3.5 — line 234, Methods §2.4
> "(median and maximum per-channel missing fraction 0.0 across plexes and conditions **in the emitted QC summary**)"

**Issue.** "emitted QC summary."
**Rewrite.** "(median and maximum per-channel missing fraction 0.0 across plexes and conditions; **Supplementary Information**)".

### 3.6 — line 474, Results §"…subtype-prior extremes" — leaked run ID + code variable
> "Because **`dct1_enrichment_score`** is a raw mean-expression difference… we re-derived the DCT2-leaning bin under three alternative subtype scores and re-ran the parent-gene enrichment and matched permutation **(run_20260709; Table 6, lower block)**."

**Issue.** A monospace code variable name (`dct1_enrichment_score`) and a raw internal run ID (`run_20260709`) inside a Results sentence.
**Rewrite.** "Because **the DCT1 enrichment score** is a raw mean-expression difference… we re-derived the DCT2-leaning bin under three alternative subtype scores and re-ran the parent-gene enrichment and matched permutation **(Table 6, lower block)**."

### 3.7 — line 492, Table 6 sub-caption — leaked run ID
> "*DCT2-leaning bin re-scored under alternative subtype scores **(run_20260709)**; DCT2 OR is Fisher…*"

**Issue.** Same `run_20260709` in a table caption.
**Rewrite.** Delete "(run_20260709)": "*DCT2-leaning bin re-scored under alternative subtype scores; DCT2 OR is Fisher…*"

### 3.8 — line 244, Methods §2.4 — inline module paths mid-sentence
> "Detectability robustness was assessed by extending the matched stratum with a per-gene missing-fraction bin and rerunning the same tests **(`src/v11/rna_protein_propagation.py`, `src/v11/observability_audit.py`)**."

**Issue.** Two source-file paths embedded in a Methods sentence. Method belongs in prose; file paths belong in Data & Code Availability (where they are already listed).
**Rewrite.** End the sentence at "…rerunning the same tests." (Optionally "…rerunning the same tests (Data and Code Availability).")

### 3.9 — line 638, Discussion §"Open novelty extensions" — inline repo doc path
> "Both, with their pre-registered nulls, are documented in **`docs/v11_layer_specificity_execution_summary_2026-06-07.md`**."

**Issue.** A repo markdown path (with an internal date-stamped filename) dropped into Discussion prose.
**Rewrite.** "Both, with their pre-registered nulls, are **specified in the project repository (Data and Code Availability)**."

### 3.10 — line 651, Author Contributions
> "…maintained the source code **and run artifacts**."

**Issue.** "run artifacts."
**Rewrite.** "…maintained the **analysis code and computational pipeline**." (or simply "…maintained the source code.")

---

## 4. MEDIUM — software/pipeline register in Methods

### 4.1 — line 336, Methods §2.10 (reproducibility sentence)
> "Analyses were implemented in **the project Python/R pipeline** with versioned configuration files, **run manifests**, tests, and scripted figure generation. Detailed **run identifiers** and secondary outputs are **archived in the repository** and companion results compendium."

**Issue.** Reads like a repo README: "the project pipeline," "run manifests," "run identifiers," "archived in the repository." The reproducibility content is worth keeping; the devops vocabulary is not.
**Rewrite.** "Analyses were implemented in **Python and R** with versioned configuration files and scripted figure generation, under automated tests. **Run identifiers and secondary outputs are available in the project repository and the companion results compendium (Data and Code Availability).**"

### 4.2 — line 196, Methods §2.3
> "Primary RRRM-2 analyses used **the project residualized expression matrix** and matched metadata **generated by the repository preprocessing pipeline**."

**Issue.** "the project … matrix," "the repository preprocessing pipeline."
**Rewrite.** "Primary RRRM-2 analyses used **a residualized expression matrix** and matched metadata **produced by a standardized preprocessing workflow** (Data and Code Availability)."

### 4.3 — line 300, Methods §2.7
> "KSEA was implemented **in the project pipeline** against the atlas rather than via a packaged tool."

**Issue.** "in the project pipeline" adds nothing and lowers register.
**Rewrite.** "KSEA was implemented **directly against the atlas** rather than via a packaged tool."

### 4.4 — "emitted" residual sweep
"emitted" appears 3× (lines 198, 234, 240 — all handled in §3.3–3.5). After those edits, confirm zero remaining: the word should not appear in the final manuscript. It is the single most reliable "software wrote this" tell.

---

## 5. MEDIUM/LOW — clumsy or notes-like phrasing

### 5.1 — line 168, Methods §2.1 — *your example (the RRRM-2 sentence)*
> "RRRM-2/OSD-771 served as the primary RNA-seq cohort **because its design is the richest**: it spans young and old female C57BL/6NTac mice **across ISS terminal dissection, live animal return, flight, habitat ground control, vivarium, and basal contexts**, supporting within-mission flight-versus-habitat-control contrasts and age stratification **that the other cohorts do not**."

**Issue.** "its design is the richest" is a colloquial value-judgment; the middle is an ungoverned pile-up of six comma-separated contexts that reads like a data-dictionary dump; "that the other cohorts do not" is clipped.
**Rewrite (journal register, facts unchanged).** "RRRM-2/OSD-771 served as the primary RNA-seq cohort because it has **the most comprehensive design**: it profiles young and old female C57BL/6NTac mice across **multiple mission contexts — ISS terminal dissection, live-animal return, flight, habitat ground control, vivarium, and basal —** and so **uniquely supports** within-mission flight-versus-habitat-control contrasts together with **age stratification not available in the other cohorts**."

### 5.2 — line 55, title page — affiliation
> `\affil{Independent research manuscript}`

**Issue.** "Independent research manuscript" is not an affiliation; it describes the document, not the author.
**Rewrite.** `\affil{Independent researcher}` (or your actual affiliation / city).

### 5.3 — lines 336, 343, 442, 512, 541, 545, 558, 657, 755 — "compendium" naming is inconsistent
Nine references to the supplementary document, in three different forms: **"companion results compendium"** (5×), **"companion compendium"** (3×), **"the compendium"** (1×).

**Issue.** A final paper names its supplement once and consistently. Three spellings read as unfinished.
**Rewrite.** Pick one label and use it everywhere. Recommended: **"the companion results compendium"** on first mention, **"the compendium"** thereafter — or, if it is really your supplement, standardize to **"Supplementary Information / Supplementary Data."**

### 5.4 — lines 328, 558, 683 — `\code{sig_info}` in prose; line 529 — `\code{brms}`
Monospace code identifiers (`sig_info`, a CMap metadata field; `brms`, an R package) set inside sentences.

**Issue.** Minor register slip; monospace tokens read as code, not method.
**Rewrite.** `sig_info` → "the CMap signature metadata" (keep the term once in parentheses if needed). `brms` → "a full Bayesian multilevel model" (drop the package name, or keep it in parentheses: "…rather than a full structural equation or Bayesian multilevel (brms) model").

### 5.5 — lines 306, 529 — "SEM" unexpanded
"SEM" is used at line 306 and 529 without expansion.
**Rewrite.** First use (line 306): "…rather than a full **structural equation model (SEM)**." Then "SEM" is fine.

---

## 6. Source-only housekeeping (invisible in the Google Doc, but marks the draft)

These do **not** appear in the rendered PDF/Doc, so they are not why the paper *reads* unfinished — but they are draft-markers in the source worth clearing before you call v12 final.

- **line 168** — trailing LaTeX comment: `% CONFIRM OSD-102/OSD-163 accessions, strains, and that both are kidney RNA-seq before finalizing.` This is an open TODO against a factual claim in Table 1 (OSD-102 = F/C57BL/6J; OSD-163 = BALB/c). Resolve the fact, then delete the comment.
- **Figure filenames** — many figures are `v11_*.pdf` (e.g. `v11_cross_layer_heatmap.pdf`, `v11_dct1_vs_dct2_evidence_ladder.pdf`) and `\graphicspath` points at `run_20260526_v11_dct1_phospho_mediation/…`. Invisible to readers, but if you ever share the source or a figure bundle, the "v11/run_2026…" names reappear. Optional: rename on final.
- **`src/v11/` code paths in Data & Code Availability (line 657)** — these are *real* file locations, so listing them there is legitimate (unlike the inline uses in §3.8–3.9). The only oddity is the "v11" directory name in a v12 paper. This is a repo-rename decision, not a prose fix — leave as-is unless you want to rename the directory and update the paths.

---

## 7. Find / replace quick table (visible-prose items)

| # | Sev | Find (current) | Replace (suggested) | Line |
|---|-----|----------------|---------------------|------|
| 1 | High | `Draft v12, July 2026` | `July 2026` | 56 |
| 2 | High | `the main v11 analysis` | `the main analysis` | 672 |
| 3 | High | `are emitted with the v11 run artifacts` | `are provided in the companion results compendium and project repository` | 198 |
| 4 | High | `frozen member lists` | `member lists` | 198 |
| 5 | High | `Raw-workbook TMT QC summaries … sheets with the v11 baseline artifacts` | `TMT quality-control summaries … tables are reported in the Supplementary Information` | 240 |
| 6 | High | `in the emitted QC summary` | `; Supplementary Information` | 234 |
| 7 | High | `` `dct1_enrichment_score` `` | `the DCT1 enrichment score` | 474 |
| 8 | High | `(run_20260709; Table…` | `(Table…` | 474 |
| 9 | High | `(run_20260709)` | *(delete)* | 492 |
| 10 | High | inline `src/v11/…py, …py` in §2.4 | *(delete; already in Data & Code Availability)* | 244 |
| 11 | High | `documented in docs/v11_layer_specificity_execution_summary_2026-06-07.md` | `specified in the project repository (Data and Code Availability)` | 638 |
| 12 | High | `and run artifacts` (Author Contributions) | `and computational pipeline` | 651 |
| 13 | Med | `the project Python/R pipeline … run manifests … archived in the repository` | `Python and R … under automated tests … available in the project repository` | 336 |
| 14 | Med | `the project residualized expression matrix … the repository preprocessing pipeline` | `a residualized expression matrix … a standardized preprocessing workflow` | 196 |
| 15 | Med | `implemented in the project pipeline against the atlas` | `implemented directly against the atlas` | 300 |
| 16 | Med | `because its design is the richest: … contexts, … that the other cohorts do not` | see §5.1 rewrite | 168 |
| 17 | Low | `Independent research manuscript` | `Independent researcher` | 55 |
| 18 | Low | `companion compendium` / `the compendium` (mixed) | standardize to one label | 343,545,558,755,541 |
| 19 | Low | `\code{sig_info}` / `\code{brms}` | de-monospace / expand | 328,529,558,683 |
| 20 | Low | `SEM` (first use) | `structural equation model (SEM)` | 306 |

---

## 8. Workflow note (important)

The Google Doc `manuscript_v12` is an **uploaded copy** of the `.docx`, which was built from `manuscript_v12.tex`. Editing the Doc in Drive and editing the source are two different things:

- **To fix the paper for real,** edit `latex_paper/manuscript_v12.tex`, rebuild the PDF/`.docx`, and re-upload — the Google Doc will **not** update on its own.
- Editing the Doc text directly in Drive would fix only that copy and would drift from the source.

I can apply every §3–§5 rewrite directly to `manuscript_v12.tex` (preserving all numbers, `\citep{}`, `\gene{}`, and math, exactly as the June audit did) on your go-ahead. The only judgment call is where the three "provided in…" pointers (§3.3, §3.4, §3.5) should aim — Supplementary Information vs. the companion compendium vs. the repository — since that depends on which of those you're treating as the paper's formal supplement.
