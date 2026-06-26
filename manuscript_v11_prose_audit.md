# Prose & Tone Audit — `manuscript_v11.tex`

**Manuscript:** `latex_paper/manuscript_v11.tex` (833 lines; newest, Jun 18 2026)
**Benchmark corpus:** `transcriptomic_texts/` — space-biology, renal-physiology, and RNA-seq/methods literature
**Date:** 2026-06-19
**Scope:** Where the manuscript's *way of presenting findings and methods* departs from the register of the surrounding literature — editorializing, authorial affect, promotional method language, rhetorical/colloquial constructions. Numbers, citations, and scientific content were **not** changed.

---

## 1. Bottom line

The manuscript is, overall, written at or above field standard: dense, careful, heavily hedged, and methodologically explicit. The tone problems are **localized and few** — they cluster almost entirely in the Results/Discussion sentences that frame the surprising DCT2 result, exactly where you flagged them. There is no systemic register problem; there is a recurring **habit of telling the reader how to feel about a result** ("unexpectedly," "provocative," "the informative part," "deliberately strong") rather than stating the result and letting the statistics carry it.

Twenty-one passages are flagged below. Six are **High** severity (clear departures from field register), nine **Medium**, six **Low**. All have been corrected in the manuscript; see §6 for the applied edits.

---

## 2. The field-norm benchmark (why these are flagged)

I characterized the register of the surrounding literature directly from the reference PDFs, with emphasis on **how authors mark a surprising or counterintuitive result**. Word-frequency across three representative papers (the central citation plus two space-biology/biology papers):

| Marker | Siew 2024 *Nat Commun* (16,145 w) | npj Microgravity review (8,780 w) | BMC Biol 2025 (5,388 w) |
|---|---|---|---|
| unexpectedly | 0 | 0 | 0 |
| surprisingly | 0 | 0 | 0 |
| strikingly | 0 | 0 | 0 |
| provocative | 0 | 0 | 0 |
| notably | 1 | 1 | 3 |
| interestingly | 1 | 3 | 0 |

**Reading of the norm.** Published work in this subfield essentially *never* editorializes surprise with "unexpectedly/surprisingly/strikingly/provocative." When a result is flagged at all, it is with a sparing, low-key marker ("notably," "interestingly") — on the order of one per several thousand words — and findings are otherwise presented with neutral observational verbs: *"we observed incidences of thrombotic microangiopathy…,"* *"…again consistent with remodelling"* (Siew 2024). The convention is to report the effect and its statistic and let the reader judge salience.

Against that baseline, the manuscript's two uses of **"unexpectedly,"** its **"the more provocative half,"** and its **"deliberately strong negative"** are genuine outliers — not because the words are wrong in English, but because they import an authorial-affect register the target journals don't use. This is the concrete basis for your instinct that the flagged sentences read as "unprofessional."

---

## 3. High-severity flags (clear departures; fix strongly advised)

### A1 — line 415, Results §"Flight-suppressed phosphosites…" (lead sentence) — *your example*
> "…they concentrate in the parent genes that define distal-nephron identity, and—**unexpectedly**—the signal that holds under the most stringent observability-matched permutation is the DCT2-leaning extreme, **not** the canonical DCT1 anchors."

**Issue.** Two problems stacked: (i) "unexpectedly" editorializes the author's surprise into a results statement; (ii) the em-dash–bracketed single adverb ("—unexpectedly—") is a dramatization device that the field register avoids in Results. The contrast "X, not Y" also reads as rhetorical emphasis.
**Field norm.** State the contrast neutrally; surprise, if relevant, belongs in the Discussion as a measured observation, not mid-Results.
**Rewrite (applied).** "…and the signal that holds under the most stringent observability-matched permutation is the DCT2-leaning extreme **rather than** the canonical DCT1 anchors."

### A2 — line 606, Conclusion (final sentence) — *your example*
> "…and would clarify the DCT2-leaning signal that, **unexpectedly**, is the most permutation-robust part of the enrichment."

**Issue.** Same editorializing adverb, now closing the paper — the most register-sensitive sentence in the manuscript.
**Rewrite (applied).** "…and would clarify the DCT2-leaning signal, **which** is the most permutation-robust part of the enrichment."

### A3 — line 555, Discussion §"The distal-nephron subtype-prior signal…"
> "The DCT2-leaning result is **the more provocative half**: the canonical NCC/SPAK/WNK anchors are DCT1-leaning, yet the DCT2-leaning bin is the one that survives observability-matched permutation."

**Issue.** "the more provocative half" is the strongest affect marker in the paper and has zero analogue in the benchmark corpus. It characterizes the result for the reader instead of stating it.
**Rewrite (applied).** "**Although** the canonical NCC/SPAK/WNK anchors are DCT1-leaning, it is the DCT2-leaning bin that survives observability-matched permutation."

### A4 — line 500, Results §"Network, cross-species…"
> "This is a **deliberately strong negative**—it shows that the phosphoproteomic signal is not merely recovering network-hub genes, and that the network method did not **manufacture the result** by routing through highly connected nodes."

**Issue.** Two markers: "deliberately strong negative" is self-praise of one's own design; "manufacture the result" is loaded, faintly defensive, and informal. Reviewers read self-congratulatory method language as a yellow flag.
**Rewrite (applied).** "This **negative control indicates** that the phosphoproteomic signal is not merely recovering network-hub genes, and that the network method did not **produce** the result by routing through highly connected nodes."

### A5 — line 318, Results §"OSD-462 shows RNA-protein mismatch"
> "**The consequence is concrete: a reader with only the RNA data would infer** that matrix-remodeling proteins accumulate in the spaceflight kidney, which the protein layer does not show."

**Issue.** Addressing a hypothetical "reader" is an essayistic/rhetorical move, not a Results register; "The consequence is concrete" is informal emphasis.
**Rewrite (applied).** "**The practical consequence is specific: RNA data alone would suggest** that matrix-remodeling proteins accumulate in the spaceflight kidney, which the protein layer does not show."

### A6 — line 557, Discussion
> "The calcium machinery is **the weakest link**: \gene{Trpv5} and \gene{Calb1} sit in the bin by expression but contribute no suppressed phosphosites…"

**Issue.** "the weakest link" is an idiom/colloquialism out of register for a Discussion of evidence strength.
**Rewrite (applied).** "The calcium machinery is **the least-supported component**: …"

---

## 4. Medium-severity flags (register slips; fix recommended)

### B1 — line 423, Results
> "…at parent-gene level (OR 1.82, q=…), and—**notably**—under the same matched parent-gene permutation that DCT1 fails…"

**Issue.** "notably" is itself within field frequency, but the em-dash dramatization "—notably—" over-emphasizes. The OR and *p* already mark salience.
**Rewrite (applied).** Remove the em-dash insertion: "…at parent-gene level (OR 1.82, q=…), and under the same matched parent-gene permutation that DCT1 fails…"

### B2 — line 576, Discussion §"Remodeling versus transporter suppression"
> "…formalize an observation that is **biologically interesting on its own**: animals with higher endothelial/remodeling scores tended to have lower NCC regulatory phosphorylation."

**Issue.** "interesting on its own" tells the reader the observation is interesting.
**Rewrite (applied).** "…formalize **a consistent observation**: animals with higher endothelial/remodeling scores tended to have lower NCC regulatory phosphorylation."

### B3 — line 529, Discussion §"Main interpretation"
> "**The divergence is the informative part**—the RNA and protein layers systematically diverge while phosphorylation aligns…"

**Issue.** Part of a repeated "informative" framing device (see B4–B5). Labeling your own finding "the informative part" is editorial.
**Rewrite (applied).** "**This divergence is itself the central finding:** the RNA and protein layers systematically diverge while phosphorylation aligns…"

### B4 — line 318, Results
> "Against that baseline **the informative result is pathway-specific**…"

**Issue.** "informative" used as an editorial signpost (one of four occurrences across the paper).
**Rewrite (applied).** "Against that baseline **the substantive result is pathway-specific**…"

### B5 — line 487, Results §"Public perturbation references…"
> "Two references that might have named an upstream mechanism were coverage-limited, **but informatively so**. … **The informative result is the negative one**…"

**Issue.** Two more "informative" signposts in two sentences.
**Rewrite (applied).** "…were coverage-limited, **but the gap is itself informative**. … **The salient result is the negative one**…" (varies the wording and breaks the repetition).

### B6 — line 500, Results
> "Three further comparisons were **factual checks**: live-return vectors attenuated…"

**Issue.** "factual checks" is vague/idiosyncratic (all analyses are presumably factual). The sentence's own closing clause calls them "boundary conditions."
**Rewrite (applied).** "Three further comparisons served as **boundary-condition checks**: …"

### B7 — line 340, Results
> "An orthogonal method sharing no inputs with the targeted analysis thus reaches the same WNK axis, **extending the result rather than restating it**."

**Issue.** "extending … rather than restating it" editorializes the contribution.
**Rewrite (applied).** "…thus reaches the same WNK axis, **providing independent support**."

### B8 — line 423, Results
> "…which is what lets a whole-kidney phosphoproteome **speak to** a specific nephron segment."

**Issue.** "speak to" is colloquial; given the whole-kidney caveat, a precise verb is also safer.
**Rewrite (applied).** "…which is what allows a whole-kidney phosphoproteome to **implicate** a specific nephron segment." (matches the title verb, "implicate").

### B9 — line 606, Conclusion (opening sentence)
> "This analysis **changes how a bulk RNA spaceflight kidney signal should be read**."

**Issue.** Grand opening claim; field register prefers a measured statement of what the results indicate.
**Rewrite (applied).** "These results **refine how a bulk RNA spaceflight kidney signal should be interpreted**."

---

## 5. Low-severity flags (optional polish; applied where unobtrusive)

| # | Line | Text | Issue | Action |
|---|---|---|---|---|
| C1 | 423 | "The DCT2-leaning bin **is the one that** strengthens where DCT1 attenuates." | "is the one that" is conversational | → "The DCT2-leaning bin **strengthens** where DCT1 attenuates." (applied) |
| C2 | 555 | "…the DCT2-leaning bin **is the one that survives**…" | same construction | folded into A3 rewrite (applied) |
| C3 | 559 | "A candidate upstream signal is now testable **rather than simply open**." | mild editorial flourish | left as-is (borderline; acceptable) — *no change* |
| C4 | 529 | "…the present analysis **adds what that endpoint alone could not show**." | faintly promotional | left as-is (defensible contribution statement) — *no change* |
| C5 | 296 | "This whole-panel agreement **is driven primarily by** the remodeling pathways." | fine; flagged only for completeness | *no change* |
| C6 | 340 | "**The strongest independent corroboration came from** the unbiased kinome-wide screen." | mild superlative; acceptable as it is quantified immediately after | *no change* |

---

## 6. Edits applied to `manuscript_v11.tex`

All rewrites in §3–§5 (those marked "applied") were made directly in `latex_paper/manuscript_v11.tex`. Each edit:

- changes only the editorializing/colloquial wording;
- preserves every number, *p*/*q* value, odds ratio, `\citep{}`, `\gene{}` macro, math, and em-dash markup elsewhere in the sentence;
- preserves sentence meaning and the scientific claim.

Net effect on the paper: "unexpectedly" → 0 occurrences (was 2); "provocative," "deliberately strong," "manufacture," "weakest link," "factual checks" → removed; the "informative" signpost reduced from 4 to 1; two em-dash adverb dramatizations removed; one reader-addressing construction neutralized. The 14 affected sentences keep their meaning and citations.

A self-documenting list of the exact before/after strings is in §3–§4 above.

---

## 7. Patterns worth keeping in mind (for future drafts)

1. **Let the statistic mark salience.** When an OR and a permutation *p* are already in the sentence, "unexpectedly/notably/provocative" is redundant and shifts register. The reference corpus does this consistently.
2. **Surprise belongs in the Discussion, stated plainly.** "Counterintuitively, the bin that survives permutation is the DCT2-leaning one" is acceptable *once* in Discussion; the same word twice (Results lead + Conclusion close) reads as a tic.
3. **Avoid characterizing your own design** ("deliberately strong," "the informative part," "manufacture the result"). State what the control shows; reviewers supply the adjectives.
4. **Don't address the reader.** "A reader with only the RNA data would infer…" → "RNA data alone would suggest…".
5. **Em-dash single-adverb insertions** ("—unexpectedly—", "—notably—") are an emphasis device; the field uses them rarely and not in Results.

None of this touches the science, the hedging discipline (which is excellent), or the methods register (which is strong throughout).
