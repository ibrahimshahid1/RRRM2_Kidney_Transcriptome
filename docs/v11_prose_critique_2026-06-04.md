# v11 prose critique & edit pass — 2026-06-04

Applied the same lens you used on the abstract/intro/methods/results to the **whole** manuscript, built a critique list, and edited prose only. No statistics, citations, figures, or claims were changed. The single place numbers were touched was de-duplication: the five-cohort meta-analysis stats appeared twice in Results — the second copy (old §3.2 close) now references the first instead of restating `p`/`I²`. All those values remain in §3.1 and the abstract.

## The problem classes (your lens, restated)

1. **Reader-constricting / interpretive-directive sentences** — telling the reader what the paper *is* or what to conclude ("the backbone of the study", "that distinction matters"). That framing is the reader's job.
2. **Self-binning / reviewer-voice** — pre-sorting the work ("incremental", "the supported X is…") instead of just stating it.
3. **Connective tics** — sentences opening with *Thus,* / *Together,* / *is therefore* as a reflex.
4. **Choppy → run-on rhythm** — a short clause followed by a 60-word data-dump, especially the intro methods list.
5. **Over-explaining method scope** — "it asks X, not whether Y" meta-commentary.
6. **Redundant restatement** — same stats / same caveat repeated across sections.

## What I changed

**Abstract**
- "Together, the public data support…" → "Across these analyses, the public data support…" (de-tic).

**Introduction**
- The hypotheses→methods sentence was a flat ~60-word run-on list. Re-grouped into the three analysis levels (RNA layer / OSD-462 anchor / DCT prior) so the methods map onto the three hypotheses, with the sensitivity/context analyses as a second sentence.

**Results**
- Deleted "This cross-cohort recurrence is the transcriptomic backbone of the study." — the recurrence stands on its own.
- "Thus, the recurrent RNA component is strongest…" → stated plainly without the connective.
- "Together, the matched multi-omic data show a *broader* layer-specific structure…" → "Across the panel, … a layer-specific structure…" (de-tic, drop "broader"/"adequately").
- Old §3.2 closing paragraph restated the full five-cohort meta-analysis (`z`, `p`, `I²`) already given in §3.1. Rewrote it as a callback that keeps only the unique Stouffer `z` values and drops the duplicated `p`/`I²`.
- "The supported conclusion is group-level concordance…" → "The data therefore support group-level concordance…".
- "The paired DCT2-leaning result changed the interpretation from…" → "The paired DCT2-leaning bin was enriched as well, broadening the pattern beyond a DCT1-only signal."
- "The supported interpretation is therefore distal-nephron subtype-prior enrichment spanning…" → "The enrichment therefore spans both DCT1-high and DCT2-leaning extremes…".
- "Thus, the enrichment follows relative phosphosite ranking…" → dropped the "Thus,".
- "The result is therefore better described as subset enrichment…" → dropped the "therefore".
- "motivating exploratory covariance-decomposition summaries" → "motivating the covariance-decomposition summaries reported here" (the subsection title already carries "Exploratory"; removed the echo).

**Discussion / Conclusion**
- Deleted "That distinction matters biologically." — the biology that follows already makes the point; lead with it. ("regulated heavily" → "regulated largely".)
- Conclusion: "Together, these findings support…" → "These findings support…".

**Methods**
- "This matrix is descriptive: it asks where curated pathway signals appear…, not whether individual genes or sites are significant." → "This matrix is descriptive, summarizing where curated pathway signals appear across molecular layers rather than testing the significance of individual genes or sites."

## Deliberately left alone

- **Honest limitations and scope statements** ("rather than", "not", whole-kidney / parent-gene-prior caveats). You explicitly value the mature/honest register; these are load-bearing, not padding. I only thinned them where a sentence restated another verbatim.
- **"exploratory"** where it flags a real analysis tier (it's a methodological category, unlike "incremental"). I removed one redundant echo, kept the rest.
- Abstract and introduction openings — your resolved versions; left as is apart from the two items above.

## Verification

- `latexmk -pdf` → exit 0, 30 pages, no errors or new warnings.
- Confirmed all targeted phrases are gone and the de-duplicated `p`/`I²` values still appear in §3.1 and the abstract.
- Connective-tic counts: `Together,` 5 → 2; `is therefore` 4 → 1.
