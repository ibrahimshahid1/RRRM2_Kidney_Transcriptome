# Next-directions strategy: working state (checkpoint)

This folder holds the inputs to the "what next" decision taken after the G1/G2 results
(`docs/CLINICAL_AXES_PODOCYTE_SENSITIVITIES_AND_RECOVERY_PERSISTENCE_2026-10-06.md`). It is committed
so the work survives a lost session or container. Read this file first when resuming.

## What has run

Workflow `wf_1207e53e-984` ran on 2026-10-06. Ten agents were planned; three finished before the
weekly usage limit stopped the rest.

| Agent | Status | Output |
|---|---|---|
| reader: project history and workstream status | done | `reader_history.md` |
| reader: leads in the G1/G2 run outputs (with post-hoc recomputations) | done | `reader_leads.md` |
| reader: external data landscape | done | `reader_external_data.md` |
| reader: prior art / novelty | failed (usage limit) | rerun with a smaller scope |
| reader: paper readiness | failed (usage limit) | largely covered by `reader_history.md`; not rerun |
| 3 proposers + 2 judges | failed (usage limit) | replaced: the coordinator drafts, then one judge reviews |

The readers' recomputations (GC/length correction, sampling-marker covariates, control-vs-control
contrasts, OSD-457 scoring) are **post hoc and unregistered**. Their scripts and TSVs lived in the
session scratchpad and are not in the repo. Any number that goes into a paper must be preregistered and
re-run through the `rrrm2` CLI.

## Budget-aware orchestration (adopted after the 2026-10-06 overrun)

The first run started 10 agents at once. Each one re-read the repo, and each downstream agent received
about 100K characters of raw digest. That used about 1.3M subagent tokens in roughly 30 minutes, then hit
the limit and lost 7 of 10 results. Rules from now on:

1. **Commit every agent result before starting the next stage.** Session scratchpads and workflow
   journals do not survive a container reset; git does.
2. **No more than 2–3 agents in flight at once.** Run stages in sequence across turns, and check the
   result between stages.
3. **Use the smallest model that fits each job:**
   - Sonnet at medium effort for read, search and extract jobs;
   - Opus only for adversarial judging and statistical design;
   - the coordinator writes the synthesis itself, because it already holds the context.
4. **Pass compact digests, not raw dumps.** The coordinator writes a digest of at most about 3K words.
   Downstream agents get that digest plus the file paths, never the full upstream JSON.
5. **Hard scope per agent.** Name the files and sources to read, cap web searches at about 15 per agent,
   and cap output with a schema and word limits.
6. **One stage per workflow.** If a limit hits, at most one stage is lost, and it can be resumed from the
   committed files.

## Remaining stages

1. A prior-art reader (one Sonnet agent, scoped) writes `reader_priorart.md`. Commit.
2. The coordinator writes `DIRECTIONS_DRAFT.md`: ranked directions, what to stop, and sequencing.
   Commit.
3. One adversarial judge (Opus) writes `DIRECTIONS_REVIEW.md`. Commit.
4. The coordinator finalizes `docs/NEXT_DIRECTIONS_2026-10-07.md` for the owner.
