# `rrrm2` — the experiment CLI

Every concluded methodology revision in this repository is registered as a
machine-readable experiment class, and can be re-run end to end from one
command. Results always land in that revision's own results folder.

## Install

```bash
venv/bin/python -m pip install -e .
```

The console script `rrrm2` is then on the venv's PATH. Without installing,
`venv/bin/python -m rrrm2 ...` is equivalent.

## The registry

One file per revision under `config/revisions/`. `config/revisions/_schema.yaml`
documents every field. A revision records its scientific question, the
conclusion it actually reached, its frozen spec config, its decision documents,
its results root, the historical locked run, and an ordered, dependency-aware
list of executable stages.

| Revision | WF | Status | What it concluded |
|---|---|---|---|
| `subtype-reference` | C | locked | Flight-blind DCT2/CNT (27) and ASDN (29) sets, frozen at SHA-256 `ed64f0e0…` |
| `v13-phospho` | C | locked | `claim_tier = neither`; DCT2/CNT non-evaluable, ASDN failed the selectivity gate |
| `clinical-axes` | D | complete | 3 of 4 axes no-go; barrier axis rejects in the opposite direction; podocyte tier leads but the structural control ranks 4th |
| `osd462-stage0` | B | complete | Condition perfectly aliased with reporter block; zero isolated canonical NCC/SPAK phosphoforms |
| `v13-compartment` | C | complete | Parent-protein annotation enrichment only; does not localise cell of origin |
| `grey60` | A | blocked | No-go; rescue explicitly retired |
| `contrast-vectors` | A | retired | Never carried to a claim; gated on Guardrail A |
| `network-rewiring` | A | retired | Negative — network features did not beat expression baselines |
| `v11-core` | – | retired | Superseded and retracted by `osd462-stage0` |

## Commands

```bash
rrrm2 list                      # every revision, most current science first
rrrm2 show v13-phospho          # question, conclusion, docs, stage graph, past runs
rrrm2 doctor                    # interpreter, Rscript, data inputs, per-revision runnability
```

```bash
rrrm2 run v13-phospho                              # all required stages
rrrm2 run v13-phospho --dry-run                    # print exact commands, execute nothing
rrrm2 run v13-phospho --with-optional              # include robustness stages
rrrm2 run clinical-axes --stage compartment-context   # that stage plus everything it needs
rrrm2 run clinical-axes --stage figures --only      # that stage alone
rrrm2 run subtype-reference --skip-r                # omit Rscript stages
```

### Where results go

```
data/results/<revision-id>/<UTC-timestamp>[_tag]/
  _logs/<stage>.log      per-stage stdout+stderr
  rrrm2_run.json         command line, exit code, duration, git commit and dirty flag
  …                      whatever the stages themselves write
data/results/<revision-id>/latest -> <most recent run>
```

`--tag mytag` suffixes the run id; `--run-id` sets it outright; `--run-dir`
writes somewhere else entirely.

### Guards

- A stage's declared dependencies are always pulled in and ordered. `--only`
  overrides this.
- `retired` and `blocked` revisions refuse to run without `--allow-retired`.
  Their conclusions are already recorded and their outputs are historical.
- Preflight checks the spec config, entry scripts, declared data inputs, and
  `Rscript` before anything executes. External data bundles under
  `data/external/` are gitignored, so preflight is the fastest way to find out
  what a fresh clone is missing.
- A failing stage stops the run and still writes the manifest, unless
  `--keep-going`.

## Gene panels

Panels live in `config/panels/*.yaml` and are resolved **from disk**, never
from the network, so a run stays deterministic and its provenance hash stays
verifiable.

```bash
rrrm2 panels list                       # every registered panel with its digest
rrrm2 panels show podocyte_core         # members, direction weights, provenance
rrrm2 panels diff podocyte_core barrier_core   # membership and direction differences
rrrm2 panels verify                     # check symbols against Ensembl REST
rrrm2 panels refresh podocyte_core      # verify, then write a dated snapshot
```

`verify` and `refresh` are the only commands that touch the network, and
neither is ever on the analysis path. `refresh` never edits a panel in place —
it writes `config/panels/_snapshots/<panel>.<date>.yaml` carrying the Ensembl
resolution alongside the membership, so promoting it is a deliberate,
reviewable act.

In analysis code, resolve a panel by id rather than pasting symbols:

```python
from rrrm2.panels import get_panel, resolve

symbols = resolve("podocyte_core")               # ['Nphs1', 'Nphs2', ...]
signed  = resolve("barrier_core", signed=True)   # {'Nphs1': -1, ...}
panel   = get_panel("podocyte_core")
panel.digest()                                   # SHA-256 over membership + direction
```

## Adding a revision

1. Write `config/revisions/<id>.yaml`; the `id` must equal the filename stem.
2. State the `question` and, once it is known, the `conclusion`.
3. Set `results_root: data/results/<id>`.
4. List stages in dependency order with `needs:`. Use `{run_dir}`,
   `{spec_config}`, `{figures_root}`, `{run_id}`, `{reference_run}`, `{repo}`
   as placeholders in `args`.
5. `venv/bin/python -m pytest tests/test_rrrm2_registry.py` — the suite checks
   that every entry script exists, every doc exists, every placeholder
   resolves, and the stage graph is acyclic.
