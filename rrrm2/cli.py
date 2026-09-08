"""Command-line entry point: list, show, run, panels, doctor."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
import sys
import textwrap

from . import panels as panels_mod
from . import runner
from .paths import REPO_ROOT, RESULTS_DIR, relative
from .registry import RegistryError, Revision, load_registry, get_revision

STATUS_ORDER = {"locked": 0, "complete": 1, "blocked": 2, "exploratory": 3, "retired": 4}


def _wrap(text: str, indent: str = "  ", width: int = 88) -> str:
    return textwrap.fill(" ".join(text.split()), width=width,
                         initial_indent=indent, subsequent_indent=indent)


# --- list -------------------------------------------------------------------

def cmd_list(args: argparse.Namespace) -> int:
    """Print every registered revision, newest science first."""
    registry = load_registry()
    revs = sorted(registry.values(),
                  key=lambda r: (STATUS_ORDER.get(r.status, 9), r.id))
    if args.json:
        print(json.dumps([{
            "id": r.id, "title": r.title, "workflow": r.workflow,
            "status": r.status, "lock_date": r.lock_date,
            "stages": [s.id for s in r.stages],
            "results_root": r.results_root,
        } for r in revs], indent=2))
        return 0

    print(f"{len(revs)} registered revisions in {relative(REPO_ROOT / 'config/revisions')}\n")
    print(f"  {'ID':<20} {'WF':<3} {'STATUS':<11} {'STAGES':>6}  TITLE")
    print(f"  {'-'*20} {'-'*3} {'-'*11} {'-'*6}  {'-'*40}")
    for r in revs:
        n = f"{len(r.ordered_stages())}/{len(r.stages)}"
        print(f"  {r.id:<20} {r.workflow:<3} {r.status:<11} {n:>6}  {r.title}")
    print("\n  STAGES is required/total; optional stages need --with-optional.")
    print("  `rrrm2 show <id>` for the question, conclusion, docs, and stage graph.")
    return 0


# --- show -------------------------------------------------------------------

def _show_revision(rev: Revision, as_json: bool) -> int:
    if as_json:
        print(json.dumps({
            "id": rev.id, "title": rev.title, "workflow": rev.workflow,
            "status": rev.status, "lock_date": rev.lock_date,
            "question": rev.question, "conclusion": rev.conclusion,
            "spec_config": rev.spec_config, "docs": list(rev.docs),
            "results_root": rev.results_root, "reference_run": rev.reference_run,
            "stages": [{"id": s.id, "title": s.title, "runner": s.runner,
                        "entry": s.entry, "needs": list(s.needs),
                        "optional": s.optional, "outputs": list(s.outputs)}
                       for s in rev.stages],
        }, indent=2))
        return 0

    print(f"\n{rev.id} — {rev.title}")
    print("=" * min(88, len(rev.id) + len(rev.title) + 3))
    lock = f", locked {rev.lock_date}" if rev.lock_date else ""
    print(f"workflow {rev.workflow}   status {rev.status}{lock}\n")

    print("QUESTION")
    print(_wrap(rev.question) + "\n")
    if rev.conclusion:
        print("CONCLUSION")
        print(_wrap(rev.conclusion) + "\n")

    print("SPEC")
    print(f"  config    {rev.spec_config or '(none)'}")
    print(f"  results   {rev.results_root}")
    print(f"  reference {rev.reference_run or '(none)'}")
    if rev.figures_root:
        print(f"  figures   {rev.figures_root}")
    print()

    if rev.docs:
        print("DOCS")
        for d in rev.docs:
            mark = " " if (REPO_ROOT / d).exists() else "?"
            print(f" {mark} {d}")
        print()

    print("STAGES")
    for st in rev.ordered_stages(with_optional=True):
        flag = "(optional)" if st.optional else ""
        needs = f"  needs: {', '.join(st.needs)}" if st.needs else ""
        print(f"  {st.id:<22} {st.runner:<8} {st.title} {flag}")
        print(f"  {'':22} {st.entry}{needs}")
    print()

    root = rev.results_root_path
    if root.exists():
        runs = sorted(p.name for p in root.iterdir() if p.is_dir())
        print(f"RUNS ({len(runs)}) under {rev.results_root}")
        for name in runs[-5:]:
            print(f"  {name}")
    else:
        print(f"RUNS  none yet under {rev.results_root}")
    print()
    return 0


def cmd_show(args: argparse.Namespace) -> int:
    """Print one revision in full."""
    return _show_revision(get_revision(args.revision), args.json)


# --- run --------------------------------------------------------------------

def cmd_run(args: argparse.Namespace) -> int:
    """Execute a revision's stages into a fresh run directory."""
    rev = get_revision(args.revision)

    if args.stage:
        stages = rev.upstream_of(args.stage, with_optional=args.with_optional)
        if args.only:
            stages = [rev.stage(args.stage)]
    else:
        stages = rev.ordered_stages(with_optional=args.with_optional)
    if args.skip_r:
        stages = [s for s in stages if s.runner != "rscript"]
    if not stages:
        print("nothing to run after filtering", file=sys.stderr)
        return 1

    check = runner.preflight(rev, stages=stages)
    for w in check.warnings:
        print(f"warning: {w}", file=sys.stderr)
    if not check.ok:
        print("\npreflight failed:", file=sys.stderr)
        for p in check.problems:
            print(f"  - {p}", file=sys.stderr)
        if not args.dry_run:
            print("\nre-run with --dry-run to print commands without executing.",
                  file=sys.stderr)
            return 2

    if rev.status in {"retired", "blocked"} and not args.allow_retired and not args.dry_run:
        print(
            f"\nrefusing to run {rev.id!r}: status is {rev.status!r}.\n"
            "Its conclusion is already recorded and its outputs are historical.\n"
            "Pass --allow-retired to reproduce it anyway.", file=sys.stderr)
        return 3

    run_id = args.run_id or runner.new_run_id(args.tag)
    if args.run_dir:
        run_dir = Path(args.run_dir).resolve()
        run_dir.mkdir(parents=True, exist_ok=True)
    else:
        run_dir = (rev.results_root_path / run_id if args.dry_run
                   else runner.allocate_run_dir(rev, run_id))

    subs = runner.substitutions(rev, run_dir, run_id)
    python = args.python or runner.default_python()

    print(f"\nrevision  {rev.id}  ({rev.status})")
    print(f"run id    {run_id}")
    print(f"run dir   {relative(run_dir)}")
    print(f"stages    {', '.join(s.id for s in stages)}\n")

    results: list[dict] = []
    for i, st in enumerate(stages, 1):
        print(f"[{i}/{len(stages)}] {st.id} — {st.title}")
        res = runner.run_stage(st, subs, python=python, run_dir=run_dir,
                               dry_run=args.dry_run, extra_args=args.extra or None)
        results.append(res)
        if res["status"] == "failed":
            print(f"\nstage {st.id!r} failed with exit code {res['returncode']}; "
                  f"log: {res['log']}", file=sys.stderr)
            if not args.keep_going:
                runner.write_run_manifest(run_dir, rev, run_id, results)
                return res["returncode"] or 1
        print()

    if not args.dry_run:
        manifest = runner.write_run_manifest(run_dir, rev, run_id, results)
        print(f"manifest  {relative(manifest)}")
        failed = [r for r in results if r["status"] == "failed"]
        print(f"done      {len(results) - len(failed)}/{len(results)} stages ok")
        return 1 if failed else 0
    return 0


# --- panels -----------------------------------------------------------------

def cmd_panels(args: argparse.Namespace) -> int:
    """List, show, verify, refresh, or diff registered gene panels."""
    registry = panels_mod.load_panels()

    if args.panels_command == "list":
        if not registry:
            print("no panels registered yet under config/panels/")
            return 0
        if args.json:
            print(json.dumps({p.id: {"species": p.species, "n_genes": len(p.genes),
                                     "frozen": p.frozen, "source": relative(p.source),
                                     "digest": p.digest()}
                              for p in registry.values()}, indent=2))
            return 0
        print(f"{len(registry)} panels\n")
        print(f"  {'PANEL':<34} {'SPECIES':<7} {'N':>4} {'FROZEN':<7} SOURCE")
        for p in sorted(registry.values(), key=lambda x: x.id):
            print(f"  {p.id:<34} {p.species:<7} {len(p.genes):>4} "
                  f"{'yes' if p.frozen else '':<7} {relative(p.source)}")
        return 0

    if args.panels_command == "show":
        p = panels_mod.get_panel(args.panel)
        print(f"\n{p.id}")
        print(f"  {p.description}")
        print(f"  species {p.species}   genes {len(p.genes)}   "
              f"frozen {'yes' if p.frozen else 'no'}")
        print(f"  source  {relative(p.source)}")
        print(f"  digest  {p.digest()}")
        if p.minimum_present is not None:
            print(f"  minimum_present {p.minimum_present}")
        print("\n  genes")
        for g, d in p.genes.items():
            print(f"    {'+' if d > 0 else '-'} {g}")
        prov = p.provenance.get("extracted_from") or []
        if prov:
            print("\n  extracted from")
            for e in prov:
                print(f"    {e.get('file')}  {e.get('symbol','')} {e.get('key','')}")
        print()
        return 0

    if args.panels_command == "diff":
        result = panels_mod.diff_panels(panels_mod.get_panel(args.left),
                                        panels_mod.get_panel(args.right))
        print(json.dumps(result, indent=2))
        return 0 if result["identical"] else 1

    if args.panels_command in {"verify", "refresh"}:
        targets = ([panels_mod.get_panel(args.panel)] if args.panel
                   else sorted(registry.values(), key=lambda x: x.id))
        if not targets:
            print("no panels registered yet under config/panels/", file=sys.stderr)
            return 1
        failures = 0
        for p in targets:
            try:
                result = panels_mod.verify_panel(p, timeout=args.timeout)
            except Exception as exc:  # network, DNS, HTTP, or payload shape
                print(f"  {p.id:<34} ERROR  {type(exc).__name__}: {exc}",
                      file=sys.stderr)
                failures += 1
                continue
            bad = result["unresolved"]
            ren = result["renamed"]
            state = "ok" if not bad and not ren else "CHECK"
            print(f"  {p.id:<34} {state:<6} "
                  f"{result['n_resolved']}/{result['n_genes']} resolved"
                  + (f", unresolved: {', '.join(bad)}" if bad else "")
                  + (f", renamed: {', '.join(r['stored']+'->'+r['current'] for r in ren)}"
                     if ren else ""))
            if bad or ren:
                failures += 1
            if args.panels_command == "refresh":
                out = panels_mod.write_snapshot(p, result)
                print(f"  {'':<34} snapshot {relative(out)}")
        if args.panels_command == "verify" and failures:
            print(f"\n{failures} panel(s) need attention.", file=sys.stderr)
        return 1 if failures else 0

    return 1


# --- doctor -----------------------------------------------------------------

def cmd_doctor(args: argparse.Namespace) -> int:
    """Report interpreter, data, and per-revision runnability."""
    print(f"repo root   {REPO_ROOT}")
    python = runner.default_python()
    print(f"python      {python}")
    print(f"Rscript     {shutil.which('Rscript') or '(not on PATH)'}")
    print(f"results     {relative(RESULTS_DIR)}"
          f"{'' if RESULTS_DIR.exists() else '  (missing)'}")

    n_panels = len(panels_mod.load_panels())
    print(f"panels      {n_panels} registered")

    print("\nrevision runnability")
    registry = load_registry()
    blocked = 0
    for rev in sorted(registry.values(), key=lambda r: r.id):
        check = runner.preflight(rev)
        if check.ok:
            print(f"  {rev.id:<20} ok")
        else:
            blocked += 1
            print(f"  {rev.id:<20} blocked")
            for p in check.problems:
                print(f"  {'':20}   - {p}")
    print(f"\n{len(registry) - blocked}/{len(registry)} revisions runnable here.")
    return 0


# --- parser -----------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    """Assemble the rrrm2 argument parser."""
    p = argparse.ArgumentParser(
        prog="rrrm2",
        description="Run any concluded revision of the RRRM-2 kidney reanalysis.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            examples:
              rrrm2 list
              rrrm2 show v13-phospho
              rrrm2 run v13-phospho --dry-run
              rrrm2 run clinical-axes --stage compartment-context
              rrrm2 panels verify
              rrrm2 doctor
        """))
    sub = p.add_subparsers(dest="command", required=True)

    pl = sub.add_parser("list", help="list every registered revision")
    pl.add_argument("--json", action="store_true")
    pl.set_defaults(func=cmd_list)

    ps = sub.add_parser("show", help="show one revision in full")
    ps.add_argument("revision")
    ps.add_argument("--json", action="store_true")
    ps.set_defaults(func=cmd_show)

    pr = sub.add_parser("run", help="run a revision into its results folder")
    pr.add_argument("revision")
    pr.add_argument("--stage", help="run this stage and everything it needs")
    pr.add_argument("--only", action="store_true",
                    help="with --stage, run only that stage")
    pr.add_argument("--with-optional", action="store_true",
                    help="include stages marked optional")
    pr.add_argument("--skip-r", action="store_true", help="omit Rscript stages")
    pr.add_argument("--dry-run", action="store_true",
                    help="print the exact commands without executing")
    pr.add_argument("--run-id", help="use this run id instead of a UTC timestamp")
    pr.add_argument("--tag", help="suffix appended to the generated run id")
    pr.add_argument("--run-dir", help="write to this directory instead of results_root")
    pr.add_argument("--python", help="interpreter for python stages")
    pr.add_argument("--keep-going", action="store_true",
                    help="continue after a failing stage")
    pr.add_argument("--allow-retired", action="store_true",
                    help="permit running a retired or blocked revision")
    pr.add_argument("--extra", nargs=argparse.REMAINDER,
                    help="extra arguments appended to every stage command")
    pr.set_defaults(func=cmd_run)

    pp = sub.add_parser("panels", help="inspect and verify the gene-panel registry")
    psub = pp.add_subparsers(dest="panels_command", required=True)
    pp_list = psub.add_parser("list", help="list registered panels")
    pp_list.add_argument("--json", action="store_true")
    pp_show = psub.add_parser("show", help="show one panel's members and provenance")
    pp_show.add_argument("panel")
    pp_ver = psub.add_parser("verify", help="check symbols against Ensembl REST")
    pp_ver.add_argument("panel", nargs="?")
    pp_ver.add_argument("--timeout", type=int, default=60)
    pp_ref = psub.add_parser("refresh",
                             help="verify and write a dated snapshot (never in place)")
    pp_ref.add_argument("panel", nargs="?")
    pp_ref.add_argument("--timeout", type=int, default=60)
    pp_diff = psub.add_parser("diff", help="compare two panels' membership")
    pp_diff.add_argument("left")
    pp_diff.add_argument("right")
    pp.set_defaults(func=cmd_panels)

    pd = sub.add_parser("doctor", help="check the environment and per-revision inputs")
    pd.set_defaults(func=cmd_doctor)

    return p


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and dispatch."""
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (RegistryError, panels_mod.PanelError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("\ninterrupted", file=sys.stderr)
        return 130
