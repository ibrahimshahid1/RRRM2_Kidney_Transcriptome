"""Preflight checks, run-directory allocation, and stage execution."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from .paths import REPO_ROOT, relative
from .registry import Revision, Stage


@dataclass
class Preflight:
    """Result of checking that a revision can run at all."""

    ok: bool
    problems: list[str]
    warnings: list[str]


def preflight(revision: Revision, *, stages: list[Stage] | None = None) -> Preflight:
    """Check interpreters, entry scripts, spec config, and declared data inputs."""
    problems: list[str] = []
    warnings: list[str] = []
    stages = stages if stages is not None else list(revision.stages)

    if revision.spec_config and not (REPO_ROOT / revision.spec_config).exists():
        problems.append(f"missing spec config: {revision.spec_config}")

    for entry in revision.requires_data:
        if not (REPO_ROOT / entry).exists():
            problems.append(
                f"missing declared data input: {entry} "
                "(external bundles are gitignored; see README 'Reproduction')"
            )

    for st in stages:
        if not st.entry_path.exists():
            problems.append(f"stage {st.id}: entry not found: {st.entry}")

    if any(st.runner == "rscript" for st in stages) and shutil.which("Rscript") is None:
        problems.append(
            "Rscript is not on PATH but this revision has R stages "
            "(re-run with --skip-r to omit them)"
        )

    if revision.status in {"retired", "blocked"}:
        warnings.append(
            f"revision status is {revision.status!r} — "
            "its outputs are historical and must not be cited as current results"
        )
    return Preflight(ok=not problems, problems=problems, warnings=warnings)


def new_run_id(tag: str | None = None) -> str:
    """Build a sortable run id, optionally suffixed with a user tag."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"{stamp}_{tag}" if tag else stamp


def allocate_run_dir(revision: Revision, run_id: str) -> Path:
    """Create ``data/results/<revision>/<run_id>/`` and refresh the `latest` link."""
    run_dir = revision.results_root_path / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    link = revision.results_root_path / "latest"
    try:
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(run_dir.name)
    except OSError:
        pass  # a filesystem without symlinks is not a reason to fail the run
    return run_dir


def substitutions(revision: Revision, run_dir: Path, run_id: str) -> dict[str, str]:
    """Placeholder values available to a stage's argument list."""
    return {
        "repo": str(REPO_ROOT),
        "run_dir": str(run_dir),
        "run_id": run_id,
        "spec_config": str(REPO_ROOT / revision.spec_config) if revision.spec_config else "",
        "results_root": str(revision.results_root_path),
        "reference_run": str(revision.reference_run_path or ""),
        "figures_root": str(REPO_ROOT / revision.figures_root) if revision.figures_root else "",
    }


def build_command(stage: Stage, subs: dict[str, str], python: str) -> list[str]:
    """Assemble the argv for one stage."""
    exe = [python] if stage.runner == "python" else ["Rscript"]
    return [*exe, str(stage.entry_path), *stage.resolved_args(subs)]


def default_python() -> str:
    """Prefer the repository venv interpreter, else the current one."""
    venv = REPO_ROOT / "venv" / "bin" / "python"
    return str(venv) if venv.exists() else sys.executable


def run_stage(
    stage: Stage,
    subs: dict[str, str],
    *,
    python: str,
    run_dir: Path,
    dry_run: bool = False,
    extra_args: list[str] | None = None,
) -> dict:
    """Execute one stage, streaming output and tee-ing it into the run directory."""
    cmd = build_command(stage, subs, python)
    if extra_args:
        cmd.extend(extra_args)
    printable = " ".join(cmd)

    if dry_run:
        print(f"  [dry-run] {printable}")
        return {"stage": stage.id, "command": cmd, "status": "dry-run"}

    env = dict(os.environ)
    env.setdefault("PYTHONPATH", str(REPO_ROOT))
    env.setdefault("MPLCONFIGDIR", str(run_dir / ".mplcache"))
    (run_dir / ".mplcache").mkdir(parents=True, exist_ok=True)

    log_dir = run_dir / "_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{stage.id}.log"

    print(f"  $ {printable}")
    started = time.time()
    with log_path.open("w") as log:
        log.write(f"# {printable}\n# started {datetime.now(timezone.utc).isoformat()}\n\n")
        log.flush()
        proc = subprocess.Popen(
            cmd, cwd=REPO_ROOT, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            log.write(line)
        code = proc.wait()
    elapsed = time.time() - started

    return {
        "stage": stage.id,
        "command": cmd,
        "status": "ok" if code == 0 else "failed",
        "returncode": code,
        "seconds": round(elapsed, 2),
        "log": relative(log_path),
    }


def write_run_manifest(run_dir: Path, revision: Revision, run_id: str,
                       results: list[dict]) -> Path:
    """Record what the CLI ran, so a run directory explains itself."""
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT,
            capture_output=True, text=True, check=False).stdout.strip() or None
        dirty = bool(subprocess.run(
            ["git", "status", "--porcelain"], cwd=REPO_ROOT,
            capture_output=True, text=True, check=False).stdout.strip())
    except OSError:
        commit, dirty = None, None

    manifest = {
        "revision": revision.id,
        "title": revision.title,
        "status": revision.status,
        "run_id": run_id,
        "run_dir": relative(run_dir),
        "spec_config": revision.spec_config,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": commit,
        "git_dirty": dirty,
        "python": sys.version.split()[0],
        "stages": results,
    }
    path = run_dir / "rrrm2_run.json"
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path
