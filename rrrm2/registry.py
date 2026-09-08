"""Load and validate the experiment-revision registry in config/revisions/."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterator

import yaml

from .paths import REPO_ROOT, REVISIONS_DIR

VALID_STATUS = frozenset({"locked", "complete", "retired", "exploratory", "blocked"})
VALID_RUNNER = frozenset({"python", "rscript"})


class RegistryError(ValueError):
    """A revision file is malformed or internally inconsistent."""


@dataclass(frozen=True)
class Stage:
    """One executable step of a revision."""

    id: str
    title: str
    runner: str
    entry: str
    args: tuple[str, ...] = ()
    needs: tuple[str, ...] = ()
    optional: bool = False
    outputs: tuple[str, ...] = ()

    @property
    def entry_path(self) -> Path:
        return REPO_ROOT / self.entry

    def resolved_args(self, subs: dict[str, str]) -> list[str]:
        """Substitute ``{run_dir}``-style placeholders into the argument list."""
        out = []
        for arg in self.args:
            try:
                out.append(arg.format(**subs))
            except KeyError as exc:
                raise RegistryError(
                    f"stage {self.id!r} uses unknown placeholder {exc.args[0]!r} "
                    f"in argument {arg!r}"
                ) from exc
        return out


@dataclass(frozen=True)
class Revision:
    """One experiment class, with its spec, docs, results root, and stages."""

    id: str
    title: str
    workflow: str
    status: str
    lock_date: str | None
    question: str
    conclusion: str | None
    spec_config: str | None
    docs: tuple[str, ...]
    results_root: str
    reference_run: str | None
    figures_root: str | None
    requires_data: tuple[str, ...] = ()
    requires_r: bool = False
    stages: tuple[Stage, ...] = field(default=())
    source: Path | None = None

    def stage(self, stage_id: str) -> Stage:
        for st in self.stages:
            if st.id == stage_id:
                return st
        known = ", ".join(s.id for s in self.stages)
        raise RegistryError(
            f"revision {self.id!r} has no stage {stage_id!r}; known stages: {known}"
        )

    def ordered_stages(self, *, with_optional: bool = False) -> list[Stage]:
        """Topologically order stages by ``needs``, dropping optional ones by default."""
        selected = [s for s in self.stages if with_optional or not s.optional]
        keep = {s.id for s in selected}
        done: set[str] = set()
        out: list[Stage] = []
        remaining = list(selected)
        while remaining:
            ready = [s for s in remaining if all(n in done or n not in keep for n in s.needs)]
            if not ready:
                stuck = ", ".join(s.id for s in remaining)
                raise RegistryError(
                    f"revision {self.id!r} has a cyclic or unsatisfiable "
                    f"stage dependency among: {stuck}"
                )
            for s in ready:
                out.append(s)
                done.add(s.id)
                remaining.remove(s)
        return out

    def upstream_of(self, stage_id: str, *, with_optional: bool = False) -> list[Stage]:
        """Return ``stage_id`` preceded by every stage it transitively needs."""
        wanted: set[str] = set()

        def walk(sid: str) -> None:
            if sid in wanted:
                return
            wanted.add(sid)
            for need in self.stage(sid).needs:
                walk(need)

        walk(stage_id)
        return [s for s in self.ordered_stages(with_optional=True)
                if s.id in wanted and (with_optional or not s.optional or s.id == stage_id)]

    @property
    def results_root_path(self) -> Path:
        return REPO_ROOT / self.results_root

    @property
    def reference_run_path(self) -> Path | None:
        return REPO_ROOT / self.reference_run if self.reference_run else None


def _stage_from_dict(rev_id: str, raw: dict[str, Any]) -> Stage:
    for key in ("id", "title", "runner", "entry"):
        if not raw.get(key):
            raise RegistryError(f"{rev_id}: stage is missing required key {key!r}")
    runner = raw["runner"]
    if runner not in VALID_RUNNER:
        raise RegistryError(
            f"{rev_id}/{raw['id']}: runner {runner!r} is not one of "
            f"{sorted(VALID_RUNNER)}"
        )
    return Stage(
        id=raw["id"],
        title=raw["title"],
        runner=runner,
        entry=raw["entry"],
        args=tuple(str(a) for a in (raw.get("args") or ())),
        needs=tuple(raw.get("needs") or ()),
        optional=bool(raw.get("optional", False)),
        outputs=tuple(raw.get("outputs") or ()),
    )


def load_revision(path: Path) -> Revision:
    """Parse and validate a single config/revisions/*.yaml file."""
    raw = yaml.safe_load(path.read_text())
    if not isinstance(raw, dict):
        raise RegistryError(f"{path}: expected a YAML mapping")

    rev_id = raw.get("id")
    if rev_id != path.stem:
        raise RegistryError(
            f"{path}: id {rev_id!r} does not match the filename stem {path.stem!r}"
        )
    status = raw.get("status")
    if status not in VALID_STATUS:
        raise RegistryError(
            f"{path}: status {status!r} is not one of {sorted(VALID_STATUS)}"
        )
    if not raw.get("results_root"):
        raise RegistryError(f"{path}: results_root is required")

    stages = tuple(_stage_from_dict(rev_id, s) for s in (raw.get("stages") or ()))
    ids = [s.id for s in stages]
    if len(ids) != len(set(ids)):
        raise RegistryError(f"{path}: duplicate stage ids")
    for st in stages:
        for need in st.needs:
            if need not in ids:
                raise RegistryError(
                    f"{path}: stage {st.id!r} needs unknown stage {need!r}"
                )

    requires = raw.get("requires") or {}
    rev = Revision(
        id=rev_id,
        title=raw.get("title") or rev_id,
        workflow=str(raw.get("workflow") or "-"),
        status=status,
        lock_date=str(raw["lock_date"]) if raw.get("lock_date") else None,
        question=(raw.get("question") or "").strip(),
        conclusion=(raw.get("conclusion") or "").strip() or None,
        spec_config=raw.get("spec_config"),
        docs=tuple(raw.get("docs") or ()),
        results_root=raw["results_root"],
        reference_run=raw.get("reference_run"),
        figures_root=raw.get("figures_root"),
        requires_data=tuple(requires.get("data") or ()),
        requires_r=bool(requires.get("r", False)),
        stages=stages,
        source=path,
    )
    rev.ordered_stages(with_optional=True)  # raises on a dependency cycle
    return rev


def iter_revision_files(directory: Path | None = None) -> Iterator[Path]:
    """Yield registry files, skipping underscore-prefixed helpers."""
    directory = directory or REVISIONS_DIR
    for path in sorted(directory.glob("*.yaml")):
        if not path.name.startswith("_"):
            yield path


def load_registry(directory: Path | None = None) -> dict[str, Revision]:
    """Load every revision, keyed by id."""
    out: dict[str, Revision] = {}
    for path in iter_revision_files(directory):
        rev = load_revision(path)
        if rev.id in out:
            raise RegistryError(f"duplicate revision id {rev.id!r}")
        out[rev.id] = rev
    return out


def get_revision(rev_id: str, directory: Path | None = None) -> Revision:
    """Look up one revision by id, with a helpful error listing valid ids."""
    registry = load_registry(directory)
    if rev_id not in registry:
        known = ", ".join(sorted(registry))
        raise RegistryError(f"unknown revision {rev_id!r}; known revisions: {known}")
    return registry[rev_id]
