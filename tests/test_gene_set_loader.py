"""Enrichr gene-set loader: fetch over HTTPS, fail loudly, never cache an empty library.

Everything here is mocked: no test touches the network, and every cache lives in
``tmp_path`` so the real ``data/processed/gene_sets/KEGG_2019_Mouse.json`` is never
read or written.
"""

from __future__ import annotations

import json
import logging
import sys
import types

import pytest
import requests

from src.enrichment import gene_set_loader as gsl

LIBRARY = "KEGG_2019_Mouse"

# Enrichr mode=text layout: term, empty description, genes; every line ends with a
# trailing tab; weighted libraries write GENE,weight; blank lines occur.
ENRICHR_TEXT = (
    "Focal adhesion\t\tActn1\tItga1\tVcl\t\n"
    "\n"
    "Tight junction\t\tCldn1\tTjp1\t\n"
    "Weighted set\t\tGeneA,1.5\tGeneB,0.2\t\n"
    "No genes at all\t\t\n"
)
PARSED = {
    "Focal adhesion": ["Actn1", "Itga1", "Vcl"],
    "Tight junction": ["Cldn1", "Tjp1"],
    "Weighted set": ["GeneA", "GeneB"],
}


class FakeResponse:
    def __init__(self, body: bytes = b"", status: int = 200, payload=None):
        self.content = body
        self.status_code = status
        self._payload = payload

    def raise_for_status(self) -> None:
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code} Client Error")

    def json(self):
        return self._payload


@pytest.fixture
def fake_requests(monkeypatch):
    """Replace ``requests.get`` and record the calls; behaviour is set per test."""
    calls: list[dict] = []
    state = {"handler": lambda url, **kw: FakeResponse(ENRICHR_TEXT.encode())}

    def fake_get(url, **kwargs):
        calls.append({"url": url, **kwargs})
        if not url.startswith("https://"):
            # what the HTTPS-only egress proxy does to plain HTTP
            raise requests.ConnectionError(f"proxy refused plain HTTP: {url}")
        return state["handler"](url, **kwargs)

    monkeypatch.setattr(requests, "get", fake_get)
    return types.SimpleNamespace(calls=calls, state=state)


def install_fake_gseapy(monkeypatch, get_library=None, get_library_name=None):
    module = types.ModuleType("gseapy")
    if get_library is not None:
        module.get_library = get_library
    if get_library_name is not None:
        module.get_library_name = get_library_name
    monkeypatch.setitem(sys.modules, "gseapy", module)
    return module


def block_gseapy(monkeypatch):
    """Importing gseapy raises ImportError (the package is not installed)."""
    monkeypatch.setitem(sys.modules, "gseapy", None)


def no_json_files(directory) -> bool:
    return not list(directory.glob("*.json"))


# --- parsing ------------------------------------------------------------------------

def test_parse_enrichr_text_matches_the_gseapy_term_to_genes_layout():
    assert gsl.parse_enrichr_text(ENRICHR_TEXT) == PARSED


def test_parse_enrichr_text_keeps_gene_order_and_duplicates():
    text = "Set\t\tB\tA\tB\tA\n"
    assert gsl.parse_enrichr_text(text) == {"Set": ["B", "A", "B", "A"]}


def test_parse_enrichr_text_of_an_empty_body_is_empty():
    assert gsl.parse_enrichr_text("") == {}
    assert gsl.parse_enrichr_text("\n\n  \n") == {}


# --- HTTPS fetch ----------------------------------------------------------------------

def test_fetch_uses_https_enrichr_text_endpoint_and_caches_json(tmp_path, fake_requests, monkeypatch):
    install_fake_gseapy(monkeypatch, get_library=lambda name: pytest.fail("gseapy must not be used"))

    library = gsl.fetch_enrichr_library(LIBRARY, tmp_path)

    assert library == PARSED
    (call,) = fake_requests.calls
    assert call["url"] == "https://maayanlab.cloud/Enrichr/geneSetLibrary"
    assert call["params"] == {"mode": "text", "libraryName": LIBRARY}
    assert call["timeout"] == gsl.ENRICHR_TIMEOUT_S
    assert json.loads((tmp_path / f"{LIBRARY}.json").read_text()) == PARSED
    assert sorted(p.name for p in tmp_path.iterdir()) == [f"{LIBRARY}.json"], "no temp files left"


def test_fetch_creates_a_missing_cache_directory(tmp_path, fake_requests):
    cache_dir = tmp_path / "does" / "not" / "exist"
    assert gsl.fetch_enrichr_library(LIBRARY, cache_dir) == PARSED
    assert (cache_dir / f"{LIBRARY}.json").exists()


def test_fetch_works_when_plain_http_is_blocked(tmp_path, fake_requests, monkeypatch):
    """gseapy 1.1.1 only speaks http://; behind an HTTPS-only proxy that fails."""
    def gseapy_http_only(name):
        raise requests.ConnectionError("proxy refused plain HTTP: http://maayanlab.cloud")

    install_fake_gseapy(monkeypatch, get_library=gseapy_http_only)
    assert gsl.fetch_enrichr_library(LIBRARY, tmp_path) == PARSED
    assert all(c["url"].startswith("https://") for c in fake_requests.calls)


def test_utf8_gene_names_survive_the_round_trip(tmp_path, fake_requests):
    fake_requests.state["handler"] = lambda url, **kw: FakeResponse("Set\t\tGène1\tGene2\n".encode())
    assert gsl.fetch_enrichr_library(LIBRARY, tmp_path) == {"Set": ["Gène1", "Gene2"]}


# --- cache ------------------------------------------------------------------------------

def test_a_populated_cache_is_used_without_any_network_call(tmp_path, fake_requests, monkeypatch):
    (tmp_path / f"{LIBRARY}.json").write_text(json.dumps(PARSED))
    install_fake_gseapy(monkeypatch, get_library=lambda name: pytest.fail("no gseapy"))

    assert gsl.fetch_enrichr_library(LIBRARY, tmp_path) == PARSED
    assert fake_requests.calls == []


@pytest.mark.parametrize("poison", ["{}", '{"Empty set": []}', "not json at all", "[]"])
def test_an_empty_or_corrupt_cache_is_ignored_and_replaced(tmp_path, fake_requests, caplog, poison):
    cache_file = tmp_path / f"{LIBRARY}.json"
    cache_file.write_text(poison)

    with caplog.at_level(logging.ERROR, logger=gsl.__name__):
        library = gsl.fetch_enrichr_library(LIBRARY, tmp_path)

    assert library == PARSED
    assert json.loads(cache_file.read_text()) == PARSED
    assert any("Ignoring" in r.message and r.levelno == logging.ERROR for r in caplog.records)


# --- fallback and failure -------------------------------------------------------------------

def test_gseapy_is_the_fallback_when_https_fails(tmp_path, fake_requests, monkeypatch):
    fake_requests.state["handler"] = lambda url, **kw: FakeResponse(status=503)
    install_fake_gseapy(monkeypatch, get_library=lambda name: {"From gseapy": ["G1", "G2"]})

    library = gsl.fetch_enrichr_library(LIBRARY, tmp_path)

    assert library == {"From gseapy": ["G1", "G2"]}
    assert json.loads((tmp_path / f"{LIBRARY}.json").read_text()) == library


def test_failure_raises_loudly_and_caches_nothing(tmp_path, fake_requests, monkeypatch, caplog):
    fake_requests.state["handler"] = lambda url, **kw: FakeResponse(status=403)
    install_fake_gseapy(monkeypatch, get_library=lambda name: (_ for _ in ()).throw(
        requests.ConnectionError("proxy refused plain HTTP")))

    with caplog.at_level(logging.ERROR, logger=gsl.__name__):
        with pytest.raises(gsl.GeneSetFetchError) as err:
            gsl.fetch_enrichr_library(LIBRARY, tmp_path)

    message = str(err.value)
    assert LIBRARY in message and "refusing to continue with an empty library" in message
    assert "https:" in message and "403" in message
    assert "gseapy:" in message and "proxy refused plain HTTP" in message
    assert no_json_files(tmp_path), "a failed fetch must not leave a cache file"
    assert any(r.levelno == logging.ERROR and LIBRARY in r.getMessage() for r in caplog.records)


def test_fetch_error_is_a_runtime_error():
    assert issubclass(gsl.GeneSetFetchError, RuntimeError)


def test_an_empty_https_body_is_a_failure_not_an_empty_library(tmp_path, fake_requests, monkeypatch):
    fake_requests.state["handler"] = lambda url, **kw: FakeResponse(b"\n\n")
    install_fake_gseapy(monkeypatch, get_library=lambda name: {})

    with pytest.raises(gsl.GeneSetFetchError, match="returned an empty library"):
        gsl.fetch_enrichr_library(LIBRARY, tmp_path)
    assert no_json_files(tmp_path)


def test_a_missing_gseapy_package_is_reported_with_an_install_hint(tmp_path, fake_requests, monkeypatch):
    fake_requests.state["handler"] = lambda url, **kw: (_ for _ in ()).throw(
        requests.ConnectionError("network unreachable"))
    block_gseapy(monkeypatch)

    with pytest.raises(gsl.GeneSetFetchError, match="pip install gseapy"):
        gsl.fetch_enrichr_library(LIBRARY, tmp_path)
    assert no_json_files(tmp_path)


def test_non_strict_mode_logs_an_error_and_still_refuses_to_cache(
        tmp_path, fake_requests, monkeypatch, caplog):
    fake_requests.state["handler"] = lambda url, **kw: FakeResponse(status=500)
    block_gseapy(monkeypatch)

    with caplog.at_level(logging.ERROR, logger=gsl.__name__):
        assert gsl.fetch_enrichr_library(LIBRARY, tmp_path, strict=False) == {}

    assert any(r.levelno == logging.ERROR for r in caplog.records)
    assert no_json_files(tmp_path)


def test_a_failed_atomic_write_leaves_no_partial_cache(tmp_path, monkeypatch):
    def broken_dump(obj, handle, **kwargs):
        handle.write("{")
        raise OSError("disk full")

    monkeypatch.setattr(gsl.json, "dump", broken_dump)
    with pytest.raises(OSError, match="disk full"):
        gsl._write_library_cache(tmp_path / f"{LIBRARY}.json", PARSED)
    assert list(tmp_path.iterdir()) == []


# --- load_gene_sets ----------------------------------------------------------------------------

def _load(tmp_path, libraries, **kwargs):
    return gsl.load_gene_sets(
        libraries=libraries,
        cache_dir=tmp_path,
        id_map_path=tmp_path / "no_id_map.json",
        include_curated=False,
        include_segment_markers=False,
        min_size=2,
        **kwargs,
    )


def test_load_gene_sets_raises_when_a_requested_library_is_unavailable(
        tmp_path, fake_requests, monkeypatch):
    fake_requests.state["handler"] = lambda url, **kw: FakeResponse(status=403)
    block_gseapy(monkeypatch)

    with pytest.raises(gsl.GeneSetFetchError):
        _load(tmp_path, [LIBRARY])


def test_load_gene_sets_non_strict_skips_the_library_and_names_it_in_the_error(
        tmp_path, fake_requests, monkeypatch, caplog):
    (tmp_path / "Cached_Lib.json").write_text(json.dumps(PARSED))
    fake_requests.state["handler"] = lambda url, **kw: FakeResponse(status=403)
    block_gseapy(monkeypatch)

    with caplog.at_level(logging.ERROR, logger=gsl.__name__):
        sets, set_to_library = _load(tmp_path, ["Cached_Lib", LIBRARY], strict=False)

    assert set(set_to_library.values()) == {"Cached_Lib"}
    assert "Cached_Lib::Focal adhesion" in sets
    assert any("NOT loaded" in r.getMessage() and LIBRARY in r.getMessage()
               for r in caplog.records)


def test_load_gene_sets_reads_a_fetched_library_and_filters_by_size(tmp_path, fake_requests):
    sets, set_to_library = _load(tmp_path, [LIBRARY], max_size=2)

    assert sets == {f"{LIBRARY}::Tight junction": ["Cldn1", "Tjp1"],
                    f"{LIBRARY}::Weighted set": ["GeneA", "GeneB"]}
    assert set(set_to_library.values()) == {LIBRARY}
    assert len(fake_requests.calls) == 1


# --- listing -----------------------------------------------------------------------------------

def test_list_enrichr_libraries_reads_the_https_statistics_endpoint(fake_requests):
    payload = {"statistics": [{"libraryName": "Reactome_2022"}, {"libraryName": "KEGG_2019_Mouse"}]}
    fake_requests.state["handler"] = lambda url, **kw: FakeResponse(payload=payload)

    assert gsl.list_enrichr_libraries() == ["KEGG_2019_Mouse", "Reactome_2022"]
    assert fake_requests.calls[0]["url"] == "https://maayanlab.cloud/Enrichr/datasetStatistics"


def test_list_enrichr_libraries_fails_loudly_when_nothing_answers(fake_requests, monkeypatch):
    fake_requests.state["handler"] = lambda url, **kw: FakeResponse(status=502)
    block_gseapy(monkeypatch)
    with pytest.raises(gsl.GeneSetFetchError, match="Could not list Enrichr libraries"):
        gsl.list_enrichr_libraries()
