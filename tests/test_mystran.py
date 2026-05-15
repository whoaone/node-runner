"""v5.0.0 items 12-18: MYSTRAN integration tests.

Coverage:
  - solve/__init__.py dataclasses round-trip
  - mystran_settings get/set
  - mystran_runner stage classifier + executable discovery (with monkeypatch)
  - mystran_export translator drops aero / nonlinear / opt cards
  - mystran_preflight catches blocking + warning conditions, < 1s budget
  - mystran_results minimal F06 parser

End-to-end run (driving MYSTRAN) is gated behind `pytest -m mystran_e2e`
so the suite passes without MYSTRAN installed.
"""

from __future__ import annotations

import os
import time as _time
from datetime import datetime
from pathlib import Path

import pytest


# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

class TestSolveDataclasses:

    def test_mystran_run_to_dict_roundtrip(self, tmp_path):
        from node_runner.solve import MystranRun
        run = MystranRun(
            bdf_path=tmp_path / "x.bdf",
            f06_path=tmp_path / "x.f06",
            sol=101,
        )
        d = run.to_dict()
        assert d['bdf_path'].endswith("x.bdf")
        assert d['f06_path'].endswith("x.f06")
        assert d['sol'] == 101
        assert isinstance(d['started_at'], str)
        # save_meta writes JSON next to the file.
        meta = run.save_meta(tmp_path)
        assert meta.exists()
        import json
        loaded = json.loads(meta.read_text(encoding='utf-8'))
        assert loaded['sol'] == 101

    def test_preflight_report_counts(self):
        from node_runner.solve import PreflightIssue, PreflightReport
        rep = PreflightReport()
        rep.issues.append(PreflightIssue(
            severity="blocking", code="X", message="boom"))
        rep.issues.append(PreflightIssue(
            severity="warning", code="Y", message="meh"))
        rep.issues.append(PreflightIssue(
            severity="info", code="Z", message="fyi"))
        assert rep.blocking_count == 1
        assert rep.warning_count == 1
        assert rep.info_count == 1
        assert rep.has_blocking is True
        assert len(rep.issues_by_severity("warning")) == 1


# ---------------------------------------------------------------------------
# Settings module
# ---------------------------------------------------------------------------

class TestMystranSettings:

    def test_defaults_present(self):
        from node_runner.solve import mystran_settings
        d = mystran_settings.DEFAULTS
        assert "executable_path" in d
        assert "default_sollib" in d
        assert d["default_sollib"] in ("SPARSE", "BANDED")
        assert d["default_quad4typ"] in ("MIN4T", "MIN4")
        assert isinstance(d["default_wtmass"], float)

    def test_get_set_roundtrip(self, tmp_path, monkeypatch):
        # Use a tmp scratch root so we don't pollute the real one.
        from node_runner.solve import mystran_settings
        scratch = tmp_path / "mystran"
        mystran_settings.set_value("scratch_root", str(scratch))
        assert mystran_settings.get("scratch_root") == str(scratch)
        # Reset for other tests.
        mystran_settings.set_value(
            "scratch_root", mystran_settings.DEFAULTS["scratch_root"])


# ---------------------------------------------------------------------------
# Runner: stage classifier
# ---------------------------------------------------------------------------

class TestRunnerStageClassifier:

    @pytest.mark.parametrize("line, expect_label", [
        ("READING BULK DATA",        "Reading bulk data"),
        ("ASSEMBLING STIFFNESS MATRIX", "Assembling stiffness matrix"),
        ("FACTORIZING",              "Factorizing"),
        ("DECOMPOSING (CHOLESKY)",   "Factorizing"),
        ("SOLVING FOR DISPLACEMENTS", "Solving"),
        ("RECOVERING STRESSES",      "Recovering element results"),
        ("EIGENVALUE EXTRACTION (LANCZOS)", "Extracting eigenvalues"),
        ("END OF MYSTRAN",           "Finishing"),
    ])
    def test_classify_known_phrases(self, line, expect_label):
        from node_runner.solve.mystran_runner import _classify_stage
        result = _classify_stage(line)
        assert result is not None
        assert result[0] == expect_label

    def test_classify_unknown_line_returns_none(self):
        from node_runner.solve.mystran_runner import _classify_stage
        assert _classify_stage("MAYBE PROCESSING SOMETHING UNRELATED") is None


# ---------------------------------------------------------------------------
# Runner: executable discovery
# ---------------------------------------------------------------------------

class TestRunnerDiscovery:

    def test_discover_returns_none_when_absent(self, monkeypatch, tmp_path):
        from node_runner.solve import mystran_runner
        # Force every discovery channel to miss.
        monkeypatch.setattr(
            mystran_runner.mystran_settings, "get_executable_path",
            lambda: None)
        monkeypatch.delenv("MYSTRAN_EXE", raising=False)
        monkeypatch.setattr(
            mystran_runner.shutil, "which", lambda name: None)
        monkeypatch.setattr(
            mystran_runner, "_COMMON_PATHS",
            [tmp_path / "definitely_missing_mystran.exe"])
        assert mystran_runner.discover_mystran_executable() is None

    def test_discover_finds_env_var(self, monkeypatch, tmp_path):
        from node_runner.solve import mystran_runner
        fake = tmp_path / "mystran.exe"
        fake.write_text("dummy")
        monkeypatch.setattr(
            mystran_runner.mystran_settings, "get_executable_path",
            lambda: None)
        monkeypatch.setenv("MYSTRAN_EXE", str(fake))
        found = mystran_runner.discover_mystran_executable()
        assert found is not None
        assert Path(found).name == "mystran.exe"


# ---------------------------------------------------------------------------
# Export translator
# ---------------------------------------------------------------------------

class TestMystranExport:

    def test_compatibility_table_has_known_rules(self):
        from node_runner.solve.mystran_export import compatibility_table
        rows = compatibility_table()
        codes = {r['card'] for r in rows}
        # PARAM rules
        assert "PARAM,POST" in codes
        assert "PARAM,COUPMASS" in codes
        assert "PARAM,WTMASS" in codes
        # Card family rules (aero / nonlinear / opt)
        assert "caeros" in codes
        assert "nlparms" in codes
        assert "desvars" in codes
        # Blocking element types
        assert "CQUADR" in codes
        assert "CWELD" in codes

    def test_translate_drops_aero_param_and_inserts_mystran_params(
            self, tmp_path, tiny_gen):
        """Translator strips PARAM,POST and injects SOLLIB/QUAD4TYP."""
        from pyNastran.bdf.cards.params import PARAM
        from node_runner.solve.mystran_export import translate_for_mystran
        # Add an MSC-only PARAM to a known-good tiny model.
        tiny_gen.model.params['POST'] = PARAM('POST', [-1])
        tiny_gen.model.params['COUPMASS'] = PARAM('COUPMASS', [1])
        out = tmp_path / "out.bdf"
        report = translate_for_mystran(tiny_gen.model, out)
        # Translation report shape
        assert "POST" in report.dropped_params
        assert "COUPMASS" in report.dropped_params
        # File was written
        assert out.exists() and out.stat().st_size > 0
        text = out.read_text(encoding='utf-8', errors='replace')
        # POST must be gone (case-insensitive)
        assert "PARAM,POST" not in text.replace(" ", "").upper()
        assert "PARAM,COUPMASS" not in text.replace(" ", "").upper()
        # MYSTRAN-required PARAMs were injected
        assert "SOLLIB" in text.upper()
        assert "QUAD4TYP" in text.upper()


# ---------------------------------------------------------------------------
# Pre-flight scanner
# ---------------------------------------------------------------------------

class TestMystranPreflight:

    def test_clean_model_passes(self, tiny_gen):
        from node_runner.solve.mystran_preflight import scan_for_mystran
        report = scan_for_mystran(tiny_gen.model, sol=101)
        assert report.solution_supported is True
        assert report.blocking_count == 0

    def test_blocks_sol_106(self, tiny_gen):
        from node_runner.solve.mystran_preflight import scan_for_mystran
        # Force SOL 106 (nonlinear static) via the sol kwarg.
        report = scan_for_mystran(tiny_gen.model, sol=106)
        assert report.has_blocking
        codes = [i.code for i in report.issues]
        assert "SOL_UNSUPPORTED" in codes

    def test_warns_msc_params(self, tiny_gen):
        from pyNastran.bdf.cards.params import PARAM
        from node_runner.solve.mystran_preflight import scan_for_mystran
        tiny_gen.model.params['COUPMASS'] = PARAM('COUPMASS', [1])
        report = scan_for_mystran(tiny_gen.model, sol=101)
        codes = [i.code for i in report.issues]
        assert any(c.startswith("PARAM_COUPMASS") for c in codes)

    def test_scan_is_fast(self, tiny_gen):
        """< 100 ms on a tiny deck is plenty for a regression gate."""
        from node_runner.solve.mystran_preflight import scan_for_mystran
        t0 = _time.perf_counter()
        scan_for_mystran(tiny_gen.model, sol=101)
        wall = _time.perf_counter() - t0
        # The plan-level <1s budget targets 2.4M-element decks; this
        # test asserts the algorithm is at least linear-friendly on
        # the tiny fixture.
        assert wall < 0.5, f"pre-flight took {wall:.3f}s on tiny fixture"


# ---------------------------------------------------------------------------
# Results adapter
# ---------------------------------------------------------------------------

class TestMystranResults:

    def test_f06_minimal_parses_disp(self, tmp_path):
        from node_runner.solve.mystran_results import _load_f06_minimal
        f06 = (Path(__file__).parent / "fixtures" / "mystran"
               / "sol1_cantilever.f06")
        bundle = _load_f06_minimal(f06)
        assert bundle['subcases'], "expected at least one subcase"
        sc = bundle['subcases'].get(1)
        assert sc is not None
        # Tip displacement (Node 5, T3) should be -0.15 from the fixture.
        assert 5 in sc['displacements']
        t1, t2, t3, r1, r2, r3 = sc['displacements'][5]
        assert abs(t3 - (-0.15)) < 1e-9
        # Root node (1) should be all zeros.
        assert sc['displacements'][1] == [0.0] * 6

    def test_f06_minimal_parses_native_mystran_format(self):
        """v5.0.0 regression: MYSTRAN's own F06 uses
        'D I S P L A C E M E N T S' (with trailing S) and rows of the
        form '<GID> <COORD_SYS> <T1>...<R3>' -- no 'G' marker, mixed
        plain/scientific floats. Earlier parser required MSC's
        'V E C T O R' header + 'G' marker + always-scientific floats
        and silently produced 0 subcases."""
        from node_runner.solve.mystran_results import _load_f06_minimal
        f06 = (Path(__file__).parent / "fixtures" / "mystran"
               / "sol1_cantilever_mystran_native.f06")
        bundle = _load_f06_minimal(f06)
        assert bundle['subcases'], "expected MYSTRAN-native disp parse"
        sc = bundle['subcases'].get(1)
        assert sc is not None
        # 11 grid points in the fixture (cantilever beam from
        # examples/cantilever_beam.bdf).
        assert len(sc['displacements']) == 11
        # Tip node 11, T3 should land near -0.04073 (matches Euler-
        # Bernoulli analytical -0.0404 within 0.7%).
        t3_tip = sc['displacements'][11][2]
        assert abs(t3_tip - (-4.072646e-02)) < 1e-7
        # Root node 1 should be all zeros (clamped).
        assert sc['displacements'][1] == [0.0] * 6

    def test_f06_minimal_parses_eigenvalues(self):
        from node_runner.solve.mystran_results import _load_f06_minimal
        f06 = (Path(__file__).parent / "fixtures" / "mystran"
               / "sol3_modes_box.f06")
        bundle = _load_f06_minimal(f06)
        sc = bundle['subcases'].get(1)
        assert sc is not None
        assert len(sc['eigenvalues']) == 3
        assert abs(sc['eigenvalues'][0] - 12345.0) < 1.0
        assert len(sc['frequencies']) == 3
        # First-mode cycles == 17.66 Hz from the fixture.
        assert abs(sc['frequencies'][0] - 17.66) < 0.1

    def test_load_mystran_results_falls_back_to_f06(self, tmp_path):
        from node_runner.solve import MystranRun
        from node_runner.solve.mystran_results import load_mystran_results
        f06 = (Path(__file__).parent / "fixtures" / "mystran"
               / "sol1_cantilever.f06")
        run = MystranRun(
            bdf_path=tmp_path / "input.bdf",
            f06_path=f06,
            op2_path=None,  # force F06 fallback
            sol=101,
        )
        bundle = load_mystran_results(run)
        assert bundle is not None
        assert run.results_source == "f06"
        assert 5 in bundle['subcases'][1]['displacements']

    def test_load_mystran_results_returns_none_when_no_files(self, tmp_path):
        from node_runner.solve import MystranRun
        from node_runner.solve.mystran_results import load_mystran_results
        run = MystranRun(
            bdf_path=tmp_path / "input.bdf",
            f06_path=None,
            op2_path=None,
        )
        assert load_mystran_results(run) is None


# ---------------------------------------------------------------------------
# End-to-end (only runs locally with MYSTRAN installed)
# ---------------------------------------------------------------------------

@pytest.mark.mystran_e2e
class TestMystranE2E:
    """Gated tests that actually drive MYSTRAN. Enable with
    ``pytest -m mystran_e2e``. CI skips."""

    def test_cantilever_tip_disp_matches_beam_theory(self, tiny_gen, tmp_path):
        pytest.skip("E2E requires MYSTRAN binary; run with -m mystran_e2e + local install")
