# Node Runner test suite

## Run everything

```powershell
.\venv\Scripts\Activate.ps1
python -m pytest tests -v
```

Add `-x` to stop at the first failure, or `-k testname` to run just one
test. Coverage report (if `pytest-cov` is installed):

```powershell
python -m pytest tests --cov=node_runner --cov-report=term-missing
```

## Layout

| File | What it covers |
|---|---|
| `conftest.py` | Shared fixtures: `tiny_gen`, `two_quads_gen`, `hex_gen` synthetic models |
| `test_engine.py` | BDF round-trip in all three field formats, format detection, free-field converter, expression evaluator (including sandbox-escape attempts), shell + solid quality metrics |
| `test_commands.py` | Every Theme A / Theme C / Convert-Units `Command` class - `execute` and `undo` paths |
| `test_ui.py` | MainWindow boots, every menu entry the v3.0.0 themes added is present, all new dialogs construct, result-browser dock is built |
| `test_b747_generator.py` | `examples/generate_b747.py` produces a strict-parseable BDF with all materials and properties intact |
| `MANUAL_TEST_PLAN.md` | Things the automated suite can't verify - visual rendering, real interactive flows, OP2 post-processing, performance on the 1M model |

## What's NOT covered automatically

- Visual rendering (Femap-style point clouds, Gaussian-mapper-related issues)
- Real interactive workflows: cross-section scrubbing, probe-on-hover, animation playback
- OP2 reader against an actual `.op2` file (we have synthetic data only)
- 1M-node performance and UI responsiveness during scene build
- Cross-platform behavior (the suite assumes Windows + the project venv)

Work through `MANUAL_TEST_PLAN.md` after the automated suite passes to
shake out anything those tests missed.
