"""PyInstaller runtime hook: register bundled .libs folders as DLL search paths.

The vtk wheel ships its native .dlls in `vtk.libs/` (delvewheel layout) and
patches `vtkmodules/__init__.py` with a `_delvewheel_patch` that calls
`os.add_dll_directory(libs_dir)` where libs_dir is computed relative to
`__file__`. Inside a PyInstaller --onefile bundle that relative path can
fall through depending on extraction layout, leaving the .pyd's static
DLL imports (vtkCommonCore-9.5.0.dll etc.) unresolved. The exe then
fails on import with "DLL load failed while importing vtkCommonCore".

Belt-and-suspenders: register every bundled `*.libs` folder with both
`os.add_dll_directory` and PATH at startup so the loader can always find
them, regardless of whether the in-package patch ran. Also covers
numpy.libs and scipy.libs.
"""

import os
import sys


def _register_libs_dirs():
    base = getattr(sys, "_MEIPASS", None)
    if not base:
        return
    add = getattr(os, "add_dll_directory", None)
    if add is None:
        return
    for entry in os.listdir(base):
        full = os.path.join(base, entry)
        if entry.endswith(".libs") and os.path.isdir(full):
            try:
                add(full)
                os.environ["PATH"] = full + os.pathsep + os.environ.get("PATH", "")
            except OSError:
                pass


_register_libs_dirs()
