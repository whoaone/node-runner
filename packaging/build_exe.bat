@echo off
REM Build Node Runner.exe using PyInstaller from the bundled venv.
REM Uses the runtime hook in packaging/rthook_vtk_libs.py to register
REM vtk.libs / numpy.libs / scipy.libs as DLL search paths so VTK loads
REM correctly inside the --onefile bundle.

setlocal
set "HERE=%~dp0"
set "ROOT=%HERE%.."
set "PY=%ROOT%\venv\Scripts\python.exe"

if not exist "%PY%" (
  echo Could not find venv Python at "%PY%". Create the venv first.
  exit /b 1
)

cd /d "%ROOT%"
"%PY%" -m PyInstaller --noconfirm --onefile --windowed ^
  --name "Node Runner" ^
  --add-data "images;images" ^
  --add-binary "venv/Lib/site-packages/vtk.libs/*;vtk.libs" ^
  --add-binary "venv/Lib/site-packages/numpy.libs/*;numpy.libs" ^
  --add-binary "venv/Lib/site-packages/scipy.libs/*;scipy.libs" ^
  --collect-all vtkmodules ^
  --collect-all pyvista ^
  --collect-all pyvistaqt ^
  --collect-all pyNastran ^
  --runtime-hook packaging/rthook_vtk_libs.py ^
  run.py

if errorlevel 1 exit /b 1
echo.
echo Built: "%ROOT%\dist\Node Runner.exe"
