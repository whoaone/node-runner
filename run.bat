@echo off
REM Node Runner launcher - pins the bundled venv Python so PySide6's Qt
REM platform plugins resolve correctly. Running plain `python run.py`
REM can pick up a different Python (Anaconda, system) where PySide6 is
REM either missing or broken, producing the "no Qt platform plugin"
REM dialog. This launcher avoids that.

setlocal
set "HERE=%~dp0"
set "PY=%HERE%venv\Scripts\python.exe"

if not exist "%PY%" (
  echo Could not find venv Python at "%PY%".
  echo Create the venv first:  python -m venv venv  ^&^&  venv\Scripts\pip install -r requirements.txt
  pause
  exit /b 1
)

"%PY%" "%HERE%run.py" %*
