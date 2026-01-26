@echo off
setlocal
cd /d "%~dp0"

REM --- Always work from this folder (where the BAT lives)
set "ROOT=%~dp0"
set "VENV=%ROOT%.venv"
set "PY=%VENV%\Scripts\python.exe"

REM --- Create venv if missing
if not exist "%PY%" (
  echo Creating venv...
  python -m venv "%VENV%"
  if errorlevel 1 goto :fail
)

REM --- Install deps INSIDE venv (critical)
echo Installing dependencies in venv...
"%PY%" -m pip install --upgrade pip >nul
"%PY%" -m pip install "websockets>=12,<13"
if errorlevel 1 goto :fail

REM --- Run server
echo Starting...
"%PY%" "%ROOT%main.py"
if errorlevel 1 goto :fail

pause
exit /b 0

:fail
echo.
echo FAILED. Check errors above.
pause
exit /b 1
