@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0"

where py >nul 2>nul
if %errorlevel%==0 (
  set "PYLAUNCH=py -3"
) else (
  where python >nul 2>nul
  if %errorlevel%==0 (
    set "PYLAUNCH=python"
  ) else (
    where winget >nul 2>nul
    if %errorlevel%==0 (
      echo Python not found. Installing via winget...
      winget install -e --id Python.Python.3.12 --source winget
      if errorlevel 1 (
        start "" "https://www.python.org/downloads/windows/"
        echo Install Python, then run this file again.
        pause
        exit /b 1
      )
    ) else (
      start "" "https://www.python.org/downloads/windows/"
      echo Install Python, then run this file again.
      pause
      exit /b 1
    )
    where py >nul 2>nul
    if %errorlevel%==0 (
      set "PYLAUNCH=py -3"
    ) else (
      where python >nul 2>nul
      if %errorlevel%==0 (
        set "PYLAUNCH=python"
      ) else (
        echo Python install did not complete. Run again after installing.
        pause
        exit /b 1
      )
    )
  )
)

set "ROOT=%~dp0"
set "VENV=%ROOT%.venv"
set "PY=%VENV%\Scripts\python.exe"

if not exist "%PY%" (
  "%PYLAUNCH%" -m venv "%VENV%"
  if errorlevel 1 goto :fail
)

"%PY%" -m pip install --upgrade pip
if errorlevel 1 goto :fail

"%PY%" -m pip install "websockets>=12,<13"
if errorlevel 1 goto :fail

set "MAIN="
if exist "%ROOT%main.py" set "MAIN=%ROOT%main.py"
if not defined MAIN (
  for /r "%ROOT%" %%F in (main.py) do (
    set "MAIN=%%F"
    goto :found
  )
)
:found
if not defined MAIN goto :fail

"%PY%" "%MAIN%"
if errorlevel 1 goto :fail

pause
exit /b 0

:fail
echo FAILED.
pause
exit /b 1
