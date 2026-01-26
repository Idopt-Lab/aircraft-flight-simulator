#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

VENV=".venv"
PY="$VENV/bin/python"

# Pick a Python
if command -v python3 >/dev/null 2>&1; then
  PYLAUNCH="python3"
elif command -v python >/dev/null 2>&1; then
  PYLAUNCH="python"
else
  echo "Python not found. Install it first:"
  echo "  brew install python"
  exit 1
fi

# Create venv if missing
if [ ! -x "$PY" ]; then
  "$PYLAUNCH" -m venv "$VENV"
fi

# Install deps
"$PY" -m pip install -U pip
"$PY" -m pip install "websockets>=12,<13"

# Run main
"$PY" main.py
