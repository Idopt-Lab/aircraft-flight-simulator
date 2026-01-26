#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

VENV="$ROOT/.venv"
PY="$VENV/bin/python"

if [[ ! -x "$PY" ]]; then
  python3 -m venv "$VENV" || python -m venv "$VENV"
fi

"$PY" -m pip install --upgrade pip
"$PY" -m pip install "websockets>=12,<13"

MAIN=""
if [[ -f "$ROOT/main.py" ]]; then
  MAIN="$ROOT/main.py"
else
  MAIN="$(find "$ROOT" -type f -name "main.py" -print -quit || true)"
fi

[[ -z "$MAIN" ]] && exit 1

"$PY" "$MAIN"
read -n 1 -s -r -p "Press any key to close..."
echo
