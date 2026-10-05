#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="ba_conjugate_workbench"
if command -v mamba >/dev/null 2>&1; then
  TOOL="mamba"
elif command -v conda >/dev/null 2>&1; then
  TOOL="conda"
else
  echo "Neither mamba nor conda was found on PATH." >&2
  exit 1
fi

if "$TOOL" env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  "$TOOL" env update -n "$ENV_NAME" -f environment.yml --prune
else
  "$TOOL" env create -f environment.yml
fi
"$TOOL" run -n "$ENV_NAME" python -m pip install -e '.[test]'
echo "Environment $ENV_NAME is ready."
