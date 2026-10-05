#!/usr/bin/env bash
set -euo pipefail

if command -v mamba >/dev/null 2>&1; then
  exec mamba run -n ba_conjugate_workbench baworkbench
elif command -v conda >/dev/null 2>&1; then
  exec conda run -n ba_conjugate_workbench baworkbench
else
  echo "Neither mamba nor conda was found on PATH." >&2
  exit 1
fi
