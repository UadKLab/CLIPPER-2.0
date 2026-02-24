#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

run_experiment() {
  local name="$1"
  shift
  echo ""
  echo "=== Running ${name} ==="
  echo "Command: python -m clipper.run $*"
  PYTHONPATH=. python -m clipper.run "$@"
  echo "=== ${name} completed ==="
}

require_file() {
  local path="$1"
  if [[ ! -f "$path" ]]; then
    echo "Required input file not found: $path" >&2
    exit 1
  fi
}

require_file "tests/HUNTER_clean_100.xlsx"
require_file "tests/cond_HUNTER.txt"
require_file "tests/tests for paper/GluC/GluC_peptide_groups_ci_1pct.xlsx"
require_file "tests/tests for paper/GluC/cond_GluC.txt"

run_experiment \
  "HUNTER clean" \
  -i "tests/HUNTER_clean_100.xlsx" \
  -cf "tests/cond_HUNTER.txt" \
  -o "ci_hunter" \
  -stat -spw -sig all -vis

run_experiment \
  "GluC" \
  -i "tests/tests for paper/GluC/GluC_peptide_groups_ci_1pct.xlsx" \
  -cf "tests/tests for paper/GluC/cond_GluC.txt" \
  -o "ci_gluc" \
  -stat -spw -sig all -vis -nx -logo all
