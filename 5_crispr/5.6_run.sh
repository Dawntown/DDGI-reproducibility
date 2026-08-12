#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# Set DDGI_PYTHON to the Python executable from an environment with DDGI installed.
python_path="${DDGI_PYTHON:-/path/to/ddgi-env/bin/python}"

"${python_path}" 5.6_run_eval_hardeff.py "${1}" > "run_eval_hardeff_${1}.log" 2>&1
