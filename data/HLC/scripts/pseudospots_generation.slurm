#!/bin/bash
#SBATCH --job-name=pseudospots
#SBATCH --output=../derived/pseudospots_generation_%j.out
#SBATCH --qos=test
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:10:00

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

mkdir -p "$SCRIPT_DIR/../derived"

if [ -d "$REPO_ROOT/.venv" ]; then
    source "$REPO_ROOT/.venv/bin/activate"
fi
python "$SCRIPT_DIR/pseudospots_generation.py"
