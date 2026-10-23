#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
STAMP="$(date +%Y%m%d_%H%M%S)"
OUT_ROOT="${ROOT_DIR}/tests/smoke_outputs/${STAMP}"
DATA_DIR="${ROOT_DIR}/tests/data"

mkdir -p "${OUT_ROOT}/HiCAT" "${OUT_ROOT}/HORmon"

echo "Pipeline: ${ROOT_DIR}"
echo "Output:   ${OUT_ROOT}"
python "${ROOT_DIR}/HORmining.py" --version
bedtools --version
dot -V

python "${ROOT_DIR}/HORmining.py" \
  -bed "${DATA_DIR}/test.asat.bed" \
  -c "${DATA_DIR}/test.clusters.tsv" \
  -O "${OUT_ROOT}/HiCAT" \
  -m HiCAT

test -s "${OUT_ROOT}/HiCAT/filter.asat.bed"
test -s "${OUT_ROOT}/HiCAT/HiCAT/HiCAT.hor.summary.xls"

python "${ROOT_DIR}/HORmining.py" \
  -bed "${DATA_DIR}/test.asat.bed" \
  -c "${DATA_DIR}/test.clusters.tsv" \
  -O "${OUT_ROOT}/HORmon" \
  -m HORmon

test -s "${OUT_ROOT}/HORmon/filter.asat.bed"
test -s "${OUT_ROOT}/HORmon/graph/graph.HORdecomposition.stat.xls"
test -s "${OUT_ROOT}/HORmon/graph/0.png"

echo "HORmining smoke tests passed"
echo "Outputs retained at: ${OUT_ROOT}"
