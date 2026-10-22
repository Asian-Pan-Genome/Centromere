#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")/.."
config_path="${1:-tests/data/test.config.yaml}"
if [[ ! -f "${config_path}" ]]; then
  echo "Missing ${config_path}. Supply a config path or use tests/data/test.config.yaml." >&2
  exit 2
fi
if grep -Fq 'executable: /path/to/RepeatMasker' "${config_path}"; then
  echo "Edit RepeatMasker executable and library paths in ${config_path} for this Linux installation." >&2
  exit 2
fi
if grep -Fq '/path/to/chm13_hsat/' "${config_path}"; then
  echo "Edit hsat23 Perl script and k-mer paths in ${config_path} for this Linux installation." >&2
  exit 2
fi

(cd tests/data && md5sum --check md5.txt)

# The workflow requires plain FASTA files for validation, indexing and reference
# anchoring. Keep only the compressed test inputs in Git; expand them locally.
expand_test_fasta() {
  local compressed="$1" uncompressed="$2" temporary
  if [[ -f "${uncompressed}" ]] && gzip -cd -- "${compressed}" | cmp -s - "${uncompressed}"; then
    return 0
  fi
  temporary=$(mktemp "${uncompressed}.tmp.XXXXXX")
  if gzip -cd -- "${compressed}" > "${temporary}"; then
    mv -f -- "${temporary}" "${uncompressed}"
  else
    rm -f -- "${temporary}"
    echo "Failed to decompress ${compressed}" >&2
    return 1
  fi
}

expand_test_fasta 'tests/data/C001-CHA-E01#Mat#chr2.fa.gz' \
  'tests/data/C001-CHA-E01#Mat#chr2.fasta'
expand_test_fasta 'tests/data/CHM13v2.chr2.fa.gz' \
  'tests/data/CHM13v2.chr2.fa'

python -m unittest discover -s tests -v
snakemake --configfile "${config_path}" --cores 16 --resources mem_mb=64000 --rerun-incomplete
