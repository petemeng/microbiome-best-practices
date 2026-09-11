#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
picrust_prefix="${project_root}/.tools/picrust2-2.6.3"
input_dir="${project_root}/data/small/picrust2-chemerin"
output_dir="${input_dir}/picrust2-out-v2.6.3"
result_dir="${project_root}/results/41-picrust2"
log_file="${result_dir}/picrust2-run.log"
complete_marker="${output_dir}/.complete"

if [[ ! -x "${picrust_prefix}/bin/picrust2_pipeline.py" ]]; then
  echo "PICRUSt2 2.6.3 was not found at ${picrust_prefix}." >&2
  exit 1
fi
if [[ -f "${complete_marker}" ]]; then
  echo "PICRUSt2 output is already complete: ${output_dir}"
  exit 0
fi
if [[ -e "${output_dir}" ]]; then
  echo "Refusing to overwrite an incomplete output directory: ${output_dir}" >&2
  exit 1
fi

mkdir -p "${result_dir}"
export PATH="${picrust_prefix}/bin:${PATH}"
export CONDA_PREFIX="${picrust_prefix}"
export PYTHONUNBUFFERED=1

{
  date --iso-8601=seconds
  picrust2_pipeline.py --version
  echo "pipeline-start"
  picrust2_pipeline.py \
    --study_fasta "${input_dir}/representative-sequences.fasta" \
    --input "${input_dir}/table.biom" \
    --output "${output_dir}" \
    --processes 8 \
    --stratified \
    --coverage \
    --remove_intermediate \
    --verbose
  find "${output_dir}" -type f -print0 | sort -z | xargs -0 sha256sum \
    > "${result_dir}/picrust2-output-sha256.tsv"
  touch "${complete_marker}"
  date --iso-8601=seconds
} > "${log_file}" 2>&1

echo "PICRUSt2 completed. Log: ${log_file}"
