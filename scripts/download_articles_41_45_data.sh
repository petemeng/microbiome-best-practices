#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cache_dir="${project_root}/data/raw/articles-41-45"
work_dir="${project_root}/.cache/articles-41-45-data"
mkdir -p "${cache_dir}" "${work_dir}"

download_checked() {
  local url="$1"
  local destination="$2"
  local expected_sha256="$3"
  if [[ ! -f "${destination}" ]] || [[ "$(sha256sum "${destination}" | cut -d' ' -f1)" != "${expected_sha256}" ]]; then
    curl -fL --retry 3 --retry-delay 2 "${url}" -o "${destination}"
  fi
  echo "${expected_sha256}  ${destination}" | sha256sum -c -
}

download_checked \
  "http://kronos.pharmacology.dal.ca/public_files/picrust/picrust2_tutorial_files/chemerin_16S.zip" \
  "${cache_dir}/chemerin_16S.zip" \
  "f57abdc069b6560f0ddb739cf9a341f5679b12f4947f4db2ad5e2cab893669e9"

download_checked \
  "https://zenodo.org/api/records/1146764/files/crc_xiang_results.tar.gz/content" \
  "${cache_dir}/crc_xiang_results.tar.gz" \
  "c16052495cc717069c670cabc9c1028750a71b7fd04ac3d28dbc1e560521d63e"
download_checked \
  "https://zenodo.org/api/records/1146764/files/crc_zhao_results.tar.gz/content" \
  "${cache_dir}/crc_zhao_results.tar.gz" \
  "2f479dab25d980f2d295244166165f5531b9d28e38896767d1acf6a08e88f91f"
download_checked \
  "https://zenodo.org/api/records/1146764/files/crc_zackular_results.tar.gz/content" \
  "${cache_dir}/crc_zackular_results.tar.gz" \
  "294305422293a2e87f9eb61d89612f49fb08b295d779b1cce7996b6d884c4498"

misurv_commit="692a2ac7079d2dbf581893b9a820c80cde0a7e31"
download_checked \
  "https://raw.githubusercontent.com/wg99526/MiSurvGit/${misurv_commit}/Data/otu.tab.txt" \
  "${cache_dir}/otu.tab.txt" \
  "2044c03631508ac8c1b799618d1e8a24ebc50cd1837eca5c22d8fecbbed9a960"
download_checked \
  "https://raw.githubusercontent.com/wg99526/MiSurvGit/${misurv_commit}/Data/tax.tab.txt" \
  "${cache_dir}/tax.tab.txt" \
  "3350c16af78eaaa3eb00e0221469a344e82cec0499f8a2f311133260030532b8"
download_checked \
  "https://raw.githubusercontent.com/wg99526/MiSurvGit/${misurv_commit}/Data/sam.dat.txt" \
  "${cache_dir}/sam.dat.txt" \
  "48eda70aaaf0ce13d1ff7508c7455e3dae2a6cb3d22aea87628e6b63070922ec"
download_checked \
  "https://raw.githubusercontent.com/wg99526/MiSurvGit/${misurv_commit}/Data/tree.tre" \
  "${cache_dir}/tree.tre" \
  "14fedf0c01814bd81e6a02f704d5f6fd19566415577ce5984a883065f737079d"

run_dir="$(mktemp -d "${work_dir}/run.XXXXXX")"
mkdir -p "${run_dir}/picrust2" "${run_dir}/microbiomehd" "${run_dir}/survival"
unzip -q "${cache_dir}/chemerin_16S.zip" -d "${run_dir}/picrust2"
for cohort in crc_xiang crc_zhao crc_zackular; do
  tar -xzf "${cache_dir}/${cohort}_results.tar.gz" -C "${run_dir}/microbiomehd"
done
cp "${cache_dir}/otu.tab.txt" "${cache_dir}/tax.tab.txt" \
  "${cache_dir}/sam.dat.txt" "${cache_dir}/tree.tre" "${run_dir}/survival/"

python3 "${project_root}/scripts/prepare_articles_41_45_data.py" \
  --picrust-source "${run_dir}/picrust2/chemerin_16S" \
  --picrust-archive "${cache_dir}/chemerin_16S.zip" \
  --microbiomehd-source "${run_dir}/microbiomehd" \
  --crc-xiang-archive "${cache_dir}/crc_xiang_results.tar.gz" \
  --crc-zhao-archive "${cache_dir}/crc_zhao_results.tar.gz" \
  --crc-zackular-archive "${cache_dir}/crc_zackular_results.tar.gz" \
  --survival-source "${run_dir}/survival" \
  --survival-commit "${misurv_commit}" \
  --output-root "${project_root}/data/small"
