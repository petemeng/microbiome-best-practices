#!/usr/bin/env bash
set -euo pipefail

project_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
raw_root="${project_root}/data/raw/articles-46-50"

download_checked() {
  local url="$1"
  local destination="$2"
  local expected_sha256="$3"

  mkdir -p "$(dirname "${destination}")"
  if [[ ! -s "${destination}" ]]; then
    curl -fL --retry 3 --retry-delay 2 "${url}" -o "${destination}"
  fi

  local observed_sha256
  observed_sha256="$(sha256sum "${destination}" | awk '{print $1}')"
  if [[ "${observed_sha256}" != "${expected_sha256}" ]]; then
    echo "SHA-256 mismatch: ${destination}" >&2
    echo "expected ${expected_sha256}" >&2
    echo "observed ${observed_sha256}" >&2
    exit 1
  fi
}

franzosa_commit="89a519d8c832008fbc6e650453e83e2f04858d02"
franzosa_base="https://raw.githubusercontent.com/borenstein-lab/microbiome-metabolome-curated-data/${franzosa_commit}/data/processed_data/FRANZOSA_IBD_2019"
franzosa_dir="${raw_root}/franzosa-ibd-2019"

download_checked "${franzosa_base}/genera.counts.tsv" "${franzosa_dir}/genera.counts.tsv" "f89aa7832be4ec75081eaa8f38ce4e80bdaf2b7c2e42f09f73e3e5eb35b58ae4"
download_checked "${franzosa_base}/metadata.tsv" "${franzosa_dir}/metadata.tsv" "f7396e3d6838b3b30f78b02bd568753757f84c956cd351966dbe654d50285376"
download_checked "${franzosa_base}/mtb.tsv" "${franzosa_dir}/mtb.tsv" "528b5e5953bd3697dd1ecf551d810d536c0679bd922e3fa3a6956c1412c6288c"
download_checked "${franzosa_base}/mtb.map.tsv" "${franzosa_dir}/mtb.map.tsv" "0dcdcce04a4e9b2b9b1632a410959baa4802ea9e14fc7c44f63bc17f699e5c65"

duran_commit="6db5e85cc5d442fd95fcdcb7250b72fa9e2ff900"
duran_base="https://raw.githubusercontent.com/ththi/Microbial-Interkingdom-Suppl/${duran_commit}/microb_interking_scripts/fig1/data"
duran_dir="${raw_root}/duran-interkingdom-2018"

download_checked "${duran_base}/otu_tab_filter_001_bac.txt" "${duran_dir}/otu_tab_filter_001_bac.txt" "ca90bd08ac7cc02f92bd50f8ac1c669aa7453ab01bd27d2d3045eadf3d57074e"
download_checked "${duran_base}/taxonomy_ref_bac.txt" "${duran_dir}/taxonomy_ref_bac.txt" "0e4060695620c8188157c35d4f97ad6bba3b8ec6b53a00ebcc597491c24ab424"
download_checked "${duran_base}/design_own_bac.txt" "${duran_dir}/design_own_bac.txt" "bd2f082b844621616060fa88b2a8c67ccf124126712a525f86fa3356b82264da"
download_checked "${duran_base}/otu_tab_filter_001_fun.txt" "${duran_dir}/otu_tab_filter_001_fun.txt" "5fbb95263f099adb7d8f051dd0da6bb4ec6950c6c012df5cad1148f1b82fad85"
download_checked "${duran_base}/taxonomy_ref_fun.txt" "${duran_dir}/taxonomy_ref_fun.txt" "1df2bf8750520a49e3b2935863641abbe51a1f280e8ee4ea5a08644afde088eb"
download_checked "${duran_base}/design_own_fun.txt" "${duran_dir}/design_own_fun.txt" "543325611c0a855e0269a09b9e652baae034301fe3824a6b3d19c47df9abe19f"
download_checked "${duran_base}/otu_tab_filter_001_oo.txt" "${duran_dir}/otu_tab_filter_001_oo.txt" "10fcfe028a57198716285e47610188cf66704c1b21ae91859cc511f940af0835"
download_checked "${duran_base}/taxonomy_ref_oo.txt" "${duran_dir}/taxonomy_ref_oo.txt" "baed37df5cd501528f12af4ef0e1823ecae096463951b9c5e8ea6df50f85b670"
download_checked "${duran_base}/design_own_oomyc.txt" "${duran_dir}/design_own_oomyc.txt" "ad132d9a741d7db1af218e9347c4588ce2d3433e4502ca62b8700a8e1a668965"

feast_commit="2f8f3df8051e0e08341f597a9f4693bfb76b3bf6"
feast_base="https://raw.githubusercontent.com/cozygene/FEAST/${feast_commit}/Data_files"
feast_dir="${raw_root}/feast-example"

download_checked "${feast_base}/otu_example.txt" "${feast_dir}/otu_example.txt" "1be44c28462fe4b67f30852698168b13dc6ce9feca2d5aaa0e5e2a081bbb30d3"
download_checked "${feast_base}/metadata_example.txt" "${feast_dir}/metadata_example.txt" "08cc3c0c10369e69ebe774b94a06c78597592228a350864df4bb541f1105bbae"

echo "Articles 46-50 raw inputs downloaded and verified under ${raw_root}"
