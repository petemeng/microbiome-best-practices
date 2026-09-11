#!/usr/bin/env bash
set -euo pipefail

microeco_version="2.0.0"
microeco_expected="454a3b71ceeea86bdd54f475a12b5bac9c3b124f663b8dc6e560827fca89ebfd"
microeco_url="https://cran.r-project.org/src/contrib/Archive/microeco/microeco_${microeco_version}.tar.gz"
microeco_target="downloads/microeco_${microeco_version}.tar.gz"

sepp_expected="e252b83d7d5fbf2a9e14e594768e3578b33b557c34e22b6abc83b324689b1360"
sepp_url="https://data.qiime2.org/classifiers/sepp-ref-dbs/sepp-refs-gg-13-8.qza"
sepp_target="downloads/sepp-refs-gg-13-8.qza"

mkdir -p downloads
curl -fL --retry 3 --output "${microeco_target}" "${microeco_url}"
printf '%s  %s\n' "${microeco_expected}" "${microeco_target}" \
  | sha256sum --check -

curl -fL --retry 3 --output "${sepp_target}" "${sepp_url}"
printf '%s  %s\n' "${sepp_expected}" "${sepp_target}" \
  | sha256sum --check -
