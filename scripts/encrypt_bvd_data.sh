#!/usr/bin/env bash
# Writes data/bvd_data.tar.gz.enc: the situation report tables of the
# Bundibugyo virus disease outbreak and the tables derived from them which
# the experiments need, encrypted with the passphrase in ~/.bvd_data_key
# (the repository secret BVD_DATA_KEY holds the same passphrase). Rerun
# after adding reports, then commit the file.
set -euo pipefail
key="${BVD_DATA_KEY_FILE:-$HOME/.bvd_data_key}"
cd "$(dirname "$0")/.."
find "data/raw/Ebola CFR" -maxdepth 1 -name 'Varified*.csv' -print0 > /tmp/bvd_files.$$
if [ -d data/raw/bvd_derived ]; then find data/raw/bvd_derived -name '*.csv' -print0 >> /tmp/bvd_files.$$; fi
COPYFILE_DISABLE=1 tar --null -T /tmp/bvd_files.$$ -czf - |
  openssl enc -aes-256-cbc -pbkdf2 -iter 200000 -md sha256 -salt -pass "file:$key" \
    -out data/bvd_data.tar.gz.enc
rm -f /tmp/bvd_files.$$
ls -l data/bvd_data.tar.gz.enc
