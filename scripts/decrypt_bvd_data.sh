#!/usr/bin/env bash
# Restores the files of data/bvd_data.tar.gz.enc under data/raw/, with the
# passphrase in the file named by BVD_DATA_KEY_FILE (default ~/.bvd_data_key).
set -euo pipefail
key="${BVD_DATA_KEY_FILE:-$HOME/.bvd_data_key}"
cd "$(dirname "$0")/.."
openssl enc -d -aes-256-cbc -pbkdf2 -iter 200000 -md sha256 -pass "file:$key" \
  -in data/bvd_data.tar.gz.enc | tar -xzf -
ls "data/raw/Ebola CFR" | wc -l
