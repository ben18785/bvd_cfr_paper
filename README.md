# The case fatality ratio and Rt of the Bundibugyo virus disease outbreak

The analysis is a targets pipeline (to come). The write-up is an Overleaf
project, cloned into `paper/` as its own repository (it is not tracked here):

    git clone https://git.overleaf.com/<project-id> paper

## Data

The situation reports are kept only encrypted, in `data/bvd_data.tar.gz.enc`.
`scripts/decrypt_bvd_data.sh` restores them under `data/raw/` with the
passphrase in `~/.bvd_data_key`; `scripts/encrypt_bvd_data.sh` rebuilds the
archive after reports are added. Neither the passphrase nor `data/raw/` is
ever committed.
