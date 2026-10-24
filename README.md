# The case fatality ratio and Rt of the Bundibugyo virus disease outbreak

The analysis is two targets pipelines, one for each model, each fitted to
the reports up to a series of dates from 1 August to 14 October:

- `_targets.R`: the data-driven model (`stan/nowcast_joint.stan`,
  `R/nowcast.R`); outputs in `outputs/data-driven/`.
- `_targets_mechanistic.R`: the mechanistic model (the package
  incidenceinflation 2.0, `R/mechanistic.R`); outputs in
  `outputs/mechanistic/`. Run with
  `targets::tar_make(script = "_targets_mechanistic.R", store = "_targets_mechanistic")`.

The fits run on GitHub Actions (`.github/workflows/data-driven.yml` and
`mechanistic.yml`), one job per chain, started by a push which changes
`runs/data-driven.txt` or `runs/mechanistic.txt` (the comments at the top of
each workflow say what these files set). The write-up is an Overleaf
project, cloned into `paper/` as its own repository (it is not tracked here):

    git clone https://git.overleaf.com/<project-id> paper

## Data

The situation reports are kept only encrypted, in `data/bvd_data.tar.gz.enc`.
`scripts/decrypt_bvd_data.sh` restores them under `data/raw/` with the
passphrase in `~/.bvd_data_key`; `scripts/encrypt_bvd_data.sh` rebuilds the
archive after reports are added. Neither the passphrase nor `data/raw/` is
ever committed.
