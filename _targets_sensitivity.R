# The sensitivity analyses of the mechanistic model set against its main fit
# (R/sensitivity.R): plots of the case fatality ratio, Rt, the eventual deaths
# and the reporting parameters of each fit.
#
#   targets::tar_make(script = "_targets_sensitivity.R", store = "_targets_sensitivity")
#
# needs the outputs of each fit (the artifact mechanistic-outputs of its run of
# mechanistic.yml) in fits_sensitivity/<fit>/, the fits named as
# sensitivity_variants ('main' for share15_tail), all fitted to the reports of
# the same date; sensitivity.yml downloads them. Writes its plots to
# outputs/sensitivity/.

library(targets)

tar_option_set(packages = c("dplyr", "ggplot2"))
tar_source("R/sensitivity.R")

dir <- "fits_sensitivity"
out <- "outputs/sensitivity"

list(
  tar_target(fit_files,
             list.files(dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE),
             format = "file"),
  tar_target(plots, {
    fit_files
    write_sensitivity_plots(dir, out)
  }, format = "file")
)
