# The comparison of nowcasts (R/comparison.R): at each fit date before the
# last report, the projections of the data-driven model, epinowcast and
# baselinenowcast of the counts which the reports one, two, ... weeks on
# showed, scored against them.
#
#   targets::tar_make(script = "_targets_comparison.R", store = "_targets_comparison")
#
# needs the reports under data/raw/
# (scripts/decrypt_bvd_data.sh) and the chains of the data-driven fits in
# fits/ (as _targets.R saves them), and writes its tables and plots to
# outputs/comparison. Each epinowcast fit (cases and deaths at one fit date)
# is a target of its own (epinowcast_<mmdd>), so that on GitHub Actions each
# runs as a job of its own, and saves its projections to fits_comparison/;
# with BVD_FITS_FROM_FILES=1 the targets take the files already there
# instead of fitting. BVD_DATES (mmdd, separated by commas) sets which fit
# dates are compared (by default all before the last report). With
# BVD_MECHANISTIC=1 the mechanistic model is compared too, from the
# projections of its chains already in fits_comparison/ (mechanistic_<mmdd>_<chain>.rds,
# made by project_mechanistic_chain on GitHub Actions) and its chains in
# fits_mechanistic/ (for its case fatality ratio).

library(targets)

tar_option_set(packages = c("dplyr", "tidyr", "ggplot2"))
tar_source("R")

from_files <- Sys.getenv("BVD_FITS_FROM_FILES", "") == "1"
with_mechanistic <- Sys.getenv("BVD_MECHANISTIC", "") == "1"
fit_dates <- as.Date(c("2026-08-01", "2026-08-09", "2026-08-19", "2026-08-26", "2026-09-02",
                       "2026-09-09", "2026-09-16", "2026-09-24", "2026-09-30", "2026-10-07"))
chosen <- trimws(strsplit(Sys.getenv("BVD_DATES", ""), ",")[[1]])
if(length(chosen[chosen != ""]) > 0)
  fit_dates <- fit_dates[format(fit_dates, "%m%d") %in% chosen]
out <- "outputs/comparison"

f_fit <- function(date) {
  label <- format(date, "%m%d")
  date <- as.character(date)
  s <- function(name) as.symbol(paste0(name, "_", label))
  path <- sprintf("fits_comparison/epinowcast_%s.rds", label)
  chains <- sprintf("fits/%s_%d.rds", label, 1:4)
  # the chains of the mechanistic fit which there are
  f_present <- function(paths) {
    paths <- paths[file.exists(paths)]
    if(length(paths) == 0)
      stop(paste("No chains of the mechanistic fit of", date))
    as.call(c(as.symbol("c"), as.list(paths)))
  }
  mechanistic <- if(with_mechanistic) list(
    tar_target_raw(paste0("mechanistic_", label),
                   f_present(sprintf("fits_comparison/mechanistic_%s_%d.rds", label, 1:4)),
                   format = "file"),
    tar_target_raw(paste0("mechanistic_chains_", label),
                   f_present(sprintf("fits_mechanistic/%s_%d.rds", label, 1:4)),
                   format = "file"),
    tar_target_raw(paste0("mechanistic_projections_", label),
                   bquote(project_mechanistic(.(s("mechanistic"))))))
  c(mechanistic, list(
    tar_target_raw(paste0("truth_", label), bquote(comparison_truth(snapshots, .(date)))),
    tar_target_raw(paste0("epinowcast_", label),
                   if(from_files) bquote(.(path)) else
                     bquote(fit_epinowcast(snapshots, .(date), .(path))),
                   format = "file"),
    tar_target_raw(paste0("chains_", label), as.call(c(as.symbol("c"), as.list(chains))), format = "file"),
    tar_target_raw(paste0("projections_", label), bquote({
      p <- c(readRDS(.(s("epinowcast")))$projections,
             fit_baselinenowcast(snapshots, .(date)),
             project_data_driven(.(s("chains")), snapshots, .(date), model_file = model_file),
             .(if(with_mechanistic) bquote(.(s("mechanistic_projections"))$projections)))
      c(p, comparison_ensemble(p))})),
    tar_target_raw(paste0("scores_", label), bquote(dplyr::bind_rows(lapply(
      .(s("projections")), comparison_score, truth = .(s("truth")))))),
    tar_target_raw(paste0("mixture_terms_", label),
                   bquote(comparison_mixture_terms(.(s("projections")), .(s("truth"))))),
    tar_target_raw(paste0("cfr_", label), bquote(dplyr::bind_rows(
      comparison_cfr_epinowcast(.(s("epinowcast"))),
      comparison_cfr_data_driven(.(s("chains")), snapshots, .(date), model_file = model_file),
      .(if(with_mechanistic) bquote(mechanistic_cfr(.(s("mechanistic_chains"))))))))))
}
f_all <- function(name) {
  as.call(c(as.symbol("list"), lapply(format(fit_dates, "%m%d"),
                                       function(l) as.symbol(paste0(name, "_", l)))))
}

c(
  list(
    tar_target(report_files,
               list.files("data/raw/Ebola CFR", pattern = "^Varified.*\\.csv$", full.names = TRUE),
               format = "file"),
    tar_target(snapshots, build_bvd_snapshots(report_files)),
    tar_target(model_file, "stan/nowcast_joint.stan", format = "file")),
  unlist(lapply(fit_dates, f_fit), recursive = FALSE),
  list(
    tar_target_raw("scores", bquote(dplyr::bind_rows(.(f_all("scores"))))),
    tar_target_raw("cfr", bquote(dplyr::bind_rows(.(f_all("cfr"))))),
    tar_target_raw("mixture_terms", bquote(dplyr::bind_rows(.(f_all("mixture_terms"))))),
    tar_target_raw("diagnostics", bquote(dplyr::bind_rows(lapply(.(f_all("epinowcast")),
                                                                  function(p) readRDS(p)$diagnostics)))),
    tar_target_raw("mechanistic_diagnostics", if(with_mechanistic)
      bquote(dplyr::bind_rows(lapply(.(f_all("mechanistic_projections")), `[[`, "diagnostics")))
      else quote(dplyr::tibble())),
    tar_target(summary_week, comparison_summary(scores, "week")),
    tar_target(summary_age, comparison_summary(
      dplyr::mutate(scores, age_band = comparison_age_band(.data$age)), "age_band")),
    tar_target(outputs, {
      dir.create(out, recursive = TRUE, showWarnings = FALSE)
      f_csv <- function(x, name) {
        path <- file.path(out, name)
        utils::write.csv(x, path, row.names = FALSE)
        path
      }
      c(f_csv(scores, "scores.csv"),
        f_csv(summary_week, "summary_by_week.csv"),
        f_csv(summary_age, "summary_by_age.csv"),
        f_csv(diagnostics, "epinowcast_diagnostics.csv"),
        f_csv(mechanistic_diagnostics, "mechanistic_diagnostics.csv"),
        f_csv(cfr, "cfr.csv"),
        f_csv(mixture_terms, "mixture_terms.csv"),
        plot_comparison_cfr(cfr, out),
        plot_comparison_weeks(summary_week, out),
        plot_comparison_pit(scores, out),
        plot_comparison_dates(scores, out))
    }, format = "file"))
)
