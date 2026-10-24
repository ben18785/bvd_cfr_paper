# The analysis with the mechanistic model (R/mechanistic.R, the package
# incidenceinflation 2.0), fitted to the reports up to each fit date: the
# dates of the data-driven model's fits (_targets.R).
#
#   targets::tar_make(script = "_targets_mechanistic.R", store = "_targets_mechanistic")
#
# needs the reports under data/raw/ (scripts/decrypt_bvd_data.sh) and writes
# its tables and plots to outputs/mechanistic/. Each fit has four chains, each
# a target of its own (fit_<mmdd>_<chain>), so that on GitHub Actions each
# chain runs as a job of its own. A chain saves its draws to
# fits_mechanistic/ and what is needed to carry it on to states_mechanistic/;
# if the latter holds the chain's state from an earlier run, the chain is
# carried on from it. With BVD_FITS_FROM_FILES=1 the chain targets take the
# files already in fits_mechanistic/ instead of fitting, and the rest of the
# pipeline (the checks of convergence, tables and plots) is made from them.
#
# BVD_ITERATIONS (default 3000) sets the iterations of each chain in all,
# BVD_TIME_LIMIT (minutes, default none) when a chain stops to be carried on
# later, BVD_RIDGE=1 adds the joint move of the fraction counted on the day of
# death and the delay to a death being recorded, and BVD_DATES (mmdd,
# separated by commas) which fits are made (by default all).

library(targets)

tar_option_set(packages = c("dplyr", "tidyr", "ggplot2"))
tar_source(c("R/data.R", "R/nowcast.R", "R/mechanistic.R"))

iterations <- as.integer(Sys.getenv("BVD_ITERATIONS", "3000"))
time_limit <- suppressWarnings(as.numeric(Sys.getenv("BVD_TIME_LIMIT", "")))
if(is.na(time_limit)) time_limit <- NULL
ridge <- Sys.getenv("BVD_RIDGE", "") == "1"
from_files <- Sys.getenv("BVD_FITS_FROM_FILES", "") == "1"
n_chains <- 4L

date_last <- as.Date("2026-10-14")
fit_dates <- as.Date(c("2026-08-01", "2026-08-09", "2026-08-19", "2026-08-26", "2026-09-02",
                       "2026-09-09", "2026-09-16", "2026-09-24", "2026-09-30", "2026-10-07",
                       "2026-10-14"))
chosen <- trimws(strsplit(Sys.getenv("BVD_DATES", ""), ",")[[1]])
if(length(chosen[chosen != ""]) > 0)
  fit_dates <- fit_dates[format(fit_dates, "%m%d") %in% chosen]
out <- "outputs/mechanistic"

# the targets of one fit date
f_fit <- function(date) {
  label <- format(date, "%m%d")
  s <- function(name) as.symbol(paste0(name, "_", label))
  chains <- lapply(seq_len(n_chains), function(k) {
    path <- sprintf("fits_mechanistic/%s_%d.rds", label, k)
    state <- sprintf("states_mechanistic/%s_%d.rds", label, k)
    # the chains a run did not make (one which failed) are left out later
    if(from_files) return(tar_target_raw(paste0("fit_", label, "_", k), bquote(.(path))))
    tar_target_raw(paste0("fit_", label, "_", k),
                   bquote(fit_mechanistic_chain(.(s("analysis")), .(k), .(path), .(state),
                                                iterations = .(iterations),
                                                time_limit = .(time_limit), ridge = .(ridge))),
                   format = "file")
  })
  chain_names <- lapply(seq_len(n_chains), function(k) as.symbol(paste0("fit_", label, "_", k)))
  c(list(tar_target_raw(paste0("analysis_", label),
                        bquote(mechanistic_analysis(snapshots, .(date))))),
    chains,
    list(
      tar_target_raw(paste0("joined_", label),
                     bquote(mechanistic_join_chains(c(..(chain_names))), splice = TRUE)),
      tar_target_raw(paste0("convergence_", label),
                     bquote(mechanistic_convergence(.(s("joined")), .(date)))),
      tar_target_raw(paste0("parameters_", label),
                     bquote(mechanistic_parameters(.(s("joined")), .(date)))),
      tar_target_raw(paste0("eventual_", label),
                     bquote(mechanistic_eventual(.(s("joined")), snapshots, .(date_last))))))
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
    tar_target(snapshots, build_bvd_snapshots(report_files))),
  unlist(lapply(fit_dates, f_fit), recursive = FALSE),
  list(
    tar_target_raw("eventual", f_all("eventual")),
    tar_target_raw("convergence", bquote(dplyr::bind_rows(.(f_all("convergence"))))),
    tar_target_raw("parameters", bquote(dplyr::bind_rows(.(f_all("parameters"))))),
    tar_target(outputs, {
      dir.create(out, recursive = TRUE, showWarnings = FALSE)
      f_csv <- function(x, name) {
        path <- file.path(out, name)
        utils::write.csv(x, path, row.names = FALSE)
        path
      }
      f_bind <- function(name) dplyr::bind_rows(lapply(eventual, `[[`, name))
      c(f_csv(convergence, "convergence.csv"),
        f_csv(parameters, "parameters.csv"),
        f_csv(f_bind("totals"), "eventual_totals.csv"),
        f_csv(f_bind("by_day"), "eventual_by_day.csv"),
        f_csv(f_bind("days"), "by_onset_day.csv"),
        write_mechanistic_parameters(parameters, f_bind("days"), out),
        plot_mechanistic_backtests(eventual, out),
        plot_mechanistic_totals(eventual, out),
        plot_mechanistic_days(eventual, out))
    }, format = "file"))
)
