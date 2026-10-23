# The analysis: the data-driven model of cases and deaths
# (stan/nowcast_joint.stan, R/nowcast.R), fitted to the reports up to each
# fit date. The last fit date is the last report; the earlier ones are those
# of the fits of the mechanistic model, at which the projections are checked
# against what the later reports showed.
#
#   targets::tar_make()
#
# needs the reports under data/raw/ (scripts/decrypt_bvd_data.sh) and writes
# its tables and plots to outputs/data-driven/. Each fit has four chains,
# each a target of its own (fit_<mmdd>_<chain>), so that on GitHub Actions
# each chain runs as a job of its own. Each chain saves its draws to fits/;
# with BVD_FITS_FROM_FILES=1 the chain targets take the files already there
# instead of fitting, and the rest of the pipeline (the projections from the
# saved draws, the checks of convergence, scores, tables and plots) is made
# from them.
#
# BVD_ITER_WARMUP and BVD_ITER_SAMPLING (default 1000 each) set the length of
# the chains, BVD_THREADS (default 4) the threads of each chain, and
# BVD_DATES (mmdd, separated by commas) which fits are made (by default all).

library(targets)

tar_option_set(packages = c("dplyr", "tidyr", "ggplot2"))
tar_source("R")

f_int <- function(name, default) as.integer(Sys.getenv(name, default))
iter_warmup <- f_int("BVD_ITER_WARMUP", "1000")
iter_sampling <- f_int("BVD_ITER_SAMPLING", "1000")
threads <- f_int("BVD_THREADS", "4")
from_files <- Sys.getenv("BVD_FITS_FROM_FILES", "") == "1"
n_chains <- 4L

date_last <- as.Date("2026-10-14")
fit_dates <- as.Date(c("2026-08-01", "2026-08-09", "2026-08-19", "2026-08-26", "2026-09-02",
                       "2026-09-09", "2026-09-16", "2026-09-24", "2026-09-30", "2026-10-07",
                       "2026-10-14"))
# BVD_DATES (fit dates as mmdd, separated by commas) keeps only those fits
chosen <- trimws(strsplit(Sys.getenv("BVD_DATES", ""), ",")[[1]])
if(length(chosen[chosen != ""]) > 0)
  fit_dates <- fit_dates[format(fit_dates, "%m%d") %in% chosen]
out <- "outputs/data-driven"

# the targets of one fit date
f_fit <- function(date) {
  label <- format(date, "%m%d")
  s <- function(name) as.symbol(paste0(name, "_", label))
  chains <- lapply(seq_len(n_chains), function(k) {
    path <- sprintf("fits/%s_%d.rds", label, k)
    command <- if(from_files) bquote(.(path)) else
      bquote(fit_nowcast_chain(.(s("data")), .(k), .(path), iter_warmup = .(iter_warmup),
                               iter_sampling = .(iter_sampling), chains = .(n_chains),
                               threads = .(threads), model_file = model))
    tar_target_raw(paste0("fit_", label, "_", k), command, format = "file")
  })
  chain_names <- lapply(seq_len(n_chains), function(k) as.symbol(paste0("fit_", label, "_", k)))
  Filter(Negate(is.null), c(list(
    tar_target_raw(paste0("analysis_", label), bquote(prepare_bvd_reports(snapshots, .(date)))),
    tar_target_raw(paste0("data_", label), bquote(nowcast_stan_data(.(s("analysis")))))),
    chains,
    list(
      tar_target_raw(paste0("joined_", label),
                     bquote(nowcast_join_chains(c(..(chain_names))), splice = TRUE)),
      tar_target_raw(paste0("convergence_", label),
                     bquote(nowcast_convergence(.(s("joined")), .(date)))),
      tar_target_raw(paste0("parameters_", label),
                     bquote(nowcast_parameters(.(s("joined")), .(date)))),
      if(date < date_last) tar_target_raw(paste0("vantage_", label),
                     bquote(nowcast_vantage(.(s("joined")), .(s("data")), .(s("analysis")),
                                            snapshots, .(date_last)))),
      tar_target_raw(paste0("eventual_", label),
                     bquote(nowcast_eventual(.(s("joined")), .(s("data")), .(s("analysis"))))))))
}
f_all <- function(name, dates = fit_dates) {
  as.call(c(as.symbol("list"), lapply(format(dates, "%m%d"),
                                       function(l) as.symbol(paste0(name, "_", l)))))
}

c(
  list(
    tar_target(report_files,
               list.files("data/raw/Ebola CFR", pattern = "^Varified.*\\.csv$", full.names = TRUE),
               format = "file"),
    tar_target(snapshots, build_bvd_snapshots(report_files)),
    tar_target(model_file, "stan/nowcast_joint.stan", format = "file"),
    # compiled once here, before the fits use it
    tar_target(model, {nowcast_model(model_file); model_file}, deployment = "main")),
  unlist(lapply(fit_dates, f_fit), recursive = FALSE),
  list(
    tar_target_raw("vantage", f_all("vantage", fit_dates[fit_dates < date_last])),
    tar_target_raw("eventual", f_all("eventual")),
    tar_target_raw("convergence", bquote(dplyr::bind_rows(.(f_all("convergence"))))),
    tar_target_raw("parameters", bquote(dplyr::bind_rows(.(f_all("parameters"))))),
    tar_target(scores, dplyr::bind_rows(lapply(vantage, nowcast_score,
                                               date_last = date_last))),
    tar_target(outputs, {
      dir.create(out, recursive = TRUE, showWarnings = FALSE)
      f_csv <- function(x, name) {
        path <- file.path(out, name)
        utils::write.csv(x, path, row.names = FALSE)
        path
      }
      by_day <- dplyr::bind_rows(lapply(eventual, `[[`, "by_day"))
      c(f_csv(convergence, "convergence.csv"),
        f_csv(parameters, "parameters.csv"),
        f_csv(scores, "scores.csv"),
        f_csv(dplyr::bind_rows(lapply(vantage, `[[`, "totals")), "vantage_totals.csv"),
        f_csv(dplyr::bind_rows(lapply(eventual, `[[`, "totals")), "eventual_totals.csv"),
        f_csv(dplyr::bind_rows(lapply(eventual, `[[`, "cfr")), "cfr.csv"),
        f_csv(by_day, "eventual_by_day.csv"),
        write_nowcast_parameters(parameters, by_day, out),
        if(length(vantage) > 0) c(plot_nowcast_backtests(vantage, date_last, out),
                                  plot_nowcast_vantage(vantage, out),
                                  unlist(lapply(vantage, plot_nowcast_vantage_days, out = out))),
        plot_nowcast_cfr(eventual, out))
    }, format = "file"))
)
