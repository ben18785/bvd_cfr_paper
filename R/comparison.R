# The comparison of nowcasts (_targets_comparison.R). At each fit date, each
# model projects the counts of cases and of deaths by onset day which later
# reports will show: the report about one, two, ... weeks on. The
# projections are scored against what those reports showed, on the onset days
# whose counts did not fall between the fit date and the report.
#
# Each model's projections take one form: a named list of 'model', 'series'
# ("cases" or "deaths"), 'date_now', 'onsets' (the onset days projected) and
# 'draws', a list with one matrix (draws by onset day) of the projected counts
# for each report, named by its date.

#' The settings of the comparison
comparison_settings <- function() {
  list(date_first = as.Date("2026-05-11"),
       date_last = as.Date("2026-10-14"),
       # the onset days scored: those of the last six weeks before the fit date
       max_age = 41L,
       weeks = 1:6,
       n_draws = 1000L,
       # epinowcast's growth rate has its default, a separate effect for each
       # day: a weekly random walk instead left 9 of the 20 fits unconverged,
       # against 2
       epinowcast = list(max_delay = 60L, chains = 4L, iter_warmup = 1000L,
                         iter_sampling = 1000L, adapt_delta = 0.95,
                         max_treedepth = 12L,
                         reference = ~ 1 + rw(week), expectation = ~ 0 + (1 | day)),
       # baselinenowcast takes only onset days whose reports are seen from the
       # day of onset, so from the first report
       baselinenowcast = list(max_delay = 28L, start = as.Date("2026-07-12"),
                              scale_factor = 3, prop_delay = 0.5))
}

#' The reports scored for a fit date: the last report on or before each
#' whole number of weeks after it
#'
#' @param snapshots as build_bvd_snapshots
#' @param date_now the fit date
#' @param settings as comparison_settings
#' @return a tibble of week and date_report
comparison_reports <- function(snapshots, date_now, settings = comparison_settings()) {
  dates <- sort(unique(snapshots$date_reported))
  dates <- dates[dates > as.Date(date_now) & dates <= settings$date_last]
  out <- dplyr::bind_rows(lapply(settings$weeks, function(k) {
    d <- dates[dates <= as.Date(date_now) + 7 * k]
    if(length(d) == 0) NULL else dplyr::tibble(week = k, date_report = max(d))
  }))
  dplyr::distinct(out, .data$date_report, .keep_all = TRUE)
}

#' The counts of each onset day in each report up to a date, carried forward
#' to the reports which do not show it
#'
#' @param snapshots as build_bvd_snapshots
#' @param column cases_reported or deaths_reported
#' @param date_to the last report taken
#' @return a tibble of date_onset, date_reported and count
comparison_counts <- function(snapshots, column, date_to) {
  snapshots %>%
    dplyr::filter(.data$date_reported <= as.Date(date_to)) %>%
    dplyr::select("date_onset", "date_reported", count = dplyr::all_of(column)) %>%
    tidyr::complete(.data$date_onset, .data$date_reported) %>%
    dplyr::arrange(.data$date_onset, .data$date_reported) %>%
    dplyr::group_by(.data$date_onset) %>%
    tidyr::fill("count", .direction = "down") %>%
    dplyr::ungroup() %>%
    dplyr::mutate(count = dplyr::coalesce(.data$count, 0L))
}

#' The reports up to a fit date as the nowcasting packages take them: a
#' record removed by a later report (up to the fit date) is taken never to
#' have been there, so that no count falls
#'
#' @param snapshots as build_bvd_snapshots
#' @param column cases_reported or deaths_reported
#' @param date_now the fit date
#' @return a tibble of date_onset, date_reported and count
comparison_monotone <- function(snapshots, column, date_now) {
  comparison_counts(snapshots, column, date_now) %>%
    dplyr::arrange(.data$date_onset, dplyr::desc(.data$date_reported)) %>%
    dplyr::group_by(.data$date_onset) %>%
    dplyr::mutate(count = cummin(.data$count)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$date_onset, .data$date_reported)
}

#' What the reports scored for a fit date showed
#'
#' @param snapshots as build_bvd_snapshots
#' @param date_now the fit date
#' @param settings as comparison_settings
#' @return a tibble by series, week, report and onset day (of the last
#'   max_age days before the fit date): age (days from onset to the fit
#'   date), now (the count at the fit date), shown (in the report), and
#'   no_fall (whether the count did not fall from the fit date to the report)
comparison_truth <- function(snapshots, date_now, settings = comparison_settings()) {
  date_now <- as.Date(date_now)
  reports <- comparison_reports(snapshots, date_now, settings)
  onsets <- seq(date_now - settings$max_age, date_now, by = "day")
  dplyr::bind_rows(lapply(c(cases = "cases_reported", deaths = "deaths_reported"), function(column) {
    counts <- comparison_counts(snapshots, column, settings$date_last) %>%
      dplyr::filter(.data$date_onset %in% onsets)
    now <- counts %>% dplyr::filter(.data$date_reported <= date_now) %>%
      dplyr::group_by(.data$date_onset) %>%
      dplyr::summarise(now = .data$count[which.max(.data$date_reported)])
    dplyr::bind_rows(lapply(seq_len(nrow(reports)), function(i) {
      counts %>%
        dplyr::filter(.data$date_reported > date_now, .data$date_reported <= reports$date_report[i]) %>%
        dplyr::left_join(now, by = "date_onset") %>%
        dplyr::mutate(now = dplyr::coalesce(.data$now, 0L)) %>%
        dplyr::group_by(.data$date_onset) %>%
        dplyr::summarise(now = .data$now[1],
                         shown = .data$count[which.max(.data$date_reported)],
                         no_fall = all(diff(c(.data$now[1], .data$count)) >= 0)) %>%
        dplyr::mutate(week = reports$week[i], date_report = reports$date_report[i])
    }))
  }), .id = "series") %>%
    dplyr::mutate(date_now = date_now, age = as.integer(date_now - .data$date_onset)) %>%
    dplyr::select("series", "date_now", "week", "date_report", "date_onset", "age",
                  "now", "shown", "no_fall")
}

#' Keeps n of the draws, evenly spaced
#' @param m a matrix of draws by onset day
#' @param n the draws kept
comparison_thin <- function(m, n) {
  m[unique(round(seq(1, nrow(m), length.out = min(n, nrow(m))))), , drop = FALSE]
}

#' Fits epinowcast to the reports up to a fit date, and projects the reports
#' scored
#'
#' The reports are placed on a daily grid of report dates. The new count of
#' a report day enters the likelihood only if the day before also had a
#' report (or it is the day of onset): otherwise it lumps several days
#' together. So the onset days before the first report enter only through
#' their later reports. The grid runs on to the last report, or to the
#' longest delay after the fit date if that is later, its days after the fit
#' date out of the likelihood, so that the fit's predictions of them are the
#' projections. The delay is lognormal, its parameters changing as a random
#' walk by week of onset, with effects of the day of the week of the report;
#' the growth of the expected counts has a separate effect for each day. The eventual
#' count of each onset day is its count in the last report up to the fit
#' date and the predicted new counts after it, up to the longest delay.
#'
#' @param snapshots as build_bvd_snapshots
#' @param date_now the fit date
#' @param path where the projections are saved
#' @param settings as comparison_settings
#' @return path; the file holds a named list: 'projections' (one for cases
#'   and one for deaths, as described above) and 'diagnostics'
fit_epinowcast <- function(snapshots, date_now, path, settings = comparison_settings()) {
  date_now <- as.Date(date_now)
  s <- settings$epinowcast
  reports <- comparison_reports(snapshots, date_now, settings)
  report_days <- sort(unique(snapshots$date_reported[snapshots$date_reported <= date_now]))
  onsets <- seq(settings$date_first, date_now, by = "day")
  scored <- seq(date_now - settings$max_age, date_now, by = "day")
  out <- lapply(c(cases = "cases_reported", deaths = "deaths_reported"), function(column) {
    counts <- comparison_monotone(snapshots, column, date_now) %>%
      dplyr::filter(.data$date_onset >= settings$date_first)
    grid <- tidyr::crossing(reference_date = onsets,
                            report_date = seq(settings$date_first,
                                              max(settings$date_last, date_now + s$max_delay), by = "day")) %>%
      dplyr::filter(.data$report_date >= .data$reference_date,
                    .data$report_date <= .data$reference_date + s$max_delay) %>%
      dplyr::left_join(dplyr::select(counts, reference_date = "date_onset",
                                     report_date = "date_reported", confirm = "count"),
                       by = c("reference_date", "report_date")) %>%
      dplyr::arrange(.data$reference_date, .data$report_date) %>%
      dplyr::group_by(.data$reference_date) %>%
      tidyr::fill("confirm", .direction = "down") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(confirm = dplyr::coalesce(.data$confirm, 0L),
                    .observed = .data$report_date <= date_now & .data$report_date %in% report_days &
                      (.data$report_date == .data$reference_date |
                         (.data$report_date - 1) %in% report_days))
    pobs <- epinowcast::enw_preprocess_data(data.table::as.data.table(grid), max_delay = s$max_delay)
    start <- Sys.time()
    fit <- epinowcast::epinowcast(
      pobs,
      reference = epinowcast::enw_reference(s$reference, distribution = "lognormal", data = pobs),
      report = epinowcast::enw_report(~ (1 | day_of_week), data = pobs),
      expectation = epinowcast::enw_expectation(s$expectation, data = pobs),
      obs = epinowcast::enw_obs(family = "negbin", observation_indicator = ".observed", data = pobs),
      fit = epinowcast::enw_fit_opts(chains = s$chains, parallel_chains = s$chains,
                                     iter_warmup = s$iter_warmup, iter_sampling = s$iter_sampling,
                                     adapt_delta = s$adapt_delta,
                                     max_treedepth = s$max_treedepth, pp = TRUE,
                                     show_messages = FALSE, refresh = 0))
    minutes <- as.numeric(difftime(Sys.time(), start, units = "mins"))
    # the predicted new count of every cell of the grid, in the grid's order
    pp <- unclass(posterior::as_draws_matrix(fit$fit[[1]]$draws("pp_obs")))
    cells <- pobs$new_confirm[[1]]
    stopifnot(ncol(pp) == nrow(cells))
    now <- counts %>% dplyr::filter(.data$date_onset %in% scored) %>%
      dplyr::group_by(.data$date_onset) %>%
      dplyr::summarise(now = .data$count[which.max(.data$date_reported)])
    now <- dplyr::left_join(dplyr::tibble(date_onset = scored), now, by = "date_onset") %>%
      dplyr::mutate(now = dplyr::coalesce(.data$now, 0L))
    draws <- lapply(reports$date_report, function(r) {
      keep <- which(as.Date(cells$reference_date) %in% scored &
                      as.Date(cells$report_date) > date_now & as.Date(cells$report_date) <= r)
      add <- matrix(0, nrow(pp), length(scored))
      if(length(keep) > 0) {
        by_onset <- t(rowsum(t(pp[, keep, drop = FALSE]),
                             as.character(as.Date(cells$reference_date[keep]))))
        add[, match(as.Date(colnames(by_onset)), scored)] <- by_onset
      }
      comparison_thin(sweep(add, 2, now$now, "+"), settings$n_draws)
    })
    names(draws) <- as.character(reports$date_report)
    # the eventual count of each onset day: its count in the last report up
    # to the fit date, and the predicted new counts after it (up to the
    # longest delay), as for the projections
    latest <- counts %>% dplyr::group_by(.data$date_onset) %>%
      dplyr::summarise(now = .data$count[which.max(.data$date_reported)])
    latest <- dplyr::coalesce(latest$now[match(onsets, latest$date_onset)], 0L)
    later <- which(as.Date(cells$report_date) > date_now)
    eventual <- matrix(0, nrow(pp), length(onsets), dimnames = list(NULL, as.character(onsets)))
    if(length(later) > 0) {
      by_onset <- t(rowsum(t(pp[, later, drop = FALSE]), as.character(as.Date(cells$reference_date[later]))))
      eventual[, colnames(by_onset)] <- by_onset
    }
    eventual <- sweep(eventual, 2, latest, "+")
    eventual <- comparison_thin(eventual[, as.character(onsets), drop = FALSE], settings$n_draws)
    list(projections = list(model = "epinowcast", series = sub("_reported", "", column),
                            date_now = date_now, onsets = scored, draws = draws),
         eventual = list(series = sub("_reported", "", column), date_now = date_now,
                         onsets = onsets, draws = eventual),
         diagnostics = dplyr::tibble(model = "epinowcast", series = sub("_reported", "", column),
                                     date_now = date_now, max_rhat = fit$max_rhat,
                                     divergent = fit$divergent_transitions,
                                     at_max_treedepth = fit$no_at_max_treedepth,
                                     minutes = minutes))
  })
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(projections = lapply(out, `[[`, "projections"),
               eventual = lapply(out, `[[`, "eventual"),
               diagnostics = dplyr::bind_rows(lapply(out, `[[`, "diagnostics"))), path)
  path
}

#' Projections of baselinenowcast at a fit date
#'
#' The chain-ladder-like baseline: the delay is estimated from the recent
#' complete reports, and its uncertainty from how well it would have done at
#' earlier dates (a negative binomial by horizon). It takes onset days from
#' the first report onwards, on a daily grid of report dates with the counts
#' carried forward, up to a maximum delay; onset days before then or older
#' than the maximum delay are projected not to change. The count of a later
#' report is that of the fit date plus the projected new counts of the delays
#' it has reached; the uncertainty fitted to the whole of what is to come is
#' used for that part of it.
#'
#' @param snapshots as build_bvd_snapshots
#' @param date_now the fit date
#' @param settings as comparison_settings
#' @param seed for the draws
#' @return a list of two projections (cases and deaths), or NULL for a series
#'   if there are too few onset days to fit
fit_baselinenowcast <- function(snapshots, date_now, settings = comparison_settings(), seed = 1) {
  date_now <- as.Date(date_now)
  s <- settings$baselinenowcast
  reports <- comparison_reports(snapshots, date_now, settings)
  scored <- seq(date_now - settings$max_age, date_now, by = "day")
  onsets <- seq(s$start, date_now, by = "day")
  out <- lapply(c(cases = "cases_reported", deaths = "deaths_reported"), function(column) {
    counts <- comparison_monotone(snapshots, column, date_now)
    wide <- counts %>%
      dplyr::filter(.data$date_onset %in% onsets, .data$date_reported >= s$start) %>%
      dplyr::right_join(tidyr::crossing(date_onset = onsets,
                                        date_reported = seq(s$start, date_now, by = "day")),
                        by = c("date_onset", "date_reported")) %>%
      dplyr::arrange(.data$date_onset, .data$date_reported) %>%
      dplyr::group_by(.data$date_onset) %>%
      tidyr::fill("count", .direction = "down") %>%
      dplyr::ungroup() %>%
      dplyr::mutate(count = dplyr::coalesce(.data$count, 0L))
    cumulative <- matrix(wide$count, length(onsets), byrow = TRUE)
    tri <- matrix(NA_real_, length(onsets), s$max_delay + 1)
    for(i in seq_along(onsets)) for(d in 0:s$max_delay) {
      j <- i + d
      if(j > ncol(cumulative)) next
      tri[i, d + 1] <- cumulative[i, j] - (if(d == 0) 0 else cumulative[i, j - 1])
    }
    rt <- baselinenowcast::as_reporting_triangle(tri, reference_dates = onsets, delays_unit = "days")
    tv <- tryCatch(baselinenowcast::allocate_reference_times(
      rt, scale_factor = s$scale_factor, prop_delay = s$prop_delay, n_min_retro_nowcasts = 2),
      error = function(e) NULL)
    if(is.null(tv)) return(NULL)
    fitted <- tryCatch({
      delay_pmf <- baselinenowcast::estimate_delay(rt, n = tv$n_history_delay)
      list(point = baselinenowcast::apply_delay(rt, delay_pmf),
           uncertainty = baselinenowcast::estimate_uncertainty_retro(
             rt, n_history_delay = tv$n_history_delay,
             n_retrospective_nowcasts = tv$n_retrospective_nowcasts))
    }, error = function(e) NULL)
    if(is.null(fitted)) return(NULL)
    set.seed(seed)
    now <- counts %>% dplyr::filter(.data$date_onset %in% scored) %>%
      dplyr::group_by(.data$date_onset) %>%
      dplyr::summarise(now = .data$count[which.max(.data$date_reported)])
    now <- dplyr::left_join(dplyr::tibble(date_onset = scored), now, by = "date_onset") %>%
      dplyr::mutate(now = dplyr::coalesce(.data$now, 0L))
    draws <- lapply(reports$date_report, function(r) {
      # the predicted new counts of the delays which report r has reached
      cut <- as.integer(r - onsets)
      f_aggregate <- function(x) vapply(seq_len(nrow(x)), function(i)
        sum(x[i, seq_len(min(ncol(x), cut[i] + 1))], na.rm = TRUE), numeric(1))
      add <- vapply(seq_len(settings$n_draws), function(k)
        baselinenowcast::sample_prediction(fitted$point, rt, fitted$uncertainty,
                                           delay_aggregator = f_aggregate)[, 1],
        numeric(length(onsets)))
      m <- matrix(0, settings$n_draws, length(scored))
      inside <- match(onsets, scored)
      m[, inside[!is.na(inside)]] <- t(add)[, !is.na(inside), drop = FALSE]
      sweep(m, 2, now$now, "+")
    })
    names(draws) <- as.character(reports$date_report)
    list(model = "baselinenowcast", series = sub("_reported", "", column),
         date_now = date_now, onsets = scored, draws = draws)
  })
  Filter(Negate(is.null), out)
}

#' Projections of the data-driven model at a fit date, from its saved chains
#'
#' @param paths the files of the chains (as fit_nowcast_chain)
#' @param snapshots as build_bvd_snapshots
#' @param date_now the fit date
#' @param settings as comparison_settings
#' @param model_file the Stan file
#' @return a list of two projections (cases and deaths)
project_data_driven <- function(paths, snapshots, date_now, settings = comparison_settings(),
                                model_file = "stan/nowcast_joint.stan") {
  date_now <- as.Date(date_now)
  reports <- comparison_reports(snapshots, date_now, settings)
  scored <- seq(date_now - settings$max_age, date_now, by = "day")
  analysis <- prepare_bvd_reports(snapshots, date_now)
  data <- nowcast_stan_data(analysis)
  fit <- nowcast_join_chains(paths)
  proj <- lapply(seq_len(nrow(reports)), function(i)
    nowcast_project(fit, data, nowcast_day(analysis, reports$date_report[i]), seed = i,
                    model_file = model_file))
  columns <- match(scored, analysis$dates$date)
  lapply(c(cases = "c", deaths = "d"), function(s) {
    draws <- lapply(proj, function(p) comparison_thin(p[[paste0("projected_", s)]][, columns, drop = FALSE],
                                                      settings$n_draws))
    names(draws) <- as.character(reports$date_report)
    list(model = "data-driven", series = if(s == "c") "cases" else "deaths",
         date_now = date_now, onsets = scored, draws = draws)
  })
}

#' Projections of one chain of a mechanistic fit to the later reports scored
#'
#' For each of n_draws of the chain's kept draws (spread evenly over them),
#' incidenceinflation::predict_reports simulates the reports scored from the
#' last report fitted: the split of each onset day's reports is drawn given
#' the draw, and the later reports from it. The chain's saved state supplies
#' what does not change between draws. Onset days not in the fit are left
#' out.
#'
#' @param fit_path the chain's draws (as fit_mechanistic_chain saves them)
#' @param state_path the chain's state file
#' @param snapshots as build_bvd_snapshots
#' @param date_now the fit date
#' @param out the file to save to
#' @param settings as comparison_settings
#' @param n_draws the draws projected
#' @return out, holding a named list of 'chain', 'date_now', 'iterations',
#'   'onsets' and 'draws' (for cases and for deaths, a list with a matrix of
#'   draws by onset day for each report, named by its date)
project_mechanistic_chain <- function(fit_path, state_path, snapshots, date_now, out,
                                      settings = comparison_settings(), n_draws = 100L) {
  date_now <- as.Date(date_now)
  fit <- readRDS(fit_path)
  variables <- readRDS(state_path)$state$variables
  # the delays of the fit's variant of the settings (see mechanistic_variant)
  family <- if(is.null(fit$variant)) NULL else mechanistic_variant(fit$variant)$delay_family
  old <- options(incidenceinflation.delay_family = if(is.null(family)) "gamma" else family)
  on.exit(options(old), add = TRUE)
  if(fit$date_now != date_now)
    stop("The chain is not of the fit date given.")
  date0 <- fit$analysis$dates$date[1]
  f_day <- function(date) as.integer(as.Date(date) - date0) + 1L
  reports <- comparison_reports(snapshots, date_now, settings)
  scored <- seq(date_now - settings$max_age, date_now, by = "day")
  scored <- scored[f_day(scored) %in% variables$deaths_state$onset_times]
  iterations <- sort(unique(fit$deaths$iteration))
  iterations <- iterations[unique(round(seq(1, length(iterations),
                                            length.out = min(n_draws, length(iterations)))))]
  tables <- c("cases", "reporting_days", "revisions", "reporting_effects", "deaths",
              "deaths_parameters")
  set.seed(1000 * fit$chain + as.integer(date_now) %% 1000)
  simulated <- lapply(iterations, function(it) {
    draw <- lapply(fit[tables], function(df) df[df$iteration == it, ])
    incidenceinflation::predict_reports(variables, draw, f_day(reports$date_report),
                                        onset_times = f_day(scored))
  })
  f_series <- function(column) {
    draws <- lapply(reports$date_report, function(date) {
      t(vapply(simulated, function(x) {
        x <- x[x$time_reported == f_day(date), ]
        as.numeric(x[[column]][match(f_day(scored), x$time_onset)])
      }, numeric(length(scored))))
    })
    stats::setNames(draws, as.character(reports$date_report))
  }
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(chain = fit$chain, date_now = date_now, iterations = iterations, onsets = scored,
               draws = list(cases = f_series("cases_reported"),
                            deaths = f_series("deaths_reported"))), out)
  out
}

#' Projections of a mechanistic fit to the later reports scored, from those
#' of its chains, with checks of their convergence
#'
#' @param paths the files of project_mechanistic_chain, one for each chain
#' @param settings as comparison_settings
#' @return a named list of 'projections' (as project_data_driven) and
#'   'diagnostics' (a row for each series: the largest R-hat and smallest
#'   bulk and tail effective sample sizes of the projected counts, over the
#'   onset days and reports whose draws are not all one value, and the share
#'   of them with an R-hat above 1.05)
project_mechanistic <- function(paths, settings = comparison_settings()) {
  chains <- lapply(paths, readRDS)
  date_now <- chains[[1]]$date_now
  onsets <- chains[[1]]$onsets
  n <- min(vapply(chains, function(x) length(x$iterations), numeric(1)))
  out <- lapply(c(cases = "cases", deaths = "deaths"), function(series) {
    reports <- names(chains[[1]]$draws[[series]])
    # draws by chain, the last n of each, for each report
    by_chain <- lapply(reports, function(r) lapply(chains, function(x)
      utils::tail(x$draws[[series]][[r]], n)))
    draws <- lapply(by_chain, function(l) comparison_thin(do.call(rbind, l), settings$n_draws))
    names(draws) <- reports
    checks <- dplyr::bind_rows(lapply(by_chain, function(l) {
      dplyr::bind_rows(lapply(seq_along(onsets), function(j) {
        m <- vapply(l, function(x) x[, j], numeric(n))
        if(length(unique(as.vector(m))) == 1) return(NULL)
        dplyr::tibble(rhat = posterior::rhat(m), ess_bulk = posterior::ess_bulk(m),
                      ess_tail = posterior::ess_tail(m))
      }))
    }))
    list(projection = list(model = "mechanistic", series = series, date_now = date_now,
                           onsets = onsets, draws = draws),
         diagnostics = dplyr::tibble(
           model = "mechanistic", series = series, date_now = date_now,
           chains = length(chains), draws_per_chain = n, cells = nrow(checks),
           max_rhat = max(checks$rhat, na.rm = TRUE),
           share_rhat_over_1.05 = mean(checks$rhat > 1.05, na.rm = TRUE),
           min_ess_bulk = min(checks$ess_bulk, na.rm = TRUE),
           min_ess_tail = min(checks$ess_tail, na.rm = TRUE)))
  })
  list(projections = unname(lapply(out, `[[`, "projection")),
       diagnostics = dplyr::bind_rows(lapply(out, `[[`, "diagnostics")))
}

#' The case fatality ratio of a mechanistic fit from the eventual counts of
#' its draws (comparison_cfr): every case recorded for an onset day is
#' eventually reported and every death recorded, and wrong records are
#' removed, so the eventual counts are the cases and fatal cases recorded for
#' each day
#'
#' @param paths the files of the fit's chains (as fit_mechanistic_chain saves
#'   them)
#' @param min_age as comparison_cfr
mechanistic_cfr <- function(paths, min_age = 21) {
  chains <- lapply(paths, readRDS)
  date_now <- chains[[1]]$date_now
  date0 <- chains[[1]]$analysis$dates$date[1]
  f_wide <- function(df, column) {
    df <- df[order(df$iteration, df$time_onset), ]
    m <- matrix(df[[column]], ncol = length(unique(df$time_onset)), byrow = TRUE)
    colnames(m) <- sort(unique(df$time_onset))
    m
  }
  counts <- lapply(chains, function(x) {
    both <- dplyr::inner_join(x$cases, x$deaths, by = c("iteration", "time_onset"))
    # the fatal cases dated by death out of each day of onset, and into each
    # day they are recorded for
    both$out <- both$deaths_true - both$deaths_recorded + both$dated_by_death
    both$cases_recorded <- both$cases_true - both$out + both$dated_by_death
    list(cases = f_wide(both, "cases_recorded"), deaths = f_wide(both, "deaths_recorded"))
  })
  n <- min(vapply(counts, function(x) nrow(x$cases), numeric(1)))
  cases <- do.call(rbind, lapply(counts, function(x) utils::tail(x$cases, n)))
  deaths <- do.call(rbind, lapply(counts, function(x) utils::tail(x$deaths, n)))
  dates <- date0 + as.integer(colnames(cases)) - 1
  comparison_cfr(cases, deaths, dates, date_now, "mechanistic", min_age)
}

#' The case fatality ratio from the eventual counts of cases and deaths by
#' onset day
#'
#' In each draw, the ratio of the eventual deaths to the eventual cases of
#' each onset week, and of all the onset days at least min_age days before
#' the fit date (pooled). Weeks whose last day is fewer than min_age days
#' before the fit date are marked as recent: most of their deaths are still
#' to come. The draws of the two series are paired by their order; for
#' epinowcast, whose fits to cases and to deaths are separate, they are
#' independent.
#'
#' @param cases,deaths draws by onset day of the eventual counts
#' @param dates the onset days
#' @param date_now the fit date
#' @param model the model's name
#' @param min_age as above
#' @return a tibble of model, date_now, week (NA for the pooled ratio),
#'   recent, and the 5%, 50% and 95% quantiles
comparison_cfr <- function(cases, deaths, dates, date_now, model, min_age = 21) {
  date_now <- as.Date(date_now)
  n <- min(nrow(cases), nrow(deaths))
  cases <- cases[seq_len(n), , drop = FALSE]
  deaths <- deaths[seq_len(n), , drop = FALSE]
  # draws with missing counts (a chain whose predictions failed) are left out
  ok <- stats::complete.cases(cases) & stats::complete.cases(deaths)
  if(!all(ok)) {
    warning(sprintf("%s at %s: %d of %d draws with missing counts left out", model,
                    format(as.Date(date_now)), sum(!ok), n), call. = FALSE)
    cases <- cases[ok, , drop = FALSE]
    deaths <- deaths[ok, , drop = FALSE]
  }
  week <- as.Date(cut(dates, "week"))
  f_week <- function(m) t(rowsum(t(m), as.character(week)))
  ratio <- f_week(deaths) / pmax(f_week(cases), 1)
  weeks <- as.Date(colnames(ratio))
  old <- as.numeric(date_now - dates) >= min_age
  pooled <- rowSums(deaths[, old, drop = FALSE]) / pmax(rowSums(cases[, old, drop = FALSE]), 1)
  q <- function(x) stats::quantile(x, c(0.05, 0.5, 0.95), names = FALSE)
  qs <- rbind(t(apply(ratio, 2, q)), q(pooled))
  recent <- c(as.numeric(date_now - (weeks + 6)) < min_age, FALSE)
  dplyr::tibble(model = model, date_now = date_now, week = c(weeks, as.Date(NA)), recent = recent,
                p05 = qs[, 1], p50 = qs[, 2], p95 = qs[, 3])
}

#' The case fatality ratio of epinowcast's separate fits (comparison_cfr)
#' @param path the file saved by fit_epinowcast
#' @param min_age as comparison_cfr
comparison_cfr_epinowcast <- function(path, min_age = 21) {
  e <- readRDS(path)$eventual
  comparison_cfr(e$cases$draws, e$deaths$draws, e$cases$onsets, e$cases$date_now,
                 "epinowcast", min_age)
}

#' The case fatality ratio of the data-driven model from the eventual counts
#' it projects (comparison_cfr), as nowcast_eventual takes them
#' @inheritParams project_data_driven
#' @param min_age as comparison_cfr
comparison_cfr_data_driven <- function(paths, snapshots, date_now, model_file = "stan/nowcast_joint.stan",
                                       min_age = 21) {
  date_now <- as.Date(date_now)
  analysis <- prepare_bvd_reports(snapshots, date_now)
  data <- nowcast_stan_data(analysis)
  fit <- nowcast_join_chains(paths)
  proj <- nowcast_project(fit, data, nowcast_day(analysis, date_now + 365), model_file = model_file)
  dates <- analysis$dates$date[seq_len(data$stan$T)]
  comparison_cfr(proj$projected_c, proj$projected_d, dates, date_now, "data-driven", min_age)
}

#' Plot of the case fatality ratio by onset week at each fit date, of
#' epinowcast's separate fits and of the data-driven model
#' @param cfr as comparison_cfr, for both models and all fit dates
#' @param out the directory
#' @return the file
plot_comparison_cfr <- function(cfr, out) {
  df <- dplyr::filter(cfr, !is.na(.data$week)) %>%
    dplyr::mutate(fit = format(.data$date_now, "%d %b"),
                  fit = factor(.data$fit, unique(.data$fit[order(.data$date_now)])))
  g <- ggplot2::ggplot(df, ggplot2::aes(.data$week, .data$p50, colour = .data$model,
                                        alpha = ifelse(.data$recent, "recent", "earlier"))) +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = .data$p05, ymax = .data$p95),
                             position = ggplot2::position_dodge(width = 3), size = 0.15) +
    ggplot2::scale_alpha_manual(values = c(earlier = 1, recent = 0.35), name = NULL) +
    ggplot2::facet_wrap(~ fit, ncol = 5) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = "Week of onset", y = "Case fatality ratio (median, 90% interval)", colour = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  comparison_save(file.path(out, "cfr_by_week.png"), g, 11, 5.5)
}

#' The log score of a count under a negative binomial with the mean and
#' variance of the draws (a Poisson if they are not overdispersed)
#' @param y the count
#' @param x the draws
comparison_log_score <- function(y, x) {
  m <- mean(x)
  v <- stats::var(x)
  if(v == 0) return(if(y == m) 0 else -Inf)
  if(v > m) stats::dnbinom(y, size = m^2 / (v - m), mu = m, log = TRUE)
  else stats::dpois(y, m, log = TRUE)
}

#' Scores of a model's projections against the reports
#'
#' For each onset day and report: the CRPS of the count and of log(count + 1),
#' the log score (comparison_log_score), the randomised probability integral
#' transform, and the 50% and 90% intervals (whether they hold the count, and
#' their widths), and whether the draws are all one value (point: a log score
#' of 0 if it is the count). Draws which are not finite or are negative are left out,
#' with a warning, and counted (n_bad, with the smallest of them, min_bad).
#'
#' @param projections as described above
#' @param truth as comparison_truth (of the same fit date)
#' @param seed for the randomised PIT
#' @return a tibble: the rows of truth for the model's series, with the scores
comparison_score <- function(projections, truth, seed = 1) {
  set.seed(seed)
  dplyr::bind_rows(lapply(names(projections$draws), function(r) {
    m <- projections$draws[[r]]
    rows <- truth %>%
      dplyr::filter(.data$series == projections$series, .data$date_report == as.Date(r),
                    .data$date_onset %in% projections$onsets)
    j <- match(rows$date_onset, projections$onsets)
    scores <- lapply(seq_len(nrow(rows)), function(k) {
      x <- m[, j[k]]
      bad <- !is.finite(x) | x < 0
      min_bad <- if(any(bad)) suppressWarnings(min(x[bad], na.rm = TRUE)) else NA_real_
      if(any(bad)) {
        warning(sprintf("%s %s, report %s, onset %s: %d draws not finite or negative (smallest %s), left out",
                        projections$model, projections$series, r, rows$date_onset[k], sum(bad),
                        format(min_bad)), call. = FALSE)
        x <- x[!bad]
        if(length(x) == 0) stop("no finite draws to score")
      }
      y <- rows$shown[k]
      q <- stats::quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95), names = FALSE)
      lower <- mean(x <= y - 1)
      dplyr::tibble(crps = scoringRules::crps_sample(y, x),
                    crps_log = scoringRules::crps_sample(log1p(y), log1p(x)),
                    log_score = comparison_log_score(y, x),
                    pit = lower + stats::runif(1) * (mean(x <= y) - lower),
                    p50 = q[3],
                    in50 = y >= q[2] & y <= q[4], in90 = y >= q[1] & y <= q[5],
                    width50 = q[4] - q[2], width90 = q[5] - q[1], n_bad = sum(bad), min_bad = min_bad,
                    point = all(x == x[1]))
    })
    dplyr::bind_cols(rows, dplyr::bind_rows(scores))
  })) %>% dplyr::mutate(model = projections$model, .before = 1)
}

#' Summaries of the scores, on the onset days without a fall
#'
#' 'common' keeps the onset days which every model scored, so that the
#' models are compared on the same counts; 'own' takes each model on all the
#' days it scored; 'uncertain' keeps those of the common days on which no
#' model's draws are all one value, since a model which projects that a count
#' will not change scores a log score of 0 when it does not, which no model
#' allowing for a change can match. The log score is summed over the onset days of each fit
#' date and report, and the sums averaged over them (log_score_sum), as well
#' as averaged over days (log_score); days given no probability (a log score
#' of -Inf) are counted (no_probability), and log_score_sum_finite leaves them
#' out.
#'
#' @param scores the scores of all the models (comparison_score)
#' @param by the columns to summarise by, besides model and series
#' @return a tibble with a column 'set' ("common", "uncertain" or "own")
comparison_summary <- function(scores, by = "week") {
  scored <- dplyr::filter(scores, .data$no_fall)
  n_models <- dplyr::n_distinct(scored$model)
  common <- scored %>%
    dplyr::group_by(.data$series, .data$date_now, .data$date_report, .data$date_onset) %>%
    dplyr::filter(dplyr::n_distinct(.data$model) == n_models) %>%
    dplyr::ungroup()
  uncertain <- common %>%
    dplyr::group_by(.data$series, .data$date_now, .data$date_report, .data$date_onset) %>%
    dplyr::filter(!any(.data$point)) %>%
    dplyr::ungroup()
  f <- function(x) {
    sums <- x %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("model", "series", by, "date_now", "date_report")))) %>%
      dplyr::summarise(log_score_sum = sum(.data$log_score),
                       log_score_sum_finite = sum(.data$log_score[is.finite(.data$log_score)]),
                       .groups = "drop") %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("model", "series", by)))) %>%
      dplyr::summarise(log_score_sum = mean(.data$log_score_sum),
                       log_score_sum_finite = mean(.data$log_score_sum_finite), .groups = "drop")
    x %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(c("model", "series", by)))) %>%
      dplyr::summarise(n = dplyr::n(), fit_dates = dplyr::n_distinct(.data$date_now),
                       crps = mean(.data$crps), crps_log = mean(.data$crps_log),
                       log_score = mean(.data$log_score),
                       no_probability = sum(!is.finite(.data$log_score)),
                       coverage50 = mean(.data$in50), coverage90 = mean(.data$in90),
                       width50 = mean(.data$width50), width90 = mean(.data$width90),
                       .groups = "drop") %>%
      dplyr::left_join(sums, by = c("model", "series", by))
  }
  dplyr::bind_rows(common = f(common), uncertain = f(uncertain), own = f(scored), .id = "set")
}

#' The age band of an onset day at the fit date
#' @param age days from onset to the fit date
comparison_age_band <- function(age) {
  cut(age, c(-1, 6, 13, 27, Inf), labels = c("0-6", "7-13", "14-27", "28+"))
}

#' Saves a plot
#' @param path the file
#' @param g the plot
#' @param width,height in inches
#' @return path
comparison_save <- function(path, g, width, height) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(path, g, width = width, height = height, dpi = 300)
  path
}

#' Plots of the scores by week after the fit date: the mean CRPS of
#' log(count + 1), and the coverage of the 50% and 90% intervals
#'
#' @param summary as comparison_summary (by week), the 'common' set
#' @param out the folder
#' @return the paths of the plots
plot_comparison_weeks <- function(summary, out) {
  x <- dplyr::filter(summary, .data$set == "common")
  g1 <- ggplot2::ggplot(x, ggplot2::aes(.data$week, .data$crps_log, colour = .data$model)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~ series, scales = "free_y") +
    ggplot2::labs(x = "Weeks after the fit date", y = "Mean CRPS of log(count + 1)", colour = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  coverage <- tidyr::pivot_longer(x, c("coverage50", "coverage90"), names_to = "interval",
                                  values_to = "coverage") %>%
    dplyr::mutate(interval = ifelse(.data$interval == "coverage50", "50% interval", "90% interval"),
                  nominal = ifelse(.data$interval == "50% interval", 0.5, 0.9))
  g2 <- ggplot2::ggplot(coverage, ggplot2::aes(.data$week, .data$coverage, colour = .data$model)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = .data$nominal), linetype = "dashed", colour = "grey40") +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_grid(interval ~ series) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "Weeks after the fit date", y = "Share of counts inside the interval", colour = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  c(comparison_save(file.path(out, "crps_by_week.png"), g1, 9, 4),
    comparison_save(file.path(out, "coverage_by_week.png"), g2, 9, 6))
}

#' Histograms of the randomised PIT of each model and series (onset days
#' without a fall, in the 'common' set)
#'
#' @param scores as comparison_score, all the models
#' @param out the folder
#' @return the path of the plot
plot_comparison_pit <- function(scores, out) {
  x <- dplyr::filter(scores, .data$no_fall)
  n_models <- dplyr::n_distinct(x$model)
  x <- x %>%
    dplyr::group_by(.data$series, .data$date_now, .data$date_report, .data$date_onset) %>%
    dplyr::filter(dplyr::n_distinct(.data$model) == n_models) %>%
    dplyr::ungroup()
  g <- ggplot2::ggplot(x, ggplot2::aes(.data$pit)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, 0.1), fill = "grey60", colour = "white") +
    ggplot2::facet_grid(model ~ series) +
    ggplot2::labs(x = "Probability integral transform", y = "Counts") +
    ggplot2::theme_bw()
  comparison_save(file.path(out, "pit.png"), g, 7, 1 + 1.8 * n_models)
}

#' The mean CRPS of log(count + 1) at each fit date, by model
#'
#' @param scores as comparison_score, all the models
#' @param out the folder
#' @return the path of the plot
plot_comparison_dates <- function(scores, out) {
  x <- scores %>%
    dplyr::filter(.data$no_fall) %>%
    dplyr::group_by(.data$model, .data$series, .data$date_now) %>%
    dplyr::summarise(crps_log = mean(.data$crps_log), .groups = "drop")
  g <- ggplot2::ggplot(x, ggplot2::aes(.data$date_now, .data$crps_log, colour = .data$model)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::facet_wrap(~ series, scales = "free_y") +
    ggplot2::labs(x = "Fit date", y = "Mean CRPS of log(count + 1)", colour = NULL) +
    ggplot2::theme_bw() + ggplot2::theme(legend.position = "bottom")
  comparison_save(file.path(out, "crps_by_fit_date.png"), g, 9, 4)
}
