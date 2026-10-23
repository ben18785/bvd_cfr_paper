# The situation reports of the Bundibugyo virus disease outbreak, 2026: the
# confirmed cases and deaths by date of symptom onset shown in each report.
#
# The raw files (data/raw/, restored by scripts/decrypt_bvd_data.sh) are not
# tracked: each is the epidemic curve of one situation report, extracted from
# its figure and verified by hand.

#' Reads one verified extraction of a situation report's epidemic curve
#'
#' Columns are read by position, since one file lacks names for the verified
#' columns: the fifth to seventh columns hold the verified total, alive and
#' dead counts.
#'
#' @param path path to a Varified_SitRep_*.csv file
#' @return a tibble with columns sitrep, date_reported, date_onset,
#'   cases_reported (all confirmed cases) and deaths_reported (those who have died)
read_bvd_sitrep <- function(path) {
  name <- basename(path)
  parts <- regmatches(name, regexec(
    "_([0-9]{3})_([0-9]{2})[_-]([0-9]{2})[_-]([0-9]{4})_", name))[[1]]
  df <- utils::read.csv(path, header = FALSE, skip = 1, stringsAsFactors = FALSE)
  dplyr::tibble(
    sitrep = as.integer(parts[2]),
    date_reported = as.Date(paste(parts[5], parts[4], parts[3], sep = "-")),
    date_onset = as.Date(df[[1]], format = "%d/%m/%Y"),
    cases_reported = as.integer(df[[5]]),
    deaths_reported = as.integer(df[[7]])) %>%
    dplyr::filter(!is.na(.data$date_onset), !is.na(.data$cases_reported),
                  .data$date_onset <= .data$date_reported)
}

#' The counts of every situation report
#'
#' The same situation report can be present in more than one file: one copy
#' is kept, having checked that they agree.
#'
#' @param paths the Varified_SitRep_*.csv files
#' @return a tibble as read_bvd_sitrep, for all reports
build_bvd_snapshots <- function(paths) {
  snapshots <- dplyr::bind_rows(lapply(sort(paths), read_bvd_sitrep))
  conflicting <- snapshots %>%
    dplyr::group_by(.data$date_onset, .data$date_reported) %>%
    dplyr::summarise(n = dplyr::n_distinct(.data$cases_reported), .groups = "drop") %>%
    dplyr::filter(.data$n > 1)
  if(nrow(conflicting) > 0)
    stop("Files for the same situation report give different case counts.")
  snapshots %>%
    dplyr::distinct(.data$date_onset, .data$date_reported, .keep_all = TRUE) %>%
    dplyr::arrange(.data$date_onset, .data$date_reported)
}

#' The reports up to a date, as onset days and report days counted from the
#' first onset day
#'
#' Counts reported on the onset day itself are taken to be genuine, so
#' cannot later fall below this: where they do, the onset-day report is
#' dropped and the onset day is treated as first seen in its next report.
#'
#' The report of date_corrected corrected the records of the onset days
#' before date_onset_corrected (it removed many records, added some and
#' changed which were shown as dead). For those onset days, the reports made
#' before it are left out.
#'
#' @param snapshots as build_bvd_snapshots
#' @param date_now the last report used
#' @param date_first the first onset day
#' @param date_corrected,date_onset_corrected the correcting report, and the
#'   first onset day it did not correct
#' @return a named list: 'snapshots' (time_onset, time_reported,
#'   cases_reported, deaths_reported), 'dates' (the date of each onset day)
#'   and date_now
prepare_bvd_reports <- function(snapshots, date_now, date_first = "2026-05-11",
                                date_corrected = "2026-09-06",
                                date_onset_corrected = "2026-06-28") {
  date_now <- as.Date(date_now)
  date_start <- as.Date(date_first)
  f_time <- function(date) as.integer(as.Date(date) - date_start) + 1L
  known <- snapshots %>%
    dplyr::filter(.data$date_onset >= date_start, .data$date_onset <= date_now,
                  .data$date_reported <= date_now)
  falls_below <- known %>%
    dplyr::group_by(.data$date_onset) %>%
    dplyr::arrange(.data$date_reported, .by_group = TRUE) %>%
    dplyr::summarise(
      drop_first = .data$date_reported[1] == .data$date_onset[1] &
        any(.data$cases_reported < .data$cases_reported[1]))
  reports <- known %>%
    dplyr::left_join(falls_below, by = "date_onset") %>%
    dplyr::filter(!(.data$drop_first & .data$date_reported == .data$date_onset)) %>%
    dplyr::transmute(time_onset = f_time(.data$date_onset),
                     time_reported = f_time(.data$date_reported),
                     cases_reported = as.integer(.data$cases_reported),
                     deaths_reported = as.integer(.data$deaths_reported)) %>%
    dplyr::arrange(.data$time_onset, .data$time_reported)
  if(any(is.na(reports$deaths_reported)))
    stop("Some reports of cases have no count of deaths.")
  if(any(reports$deaths_reported > reports$cases_reported))
    stop("Some reports show more deaths than cases.")
  if(as.Date(date_corrected) <= date_now & as.Date(date_onset_corrected) > date_start)
    reports <- dplyr::filter(reports, !(.data$time_onset < f_time(date_onset_corrected) &
                                          .data$time_reported < f_time(date_corrected)))
  list(snapshots = reports,
       dates = dplyr::tibble(date = seq(date_start, date_now, by = 1),
                             time_onset = seq_along(seq(date_start, date_now, by = 1))),
       date_now = date_now)
}
