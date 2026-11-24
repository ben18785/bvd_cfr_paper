# The sensitivity analyses of the mechanistic model set against its main fit
# (share15_tail), all fitted to the reports of the same date: plots of the
# case fatality ratio, Rt, the eventual deaths and the reporting parameters,
# made from the outputs of each fit (the artifact mechanistic-outputs of its
# run of mechanistic.yml; _targets_sensitivity.R)

#' The fits compared and their labels: the main fit first, then each
#' sensitivity analysis (mechanistic_variant), with the change it makes
sensitivity_variants <- function() {
  c(main = "Main fit",
    s1_cfr_loose = "S1: ~5 CFR changes expected (vs 1)",
    s2_misdating5 = "S2: misdated fraction fixed at 0.05",
    s3_backlog_tight = "S3: tighter backlog priors (sd 0.5 vs 1)",
    s4_hold14 = "S4: backlog hold 14 days (vs 21)",
    s4_hold28 = "S4: backlog hold 28 days (vs 21)",
    s10_bundibugyo = "S10: onset to death mean 11.7 d (vs 9.3)",
    s10_short = "S10: onset to death mean 7 d (vs 9.3)",
    s12_fixed_june = "S12: recording delay fixed to 30 Jun (vs 31 Jul)",
    s12_share75 = "S12: share reported at death 0.75 (vs 1)")
}

#' Reads one table of the outputs of every fit
#'
#' @param dir the folder holding a folder of outputs for each fit, named as
#'   sensitivity_variants (those missing are left out, with a warning)
#' @param file the table (by_onset_day.csv, eventual_by_day.csv,
#'   eventual_totals.csv, parameters.csv or convergence.csv)
#' @return a tibble with the fit ('variant') and its label, a factor in the
#'   order of sensitivity_variants
read_sensitivity_table <- function(dir, file) {
  variants <- sensitivity_variants()
  present <- names(variants)[file.exists(file.path(dir, names(variants), file))]
  if(!"main" %in% present) stop("no outputs of the main fit in ", dir)
  if(length(present) < length(variants))
    warning("no outputs of ", paste(setdiff(names(variants), present), collapse = ", "))
  dplyr::bind_rows(lapply(present, function(v)
    utils::read.csv(file.path(dir, v, file)) %>%
      dplyr::mutate(variant = v, label = factor(variants[[v]], levels = variants))))
}

sensitivity_colours <- function() c(main = "grey45", variant = "#D55E00")

sensitivity_theme <- function() {
  ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(hjust = 0, size = 8.5))
}

#' Plots a quantity by onset day for each sensitivity analysis, a panel each,
#' against the main fit
#'
#' @param df rows of the quantity (date, p05, p50, p95, variant, label)
#' @param ylab the label of the y axis
#' @param title the title
#' @return a ggplot
plot_sensitivity_facets <- function(df, ylab, title) {
  cols <- sensitivity_colours()
  df <- dplyr::mutate(df, date = as.Date(.data$date))
  main <- dplyr::select(dplyr::filter(df, .data$variant == "main"), -"variant", -"label")
  other <- dplyr::filter(df, .data$variant != "main")
  ggplot2::ggplot(other, ggplot2::aes(.data$date)) +
    ggplot2::geom_ribbon(data = main, ggplot2::aes(ymin = .data$p05, ymax = .data$p95),
                         fill = cols[["main"]], alpha = 0.3) +
    ggplot2::geom_line(data = main, ggplot2::aes(y = .data$p50), colour = cols[["main"]]) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p05, ymax = .data$p95),
                         fill = cols[["variant"]], alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = .data$p50), colour = cols[["variant"]]) +
    ggplot2::facet_wrap(~label, ncol = 3) +
    ggplot2::scale_x_date(date_labels = "%b", date_breaks = "1 month") +
    ggplot2::labs(x = "Onset date", y = ylab, title = title,
                  subtitle = "Orange: sensitivity analysis; grey: main fit (medians and 90% intervals)") +
    sensitivity_theme()
}

#' Medians and 90% intervals of a quantity, one row for each fit, with the
#' main fit's as a grey band
#'
#' @param df rows of the quantity (p05, p50, p95, variant, label), and the
#'   facet if facet is given
#' @param xlab the label of the x axis
#' @param facet the column of the facets, or NULL
#' @param ref a value marked by a dashed line, or NULL
#' @return a ggplot
plot_sensitivity_forest <- function(df, xlab, facet = NULL, ref = NULL) {
  cols <- sensitivity_colours()
  df <- dplyr::mutate(df, label = factor(.data$label, levels = rev(levels(.data$label))))
  mains <- dplyr::filter(df, .data$variant == "main")
  g <- ggplot2::ggplot(df, ggplot2::aes(y = .data$label, x = .data$p50, xmin = .data$p05,
                                        xmax = .data$p95, colour = .data$variant == "main")) +
    ggplot2::geom_rect(data = mains, ggplot2::aes(xmin = .data$p05, xmax = .data$p95),
                       ymin = -Inf, ymax = Inf, inherit.aes = FALSE, fill = cols[["main"]],
                       alpha = 0.15) +
    ggplot2::geom_vline(data = mains, ggplot2::aes(xintercept = .data$p50),
                        colour = cols[["main"]], linewidth = 0.3)
  if(!is.null(ref)) g <- g + ggplot2::geom_vline(xintercept = ref, linetype = 2, linewidth = 0.3)
  g <- g + ggplot2::geom_pointrange(size = 0.2) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = cols[["variant"]]),
                                 guide = "none") +
    ggplot2::labs(x = xlab, y = NULL) + sensitivity_theme()
  if(!is.null(facet))
    g <- g + ggplot2::facet_wrap(ggplot2::vars(.data[[facet]]), nrow = 1, scales = "free_x")
  g
}

#' The headline estimates: the case fatality ratio on 15 June, the eventual
#' deaths and cases, and Rt for four onset days
#'
#' @param days the rows of by_onset_day.csv of every fit
#' @param totals the rows of eventual_totals.csv of every fit
#' @return a patchwork of ggplots
plot_sensitivity_summary <- function(days, totals) {
  cfr <- dplyr::filter(days, .data$quantity == "Case fatality ratio", .data$date == "2026-06-15")
  rt_dates <- as.Date(c("2026-06-01", "2026-07-01", "2026-08-01", "2026-09-01"))
  rt <- dplyr::filter(days, .data$quantity == "Rt", as.Date(.data$date) %in% rt_dates) %>%
    dplyr::mutate(when = factor(format(as.Date(.data$date), "%d %b"),
                                levels = format(rt_dates, "%d %b")))
  blank_y <- ggplot2::theme(axis.text.y = ggplot2::element_blank())
  top <- patchwork::wrap_plots(
    plot_sensitivity_forest(cfr, "Case fatality ratio (15 June)"),
    plot_sensitivity_forest(dplyr::filter(totals, .data$series == "deaths"), "Eventual deaths") + blank_y,
    plot_sensitivity_forest(dplyr::filter(totals, .data$series == "cases"), "Eventual cases") + blank_y,
    nrow = 1)
  bottom <- plot_sensitivity_forest(rt, expression(R[t] ~ "by onset date"), facet = "when", ref = 1)
  patchwork::wrap_plots(top, bottom, ncol = 1) +
    patchwork::plot_annotation(title = "Sensitivity analyses: headline estimates",
                               subtitle = "Medians and 90% intervals; grey band and line: main fit")
}

#' The parameters of the reports compared across fits: the delays, the heavy
#' tail of the reporting delay, the wrong records, the fraction of fatal
#' cases dated by death and the backlog of deaths
#'
#' @param pars the rows of parameters.csv of every fit
#' @return a ggplot
plot_sensitivity_parameters <- function(pars) {
  keep <- c(
    "reporting_changepoints.centre_mean" = "Case reporting delay: mean, centre (days)",
    "reporting_changepoints.centre_sd" = "Case reporting delay: sd, centre (days)",
    "reporting_effects.tail_weight" = "Case delay: weight of heavy tail",
    "reporting_effects.tail_mean" = "Case delay: mean of heavy tail (days)",
    "revisions.eta" = "Wrong records per genuine record (eta)",
    "revisions.rho" = "Wrong records withdrawn daily (rho)",
    "overdispersion.overdispersion" = "Overdispersion of cases",
    "deaths_parameters.delay_centre_mean" = "Death recording delay: mean, centre (days)",
    "deaths_parameters.delay_centre_sd" = "Death recording delay: sd, centre (days)",
    "deaths_parameters.dated_fraction" = "Fatal cases dated by death (fraction)",
    "deaths_parameters.wrong_dead" = "Wrong records shown dead",
    "deaths_parameters.backlog_share_rw_sd" = "Backlog share: weekly step sd",
    "deaths_parameters.backlog_rate_rw_sd" = "Backlog rate: weekly step sd")
  df <- pars %>%
    dplyr::mutate(key = paste(.data$part, .data$parameter, sep = ".")) %>%
    dplyr::filter(.data$key %in% names(keep)) %>%
    dplyr::mutate(name = factor(keep[.data$key], levels = keep),
                  p05 = .data$q5, p50 = .data$median, p95 = .data$q95)
  plot_sensitivity_forest(df, NULL) +
    ggplot2::facet_wrap(~name, scales = "free_x", ncol = 4) +
    ggplot2::labs(title = "Sensitivity analyses: reporting parameters",
                  subtitle = "Medians and 90% intervals; grey band and line: main fit") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7),
                   axis.text.x = ggplot2::element_text(size = 7))
}

#' Writes the plots of the sensitivity analyses
#'
#' @param dir the folder of the fits' outputs (read_sensitivity_table)
#' @param out the folder of the plots
#' @return the paths of the files written
write_sensitivity_plots <- function(dir, out) {
  days <- read_sensitivity_table(dir, "by_onset_day.csv")
  by_day <- read_sensitivity_table(dir, "eventual_by_day.csv")
  totals <- read_sensitivity_table(dir, "eventual_totals.csv")
  pars <- read_sensitivity_table(dir, "parameters.csv")
  conv <- read_sensitivity_table(dir, "convergence.csv")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  f_save <- function(g, name, width = 9, height = 7.5) {
    path <- file.path(out, name)
    ggplot2::ggsave(path, g, width = width, height = height, dpi = 150)
    path
  }
  f_days <- function(q) dplyr::filter(days, .data$quantity == q)
  status <- file.path(out, "sensitivity_convergence.csv")
  utils::write.csv(dplyr::select(conv, "variant", "date_now", "iterations", "rhat_max", "worst"),
                   status, row.names = FALSE)
  c(status,
    f_save(plot_sensitivity_summary(days, totals), "sensitivity_summary.png", 11, 7.5),
    f_save(plot_sensitivity_parameters(pars), "sensitivity_reporting_parameters.png", 13, 8.5),
    f_save(plot_sensitivity_facets(f_days("Case fatality ratio"), "Case fatality ratio",
                                   "Case fatality ratio by onset date"),
           "sensitivity_cfr_by_onset.png"),
    f_save(plot_sensitivity_facets(f_days("Rt"), expression(R[t]), "Reproduction number") +
             ggplot2::geom_hline(yintercept = 1, linetype = 2, linewidth = 0.3) +
             ggplot2::scale_y_log10(breaks = c(0.5, 0.75, 1, 1.5, 2, 4, 8)),
           "sensitivity_rt_by_onset.png"),
    f_save(plot_sensitivity_facets(dplyr::filter(by_day, .data$series == "deaths"),
                                   "Eventual deaths", "Eventual deaths by onset date"),
           "sensitivity_eventual_deaths_by_onset.png"),
    f_save(plot_sensitivity_facets(f_days("Mean delay to a case being reported (days)"), "Days",
                                   "Mean delay to a case being reported"),
           "sensitivity_case_delay_by_onset.png"),
    f_save(plot_sensitivity_facets(f_days("Mean delay to a death being recorded (days)"), "Days",
                                   "Mean delay to a death being recorded"),
           "sensitivity_death_delay_by_onset.png"))
}
