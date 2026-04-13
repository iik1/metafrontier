utils::globalVariables(".data")

#' Autoplot Method for Metafrontier Objects
#'
#' Creates diagnostic and summary plots for metafrontier objects
#' using \pkg{ggplot2}.
#'
#' @param object a \code{"metafrontier"} object.
#' @param which character. Which plot to produce:
#'   \code{"tgr"} (default) for TGR density distributions by group
#'   (degenerate groups with zero-variance TGR are annotated and excluded
#'   from the density),
#'   \code{"efficiency"} for the TE/TGR/TE* decomposition (boxplots),
#'   \code{"decomposition"} for grouped bars of mean TE, TGR, and TE*,
#'   \code{"frontier"} for metafrontier vs group frontiers scatter plot.
#' @param ... additional arguments (currently unused).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
#'   fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#'                       group = "group")
#'   ggplot2::autoplot(fit, which = "tgr")
#'   ggplot2::autoplot(fit, which = "decomposition")
#' }
#' }
#'
#' @export
autoplot.metafrontier <- function(object,
                                  which = c("tgr", "efficiency",
                                            "decomposition", "frontier"),
                                  ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot. ",
         "Install it with install.packages('ggplot2').", call. = FALSE)
  }

  which <- match.arg(which)

  df <- data.frame(
    group = object$group_vec,
    te_group = object$te_group,
    tgr = object$tgr,
    te_meta = object$te_meta,
    stringsAsFactors = FALSE
  )

  if (which == "tgr") {
    # Filter out degenerate groups (near-zero variance, e.g. TGR = 1 exactly)
    group_sds <- tapply(df$tgr, df$group, stats::sd, na.rm = TRUE)
    degenerate <- names(group_sds)[group_sds < 1e-8]
    df_plot <- df[!df$group %in% degenerate, , drop = FALSE]

    p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$tgr,
                                                fill = .data$group)) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::labs(x = "Technology Gap Ratio",
                    y = "Density",
                    title = "TGR Distribution by Group",
                    fill = "Group") +
      ggplot2::theme_minimal()

    if (length(degenerate) > 0) {
      label <- paste0(degenerate, ": TGR = ",
                      round(group_sds[degenerate] + 1, 4),
                      collapse = "; ")
      p <- p + ggplot2::annotate("text", x = 1, y = Inf, label = label,
                                  hjust = 1, vjust = 1.5, size = 3,
                                  colour = "grey40")
    }

  } else if (which == "efficiency") {
    df_long <- data.frame(
      group = rep(df$group, 3),
      type = rep(c("TE (group)", "TGR", "TE*"), each = nrow(df)),
      value = c(df$te_group, df$tgr, df$te_meta),
      stringsAsFactors = FALSE
    )
    df_long$type <- factor(df_long$type,
                           levels = c("TE (group)", "TGR", "TE*"))

    p <- ggplot2::ggplot(df_long,
                         ggplot2::aes(x = .data$group,
                                      y = .data$value,
                                      fill = .data$type)) +
      ggplot2::geom_boxplot(position = "dodge") +
      ggplot2::labs(x = "Group", y = "Value",
                    title = "Efficiency Decomposition",
                    fill = "Component") +
      ggplot2::theme_minimal()

  } else if (which == "decomposition") {
    # Grouped bars showing TE, TGR, TE* (multiplicative decomposition)
    df_summary <- do.call(rbind, lapply(unique(df$group), function(g) {
      idx <- df$group == g
      data.frame(
        group = g,
        component = c("TE", "TGR", "TE*"),
        value = c(mean(df$te_group[idx]),
                  mean(df$tgr[idx]),
                  mean(df$te_meta[idx])),
        stringsAsFactors = FALSE
      )
    }))
    df_summary$component <- factor(df_summary$component,
                                    levels = c("TE", "TGR", "TE*"))

    p <- ggplot2::ggplot(df_summary,
                         ggplot2::aes(x = .data$group,
                                      y = .data$value,
                                      fill = .data$component)) +
      ggplot2::geom_col(position = ggplot2::position_dodge(0.8),
                        width = 0.7, alpha = 0.85) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed",
                          colour = "grey40") +
      ggplot2::labs(x = "Group", y = "Mean Efficiency",
                    title = "Efficiency Decomposition (TE* = TE x TGR)",
                    fill = "Component") +
      ggplot2::coord_cartesian(ylim = c(0, 1.05)) +
      ggplot2::theme_minimal()

  } else if (which == "frontier") {
    if (is.null(object$meta_frontier) || is.null(object$group_frontier)) {
      stop("Frontier plot requires SFA-based metafrontier.", call. = FALSE)
    }
    df$meta_frontier <- object$meta_frontier
    df$group_frontier <- object$group_frontier

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$group_frontier,
                                           y = .data$meta_frontier,
                                           color = .data$group)) +
      ggplot2::geom_point(alpha = 0.6) +
      ggplot2::geom_abline(slope = 1, intercept = 0,
                           linetype = "dashed", color = "grey40") +
      ggplot2::labs(x = "Group Frontier",
                    y = "Metafrontier",
                    title = "Metafrontier vs Group Frontiers",
                    color = "Group") +
      ggplot2::theme_minimal()
  }

  p
}


#' Autoplot Method for Malmquist Meta Objects
#'
#' @param object a \code{"malmquist_meta"} object.
#' @param which character. Which plot:
#'   \code{"decomposition"} (default), \code{"tgr_evolution"},
#'   or \code{"mpi_trend"}.
#' @param ... additional arguments.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   panels <- lapply(1:3, function(t) {
#'     sim <- simulate_metafrontier(n_groups = 2, n_per_group = 30,
#'                                  seed = 42 + t)
#'     sim$data$time <- t
#'     sim$data$id <- seq_len(nrow(sim$data))
#'     sim$data
#'   })
#'   pdata <- do.call(rbind, panels)
#'   malm <- malmquist_meta(log_y ~ log_x1 + log_x2, data = pdata,
#'                          group = "group", time = "time")
#'   ggplot2::autoplot(malm, which = "decomposition")
#' }
#' }
#'
#' @export
autoplot.malmquist_meta <- function(object,
                                    which = c("decomposition",
                                              "tgr_evolution",
                                              "mpi_trend"),
                                    ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.", call. = FALSE)
  }

  which <- match.arg(which)
  m <- object$malmquist

  if (which == "decomposition") {
    df_mean <- do.call(rbind, lapply(unique(m$group), function(g) {
      idx <- m$group == g
      data.frame(
        group = g,
        component = c("TEC", "TGC", "TC*"),
        value = c(mean(m$TEC[idx], na.rm = TRUE),
                  mean(m$TGC[idx], na.rm = TRUE),
                  mean(m$TC[idx], na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }))

    p <- ggplot2::ggplot(df_mean,
                         ggplot2::aes(x = .data$group,
                                      y = .data$value,
                                      fill = .data$component)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
      ggplot2::labs(x = "Group", y = "Mean Value",
                    title = "Malmquist Decomposition",
                    fill = "Component") +
      ggplot2::theme_minimal()

  } else if (which == "tgr_evolution") {
    tgr_df <- object$tgr
    tgr_long <- data.frame(
      group = rep(tgr_df$group, 2),
      period = c(tgr_df$period_from, tgr_df$period_to),
      TGR = c(tgr_df$TGR_from, tgr_df$TGR_to),
      stringsAsFactors = FALSE
    )
    tgr_mean <- aggregate(TGR ~ group + period, data = tgr_long, mean,
                          na.rm = TRUE)

    p <- ggplot2::ggplot(tgr_mean,
                         ggplot2::aes(x = .data$period,
                                      y = .data$TGR,
                                      color = .data$group)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::labs(x = "Period", y = "Mean TGR",
                    title = "TGR Evolution",
                    color = "Group") +
      ggplot2::theme_minimal()

  } else {
    mpi_mean <- aggregate(MPI ~ group + period_from, data = m, mean,
                          na.rm = TRUE)

    p <- ggplot2::ggplot(mpi_mean,
                         ggplot2::aes(x = .data$period_from,
                                      y = .data$MPI,
                                      color = .data$group)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed") +
      ggplot2::labs(x = "Period", y = "Mean MPI",
                    title = "MPI Trend",
                    color = "Group") +
      ggplot2::theme_minimal()
  }

  p
}


#' Autoplot Method for Bootstrap TGR Objects
#'
#' @param object a \code{"boot_tgr"} object.
#' @param which character. \code{"distribution"} (default) or
#'   \code{"ci"}.
#' @param ... additional arguments.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   sim <- simulate_metafrontier(n_groups = 2, n_per_group = 50, seed = 42)
#'   fit <- metafrontier(log_y ~ log_x1 + log_x2, data = sim$data,
#'                       group = "group", meta_type = "stochastic")
#'   boot <- boot_tgr(fit, R = 50, seed = 1, progress = FALSE)
#'   ggplot2::autoplot(boot)
#' }
#' }
#'
#' @export
autoplot.boot_tgr <- function(object,
                              which = c("distribution", "ci"),
                              ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for autoplot.", call. = FALSE)
  }

  which <- match.arg(which)

  if (which == "distribution") {
    # Distribution of mean TGR per group across bootstrap reps
    frames <- list()
    for (g in object$groups) {
      idx <- which(object$group_vec == g)
      group_means <- apply(object$tgr_boot[, idx, drop = FALSE], 1, mean)
      frames[[g]] <- data.frame(
        group = g,
        mean_tgr = group_means,
        stringsAsFactors = FALSE
      )
    }
    df <- do.call(rbind, frames)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mean_tgr,
                                           fill = .data$group)) +
      ggplot2::geom_histogram(bins = 30, alpha = 0.6,
                              position = "identity") +
      ggplot2::labs(x = "Bootstrap Mean TGR",
                    y = "Count",
                    title = "Bootstrap TGR Distribution",
                    fill = "Group") +
      ggplot2::theme_minimal()

  } else {
    ci <- object$ci_group
    p <- ggplot2::ggplot(ci,
                         ggplot2::aes(x = .data$Group,
                                      y = .data$Mean_TGR)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = ci[, 3], ymax = ci[, 4]),
        width = 0.2
      ) +
      ggplot2::labs(x = "Group", y = "Mean TGR",
                    title = "Bootstrap CI for Mean TGR") +
      ggplot2::theme_minimal()
  }

  p
}
