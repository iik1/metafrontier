#' Plot a Metafrontier Object
#'
#' Produces diagnostic and summary plots for a fitted metafrontier
#' model.
#'
#' @param x a \code{"metafrontier"} object.
#' @param which character. Type of plot to produce:
#'   \code{"tgr"} (default) for TGR distributions by group,
#'   \code{"efficiency"} for TE* vs TE scatter coloured by group,
#'   \code{"decomposition"} for side-by-side boxplots of TE, TGR,
#'   and TE* by group, or \code{"frontier"} for group frontier
#'   vs metafrontier (only for single-input SFA models).
#' @param ... additional graphical parameters passed to base
#'   plotting functions.
#'
#' @return Invisibly returns \code{x}.
#'
#' @examples
#' set.seed(42)
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2,
#'                     data = sim$data, group = "group")
#' plot(fit, which = "tgr")
#' plot(fit, which = "decomposition")
#'
#' @export
plot.metafrontier <- function(x,
                              which = c("tgr", "efficiency",
                                        "decomposition", "frontier"),
                              ...) {
  which <- match.arg(which)

  groups <- x$groups
  group_vec <- x$group_vec
  n_groups <- length(groups)

  # Colour palette
  cols <- if (n_groups <= 8) {
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
      "#FF7F00", "#FFFF33", "#A65628", "#F781BF")[seq_len(n_groups)]
  } else {
    grDevices::rainbow(n_groups)
  }

  if (which == "tgr") {
    .plot_tgr(x, groups, group_vec, cols, ...)

  } else if (which == "efficiency") {
    .plot_efficiency(x, groups, group_vec, cols, ...)

  } else if (which == "decomposition") {
    .plot_decomposition(x, groups, group_vec, cols, ...)

  } else if (which == "frontier") {
    .plot_frontier(x, groups, group_vec, cols, ...)
  }

  invisible(x)
}


#' @noRd
.plot_tgr <- function(x, groups, group_vec, cols, ...) {
  tgr_list <- technology_gap_ratio(x, by_group = TRUE)

  xlim <- range(unlist(tgr_list), na.rm = TRUE)
  xlim[1] <- max(0, xlim[1] - 0.05)
  xlim[2] <- min(1, xlim[2] + 0.05)

  dens_list <- lapply(tgr_list, density, from = xlim[1], to = xlim[2])
  ylim <- c(0, max(sapply(dens_list, function(d) max(d$y))) * 1.1)

  plot(NULL, xlim = xlim, ylim = ylim,
       xlab = "Technology Gap Ratio (TGR)",
       ylab = "Density",
       main = "Technology Gap Ratio by Group",
       ...)

  for (i in seq_along(groups)) {
    lines(dens_list[[i]], col = cols[i], lwd = 2)
  }

  legend("topleft", legend = groups, col = cols, lwd = 2, bty = "n")
}


#' @noRd
.plot_efficiency <- function(x, groups, group_vec, cols, ...) {
  plot(x$te_group, x$te_meta,
       col = cols[as.integer(group_vec)],
       pch = 16, cex = 0.8,
       xlab = "Group Efficiency (TE)",
       ylab = "Metafrontier Efficiency (TE*)",
       main = "TE* vs TE by Group",
       xlim = c(0, 1), ylim = c(0, 1),
       ...)
  abline(0, 1, lty = 2, col = "grey50")
  legend("topleft", legend = groups, col = cols, pch = 16, bty = "n")
}


#' @noRd
.plot_decomposition <- function(x, groups, group_vec, cols, ...) {
  old_par <- par(mfrow = c(1, 3), mar = c(5, 4, 3, 1))
  on.exit(par(old_par))

  boxplot(x$te_group ~ group_vec, col = cols,
          main = "Group Efficiency (TE)",
          ylab = "TE", xlab = "Group",
          ylim = c(0, 1), ...)

  boxplot(x$tgr ~ group_vec, col = cols,
          main = "Technology Gap Ratio (TGR)",
          ylab = "TGR", xlab = "Group",
          ylim = c(0, 1), ...)

  boxplot(x$te_meta ~ group_vec, col = cols,
          main = "Metafrontier Efficiency (TE*)",
          ylab = "TE*", xlab = "Group",
          ylim = c(0, 1), ...)
}


#' @noRd
.plot_frontier <- function(x, groups, group_vec, cols, ...) {
  if (is.null(x$meta_coef) || length(x$meta_coef) > 2) {
    message("Frontier plot is only available for single-input SFA models.")
    return(invisible(NULL))
  }

  X_col <- if (ncol(x$data) > 0) {
    x_name <- names(x$meta_coef)[2]
    x$data[[x_name]]
  } else {
    return(invisible(NULL))
  }

  y <- model.response(model.frame(x$formula, data = x$data))

  plot(X_col, y,
       col = cols[as.integer(group_vec)],
       pch = 16, cex = 0.6,
       xlab = names(x$meta_coef)[2],
       ylab = "log(y)",
       main = "Group Frontiers and Metafrontier",
       ...)

  # Plot group frontiers
  x_seq <- seq(min(X_col), max(X_col), length.out = 100)
  for (i in seq_along(groups)) {
    beta_g <- x$group_coef[[groups[i]]]
    y_g <- beta_g[1] + beta_g[2] * x_seq
    lines(x_seq, y_g, col = cols[i], lwd = 2, lty = 2)
  }

  # Plot metafrontier
  y_meta <- x$meta_coef[1] + x$meta_coef[2] * x_seq
  lines(x_seq, y_meta, col = "black", lwd = 3)

  legend("topleft",
         legend = c(groups, "Metafrontier"),
         col = c(cols, "black"),
         lwd = c(rep(2, length(groups)), 3),
         lty = c(rep(2, length(groups)), 1),
         bty = "n")
}
