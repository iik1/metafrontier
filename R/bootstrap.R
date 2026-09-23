#' Bootstrap Confidence Intervals for the Technology Gap Ratio
#'
#' Computes bootstrap confidence intervals for TGR estimates from
#' a fitted metafrontier model. Supports both parametric (residual
#' resampling) and nonparametric (case resampling) bootstraps.
#'
#' The parametric bootstrap keeps the design fixed and redraws the
#' response, so column \code{i} of \code{tgr_boot} is the TGR of DMU
#' \code{i} in every replicate. It returns observation-level intervals
#' (\code{ci}) as well as group-level intervals.
#'
#' The nonparametric bootstrap resamples rows with replacement within
#' each group, stacks the resampled rows group by group in the order of
#' \code{groups}, and refits the model. Column \code{j} of
#' \code{tgr_boot} is then the \code{j}-th resampled unit, a random
#' draw from group \code{boot_group[j]}, not DMU \code{j}. The
#' group-level intervals are computed over these position blocks and
#' do not depend on the row order of the data. No observation-level
#' intervals are returned (\code{ci} is \code{NULL}): a column is not
#' tied to any DMU, and evaluating the original DMUs against each
#' bootstrap frontier instead is, for DEA, the naive bootstrap, which
#' is inconsistent (Kneip, Simar and Wilson, 2008). For SFA fits,
#' observation-level intervals are available from the parametric
#' bootstrap.
#'
#' Group-level intervals are percentile intervals of the within-group
#' mean (\code{ci_group}) and median (\code{ci_group_median}) of the
#' TGR across replicates. Rows dropped from the fit because of missing
#' values are left out of both bootstraps.
#'
#' @param object a \code{"metafrontier"} object.
#' @param R integer. Number of bootstrap replications (default 999).
#' @param type character. \code{"parametric"} resamples from estimated
#'   error distributions; \code{"nonparametric"} resamples rows within
#'   groups with replacement.
#' @param level numeric. Confidence level (default 0.95).
#' @param ci_type character. \code{"percentile"} (default) or
#'   \code{"bca"} (bias-corrected and accelerated). Applies to the
#'   observation-level intervals of the parametric bootstrap; the
#'   group-level intervals are always percentile intervals.
#' @param seed optional integer seed for reproducibility.
#' @param progress logical. Show progress bar (default \code{TRUE}).
#' @param ncores integer. Number of CPU cores for parallel bootstrap
#'   (default 1, sequential). Requires the \code{parallel} package.
#' @param ... additional arguments passed to \code{\link{metafrontier}}.
#'
#' @return An object of class \code{"boot_tgr"} containing:
#'   \describe{
#'     \item{tgr_boot}{R x n matrix of bootstrapped TGR values. For the
#'       parametric bootstrap, column \code{i} is DMU \code{i}; for the
#'       nonparametric bootstrap, column \code{j} is the \code{j}-th
#'       resampled unit (see Details).}
#'     \item{tgr_original}{original TGR estimates}
#'     \item{ci}{n x 2 matrix of observation-level confidence intervals
#'       (parametric bootstrap), or \code{NULL} (nonparametric
#'       bootstrap)}
#'     \item{ci_group}{data frame of group-level mean TGR intervals}
#'     \item{ci_group_median}{data frame of group-level median TGR
#'       intervals}
#'     \item{R_effective}{number of successful replications}
#'     \item{R}{requested number of replications}
#'     \item{type}{bootstrap type used}
#'     \item{ci_type}{CI type used}
#'     \item{level}{confidence level}
#'     \item{group_vec}{group of each element of \code{tgr_original}}
#'     \item{boot_group}{group of each column of \code{tgr_boot}}
#'     \item{groups}{group labels}
#'   }
#'
#' @references
#' Kneip, A., Simar, L. and Wilson, P.W. (2008). Asymptotics and
#' consistent bootstraps for DEA estimators in nonparametric frontier
#' models. \emph{Econometric Theory}, 24(6), 1663--1697.
#' \doi{10.1017/S0266466608080651}
#'
#' @examples
#' \donttest{
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 100,
#'                              seed = 42)
#' fit <- metafrontier(log_y ~ log_x1 + log_x2,
#'                     data = sim$data, group = "group",
#'                     meta_type = "stochastic")
#' boot <- boot_tgr(fit, R = 50, seed = 1)
#' print(boot)
#' confint(boot)
#'
#' # Nonparametric bootstrap: group-level intervals only
#' boot_np <- boot_tgr(fit, R = 50, type = "nonparametric", seed = 1)
#' boot_np$ci_group
#' boot_np$ci_group_median
#' }
#'
#' @export
boot_tgr <- function(object, R = 999,
                     type = c("parametric", "nonparametric"),
                     level = 0.95,
                     ci_type = c("percentile", "bca"),
                     seed = NULL,
                     progress = TRUE,
                     ncores = 1L,
                     ...) {

  type <- match.arg(type)
  ci_type <- match.arg(ci_type)

  if (!inherits(object, "metafrontier")) {
    stop("'object' must be a fitted metafrontier model.", call. = FALSE)
  }

  # The parametric bootstrap draws new noise and inefficiency terms
  # from the estimated error distributions, which requires a
  # distributional model. DEA is nonparametric and provides no such
  # model, so only case resampling is valid for DEA fits.
  if (type == "parametric" && object$method == "dea") {
    stop("Parametric bootstrap is not available for DEA metafrontiers. ",
         "Use type = 'nonparametric'.", call. = FALSE)
  }

  if (!is.numeric(R) || length(R) != 1L || R < 1L) {
    stop("'R' must be a positive integer (number of bootstrap replications).",
         call. = FALSE)
  }
  R <- as.integer(R)

  if (!is.null(seed)) set.seed(seed)

  n <- length(object$tgr)
  tgr_boot <- matrix(NA_real_, nrow = R, ncol = n)
  n_fail <- 0L

  if (ncores > 1L) {
    # Parallel bootstrap
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("parallel package not available. Falling back to sequential.",
              call. = FALSE)
      ncores <- 1L
    }
  }

  if (ncores > 1L) {
    if (progress) message("Running ", R, " bootstrap replicates on ", ncores, " cores...")

    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Export required functions and objects to workers
    parallel::clusterExport(cl, c("object", "type"), envir = environment())
    parallel::clusterEvalQ(cl, library(metafrontier))

    # Set seeds for reproducibility
    if (!is.null(seed)) {
      parallel::clusterSetRNGStream(cl, seed)
    }

    boot_results <- parallel::parLapply(cl, seq_len(R), function(b) {
      tryCatch(
        .boot_one_replicate(object, type, ...),
        error = function(e) NULL
      )
    })

    for (b in seq_len(R)) {
      if (!is.null(boot_results[[b]])) {
        tgr_boot[b, ] <- boot_results[[b]]
      } else {
        n_fail <- n_fail + 1L
      }
    }
  } else {
    # Sequential bootstrap
    if (progress) {
      pb <- utils::txtProgressBar(min = 0, max = R, style = 3)
    }

    for (b in seq_len(R)) {
      boot_result <- tryCatch(
        .boot_one_replicate(object, type, ...),
        error = function(e) NULL
      )

      if (!is.null(boot_result)) {
        tgr_boot[b, ] <- boot_result
      } else {
        n_fail <- n_fail + 1L
      }

      if (progress) utils::setTxtProgressBar(pb, b)
    }

    if (progress) close(pb)
  }

  R_effective <- R - n_fail
  if (n_fail > 0) {
    fail_pct <- round(100 * n_fail / R, 1)
    if (fail_pct > 10) {
      warning(n_fail, " of ", R, " bootstrap replications failed (",
              fail_pct, "%). Results may be unreliable.", call. = FALSE)
    } else {
      message(n_fail, " of ", R, " replications failed; ",
              R_effective, " successful.")
    }
  }

  # Remove failed rows
  tgr_boot <- tgr_boot[!is.na(tgr_boot[, 1]), , drop = FALSE]

  # Group of each original TGR, and of each column of tgr_boot. The
  # parametric bootstrap keeps the design, so column i is DMU i. The
  # nonparametric bootstrap stacks the resampled rows group by group,
  # so column j is a draw from the j-th position block. Rows dropped by
  # na.action in an SFA fit have no TGR.
  tgr_group <- if (!is.null(object$valid_rows)) {
    object$group_vec[object$valid_rows]
  } else {
    object$group_vec
  }
  boot_group <- if (type == "nonparametric") {
    rep(object$groups, times = as.vector(table(tgr_group)))
  } else {
    as.character(tgr_group)
  }

  # Observation-level CIs are only defined when columns are DMUs
  alpha <- (1 - level) / 2
  ci <- if (type == "parametric") {
    .boot_ci(tgr_boot, object$tgr, alpha, ci_type)
  } else {
    NULL
  }

  # Group-level mean and median TGR CIs
  ci_group <- .boot_ci_group(tgr_boot, object$tgr, tgr_group, boot_group,
                             object$groups, alpha, stat = "mean")
  ci_group_median <- .boot_ci_group(tgr_boot, object$tgr, tgr_group,
                                    boot_group, object$groups, alpha,
                                    stat = "median")

  out <- list(
    tgr_boot = tgr_boot,
    tgr_original = object$tgr,
    ci = ci,
    ci_group = ci_group,
    ci_group_median = ci_group_median,
    R_effective = R_effective,
    R = R,
    type = type,
    ci_type = ci_type,
    level = level,
    group_vec = tgr_group,
    boot_group = boot_group,
    groups = object$groups
  )
  class(out) <- "boot_tgr"
  out
}


# ---------- Internal: single bootstrap replicate ----------

.boot_one_replicate <- function(object, type, ...) {
  formula <- object$formula
  groups <- object$groups

  # Only rows used in the original fit (an SFA fit drops rows with
  # missing values) enter the bootstrap, so the data line up with the
  # group models and every replicate returns one TGR per fitted row
  rows <- if (!is.null(object$valid_rows)) {
    object$valid_rows
  } else {
    seq_len(nrow(object$data))
  }
  data <- object$data[rows, , drop = FALSE]
  group_vec <- object$group_vec[rows]

  if (type == "parametric") {
    # Parametric: keep the design fixed and regenerate the response by
    # drawing new noise (v) and inefficiency (u) terms from the fitted
    # group-specific error distributions
    boot_data <- .parametric_resample(object, data, group_vec)
  } else {
    # Nonparametric: case resampling, i.e. resample whole rows with
    # replacement within each group so group sizes are preserved
    boot_rows <- .nonparametric_resample(group_vec, groups)
    boot_data <- data[boot_rows, , drop = FALSE]
    rownames(boot_data) <- NULL
    group_vec <- group_vec[boot_rows]
  }

  # Re-fit the metafrontier with the fitted group labels rather than a
  # column name: the fit may have been given 'group' as a vector, and
  # the data may hold an unrelated column called "group"
  boot_fit <- metafrontier(
    formula = formula,
    data = boot_data,
    group = group_vec,
    method = object$method,
    meta_type = object$meta_type,
    dist = if (object$method == "sfa" &&
               !is.null(object$group_models[[1]]$dist)) {
      object$group_models[[1]]$dist
    } else "hnormal",
    orientation = if (object$method == "dea" &&
                      !is.null(object$orientation)) {
      object$orientation
    } else "output",
    rts = if (object$method == "dea" &&
              !is.null(object$rts)) {
      object$rts
    } else "crs",
    # Forward the estimation choices of the original fit so the
    # bootstrap distribution reflects the reported point estimates.
    type = if (!is.null(object$type)) object$type else "radial",
    direction = if (!is.null(object$direction)) {
      object$direction
    } else "proportional",
    estimator = if (!is.null(object$estimator)) {
      object$estimator
    } else "bc88",
    objective = if (!is.null(object$objective)) {
      object$objective
    } else "lp",
    engine = if (!is.null(object$engine)) object$engine else "internal",
    ...
  )

  boot_fit$tgr
}


# data and group_vec are restricted to the rows used in the fit, so the
# response, the group labels and each group's fitted frontier line up
.parametric_resample <- function(object, data, group_vec) {
  groups <- object$groups
  formula <- object$formula

  if (inherits(formula, "Formula")) {
    f_base <- formula(formula, rhs = 1)
  } else {
    f_base <- formula
  }

  mf <- model.frame(f_base, data = data, na.action = na.omit)
  y <- model.response(mf)
  boot_y <- y  # will be overwritten per group

  for (g in groups) {
    idx <- which(group_vec == g)
    gm <- object$group_models[[g]]
    sigma_v <- gm$sigma_v
    sigma_u <- gm$sigma_u
    n_g <- length(idx)

    # Resample from estimated distributions
    v_new <- rnorm(n_g, mean = 0, sd = sigma_v)

    # Draw from the fitted inefficiency distribution
    dist_g <- gm$dist
    if (is.null(dist_g)) dist_g <- "hnormal"

    u_new <- switch(dist_g,
      hnormal = abs(rnorm(n_g, mean = 0, sd = sigma_u)),
      exponential = rexp(n_g, rate = 1 / sigma_u),
      tnormal = {
        # Truncated normal: draw from N(mu, sigma_u^2) truncated at 0
        mu_val <- if (!is.null(gm$mu_vec)) mean(gm$mu_vec)
                  else if ("mu" %in% names(gm$all_params)) gm$all_params["mu"]
                  else 0
        # Simple rejection sampling for truncated normal
        raw <- rnorm(n_g * 3, mean = mu_val, sd = sigma_u)
        raw <- raw[raw >= 0]
        if (length(raw) < n_g) {
          # Fallback: use abs() if rejection sampling doesn't give enough
          raw <- c(raw, abs(rnorm(n_g, mean = mu_val, sd = sigma_u)))
        }
        raw[seq_len(n_g)]
      },
      abs(rnorm(n_g, mean = 0, sd = sigma_u))  # default fallback
    )

    # Fitted frontier value = X %*% beta
    fitted_g <- gm$fitted
    boot_y[idx] <- fitted_g + v_new - u_new
  }

  # Replace response in data
  resp_name <- all.vars(f_base)[1]
  boot_data <- data
  boot_data[[resp_name]] <- boot_y
  boot_data
}


# Row indices of a case resample, stacked group by group
.nonparametric_resample <- function(group_vec, groups) {
  boot_rows <- integer(0)
  for (g in groups) {
    idx <- which(group_vec == g)
    boot_idx <- sample(idx, length(idx), replace = TRUE)
    boot_rows <- c(boot_rows, boot_idx)
  }
  boot_rows
}


# ---------- Internal: CI computation ----------

.boot_ci <- function(tgr_boot, tgr_orig, alpha, ci_type) {
  n <- ncol(tgr_boot)
  ci <- matrix(NA_real_, nrow = n, ncol = 2)

  if (ci_type == "percentile") {
    for (i in seq_len(n)) {
      vals <- tgr_boot[, i]
      vals <- vals[is.finite(vals)]
      if (length(vals) >= 2) {
        ci[i, ] <- quantile(vals, probs = c(alpha, 1 - alpha))
      }
    }
  } else {
    # BCa
    for (i in seq_len(n)) {
      vals <- tgr_boot[, i]
      vals <- vals[is.finite(vals)]
      if (length(vals) < 2) next

      # Bias correction
      z0 <- qnorm(mean(vals < tgr_orig[i]))

      # Acceleration (jackknife)
      n_boot <- length(vals)
      theta_dot <- mean(vals)
      diffs <- theta_dot - vals
      a <- sum(diffs^3) / (6 * (sum(diffs^2))^1.5)

      # Adjusted quantiles
      z_alpha <- qnorm(alpha)
      z_1alpha <- qnorm(1 - alpha)

      a1 <- pnorm(z0 + (z0 + z_alpha) / (1 - a * (z0 + z_alpha)))
      a2 <- pnorm(z0 + (z0 + z_1alpha) / (1 - a * (z0 + z_1alpha)))

      ci[i, ] <- quantile(vals, probs = c(a1, a2))
    }
  }

  colnames(ci) <- paste0(format(100 * c(alpha, 1 - alpha),
                                trim = TRUE, digits = 3), "%")
  ci
}


# Percentile intervals for a within-group statistic of the TGR.
# tgr_group labels the original estimates; boot_group labels the
# columns of tgr_boot (these differ for the nonparametric bootstrap).
.boot_ci_group <- function(tgr_boot, tgr_orig, tgr_group, boot_group,
                           groups, alpha, stat = c("mean", "median")) {
  stat <- match.arg(stat)
  stat_fun <- match.fun(stat)

  result <- data.frame(
    Group = groups,
    Estimate = NA_real_,
    Lower = NA_real_,
    Upper = NA_real_,
    stringsAsFactors = FALSE
  )

  for (j in seq_along(groups)) {
    g <- groups[j]
    result$Estimate[j] <- stat_fun(tgr_orig[tgr_group == g])

    # Group statistic per bootstrap replicate
    group_stats <- apply(tgr_boot[, boot_group == g, drop = FALSE], 1,
                         stat_fun)
    group_stats <- group_stats[is.finite(group_stats)]

    if (length(group_stats) >= 2) {
      result$Lower[j] <- quantile(group_stats, probs = alpha)
      result$Upper[j] <- quantile(group_stats, probs = 1 - alpha)
    }
  }

  names(result)[2:4] <- c(
    if (stat == "mean") "Mean_TGR" else "Median_TGR",
    paste0(format(100 * c(alpha, 1 - alpha), trim = TRUE, digits = 3), "%")
  )
  result
}


# ---------- S3 methods ----------

#' @export
print.boot_tgr <- function(x, digits = 4, ...) {
  cat("\nBootstrap TGR Confidence Intervals\n")
  cat("----------------------------------\n")
  cat("Type:          ", x$type, "\n")
  cat("CI method:     ", x$ci_type, "\n")
  cat("Replications:  ", x$R_effective, "/", x$R, "\n")
  cat("Level:         ", x$level, "\n\n")

  titles <- c(ci_group = "Group-level mean TGR:",
              ci_group_median = "Group-level median TGR:")
  for (tab in names(titles)) {
    cat(titles[[tab]], "\n", sep = "")
    ci_print <- x[[tab]]
    num_cols <- sapply(ci_print, is.numeric)
    ci_print[num_cols] <- lapply(ci_print[num_cols], round, digits = digits)
    print(ci_print, row.names = FALSE)
    cat("\n")
  }
  invisible(x)
}


#' @export
confint.boot_tgr <- function(object, parm, level, ...) {
  if (is.null(object$ci)) {
    stop("Observation-level intervals are not available for the ",
         "nonparametric bootstrap: each column of 'tgr_boot' is a ",
         "resampled unit, not a fixed DMU. Use the group-level intervals ",
         "in '$ci_group' and '$ci_group_median'.", call. = FALSE)
  }
  if (!missing(level) && level != object$level) {
    warning("Recomputing CI at a different level requires re-running ",
            "boot_tgr(). Returning CI at the original level = ",
            object$level, ".", call. = FALSE)
  }
  ci <- object$ci
  if (!missing(parm)) {
    ci <- ci[parm, , drop = FALSE]
  }
  ci
}


#' @export
plot.boot_tgr <- function(x, which = c("distribution", "ci"),
                          group = NULL, ...) {
  which <- match.arg(which)

  if (which == "distribution") {
    # Histogram of mean TGR per group across bootstrap reps
    groups_plot <- if (!is.null(group)) group else x$groups
    n_groups <- length(groups_plot)
    old_par <- graphics::par(mfrow = c(1, n_groups))
    on.exit(graphics::par(old_par))

    for (g in groups_plot) {
      group_means <- apply(x$tgr_boot[, x$boot_group == g, drop = FALSE],
                           1, mean)
      graphics::hist(group_means, main = paste("TGR:", g),
                     xlab = "Mean TGR", col = "lightblue", border = "white")
      graphics::abline(v = mean(x$tgr_original[x$group_vec == g]),
                       col = "red", lwd = 2)
      # Dashed lines at the group-level CI bounds (columns 3:4 of
      # ci_group; names are level-dependent, e.g. "2.5%"/"97.5%")
      if (!is.null(x$ci_group)) {
        bounds <- unlist(x$ci_group[x$ci_group$Group == g, 3:4])
        graphics::abline(v = bounds, col = "red", lwd = 1, lty = 2)
      }
    }
  } else {
    # CI plot per group
    ci <- x$ci_group
    n_g <- nrow(ci)
    graphics::plot(seq_len(n_g), ci$Mean_TGR,
                   ylim = range(ci[, 3:4], na.rm = TRUE),
                   xaxt = "n", xlab = "Group", ylab = "Mean TGR",
                   pch = 19, main = "Bootstrap CI for Mean TGR")
    graphics::axis(1, at = seq_len(n_g), labels = ci$Group)
    graphics::segments(seq_len(n_g), ci[, 3], seq_len(n_g), ci[, 4],
                       lwd = 2)
  }
  invisible(x)
}
