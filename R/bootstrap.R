#' Bootstrap Confidence Intervals for the Technology Gap Ratio
#'
#' Computes bootstrap confidence intervals for TGR estimates from
#' a fitted metafrontier model. Supports both parametric (residual
#' resampling) and nonparametric (case resampling) bootstraps.
#'
#' @param object a \code{"metafrontier"} object.
#' @param R integer. Number of bootstrap replications (default 999).
#' @param type character. \code{"parametric"} resamples from estimated
#'   error distributions; \code{"nonparametric"} resamples rows within
#'   groups with replacement.
#' @param level numeric. Confidence level (default 0.95).
#' @param ci_type character. \code{"percentile"} (default) or
#'   \code{"bca"} (bias-corrected and accelerated).
#' @param seed optional integer seed for reproducibility.
#' @param progress logical. Show progress bar (default \code{TRUE}).
#' @param ncores integer. Number of CPU cores for parallel bootstrap
#'   (default 1, sequential). Requires the \code{parallel} package.
#' @param ... additional arguments passed to \code{\link{metafrontier}}.
#'
#' @return An object of class \code{"boot_tgr"} containing:
#'   \describe{
#'     \item{tgr_boot}{R x n matrix of bootstrapped TGR values}
#'     \item{tgr_original}{original TGR estimates}
#'     \item{ci}{n x 2 matrix of observation-level confidence intervals}
#'     \item{ci_group}{data frame of group-level mean TGR intervals}
#'     \item{R_effective}{number of successful replications}
#'     \item{R}{requested number of replications}
#'     \item{type}{bootstrap type used}
#'     \item{ci_type}{CI type used}
#'     \item{level}{confidence level}
#'   }
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

  # Compute CIs
  alpha <- (1 - level) / 2
  ci <- .boot_ci(tgr_boot, object$tgr, alpha, ci_type)

  # Group-level mean TGR CIs
  ci_group <- .boot_ci_group(tgr_boot, object$tgr,
                             object$group_vec, object$groups,
                             alpha, ci_type)

  out <- list(
    tgr_boot = tgr_boot,
    tgr_original = object$tgr,
    ci = ci,
    ci_group = ci_group,
    R_effective = R_effective,
    R = R,
    type = type,
    ci_type = ci_type,
    level = level,
    group_vec = object$group_vec,
    groups = object$groups
  )
  class(out) <- "boot_tgr"
  out
}


# ---------- Internal: single bootstrap replicate ----------

.boot_one_replicate <- function(object, type, ...) {
  data <- object$data
  formula <- object$formula
  group_vec <- object$group_vec
  groups <- object$groups

  if (type == "parametric") {
    # Parametric: resample residuals from estimated distributions
    boot_data <- .parametric_resample(object)
  } else {
    # Nonparametric: resample rows within groups
    boot_data <- .nonparametric_resample(data, group_vec, groups)
  }

  # Re-fit the metafrontier using the original group column name
  group_col <- if (!is.null(object$group_col)) object$group_col else "group"

  boot_fit <- metafrontier(
    formula = formula,
    data = boot_data,
    group = group_col,
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
    ...
  )

  boot_fit$tgr
}


.parametric_resample <- function(object) {
  data <- object$data
  group_vec <- object$group_vec
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


.nonparametric_resample <- function(data, group_vec, groups) {
  boot_rows <- integer(0)
  for (g in groups) {
    idx <- which(group_vec == g)
    boot_idx <- sample(idx, length(idx), replace = TRUE)
    boot_rows <- c(boot_rows, boot_idx)
  }
  boot_data <- data[boot_rows, , drop = FALSE]
  rownames(boot_data) <- NULL
  boot_data
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


.boot_ci_group <- function(tgr_boot, tgr_orig, group_vec, groups,
                           alpha, ci_type) {
  result <- data.frame(
    Group = groups,
    Mean_TGR = NA_real_,
    Lower = NA_real_,
    Upper = NA_real_,
    stringsAsFactors = FALSE
  )

  for (j in seq_along(groups)) {
    g <- groups[j]
    idx <- which(group_vec == g)
    result$Mean_TGR[j] <- mean(tgr_orig[idx])

    # Mean TGR per bootstrap replicate
    group_means <- apply(tgr_boot[, idx, drop = FALSE], 1, mean)
    group_means <- group_means[is.finite(group_means)]

    if (length(group_means) >= 2) {
      result$Lower[j] <- quantile(group_means, probs = alpha)
      result$Upper[j] <- quantile(group_means, probs = 1 - alpha)
    }
  }

  names(result)[3:4] <- paste0(format(100 * c(alpha, 1 - alpha),
                                      trim = TRUE, digits = 3), "%")
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

  cat("Group-level mean TGR:\n")
  ci_print <- x$ci_group
  num_cols <- sapply(ci_print, is.numeric)
  ci_print[num_cols] <- lapply(ci_print[num_cols], round, digits = digits)
  print(ci_print, row.names = FALSE)
  cat("\n")
  invisible(x)
}


#' @export
confint.boot_tgr <- function(object, parm, level, ...) {
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
      idx <- which(x$group_vec == g)
      group_means <- apply(x$tgr_boot[, idx, drop = FALSE], 1, mean)
      graphics::hist(group_means, main = paste("TGR:", g),
                     xlab = "Mean TGR", col = "lightblue", border = "white")
      graphics::abline(v = mean(x$tgr_original[idx]), col = "red", lwd = 2)
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
