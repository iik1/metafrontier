#' Test Poolability of Group Frontiers
#'
#' Tests the null hypothesis that all groups share a common frontier
#' (i.e., the metafrontier coincides with all group frontiers) against
#' the alternative that group-specific frontiers differ. For SFA-based
#' metafrontiers a likelihood ratio test is used; for DEA-based
#' metafrontiers a permutation test is used.
#'
#' @param object a fitted \code{"metafrontier"} object with
#'   \code{method = "sfa"} or \code{method = "dea"}.
#' @param B integer. Number of permutation replicates for the DEA
#'   permutation test (default 199). Ignored for SFA objects.
#' @param seed integer or \code{NULL}. Random seed for the DEA
#'   permutation test, for reproducibility. Ignored for SFA objects.
#' @param ... additional arguments (currently unused).
#'
#' @return A list of class \code{"htest"} with components:
#'   \describe{
#'     \item{statistic}{the test statistic (LR statistic for SFA; the
#'       mean technology gap, \eqn{\bar{S} = \mathrm{mean}(1 - TGR)},
#'       for DEA)}
#'     \item{parameter}{degrees of freedom (SFA) or the effective
#'       number of permutation replicates (DEA)}
#'     \item{p.value}{p-value of the test}
#'     \item{method}{description of the test}
#'   }
#'
#' @details
#' \strong{Likelihood ratio test (SFA).} The LR statistic is:
#' \deqn{LR = -2 [LL_{pooled} - \sum_j LL_j]}
#' where \eqn{LL_{pooled}} is the log-likelihood of the pooled (single
#' frontier) model and \eqn{LL_j} are the group-specific
#' log-likelihoods. Under H0, the statistic follows a chi-squared
#' distribution with degrees of freedom equal to
#' \eqn{df = k_{groups} - k_{pooled}}, where \eqn{k_{groups}} is the
#' total number of parameters across all group-specific models and
#' \eqn{k_{pooled}} is the number of parameters in the pooled model.
#' For \eqn{J} groups each with \eqn{p} frontier parameters plus
#' distributional parameters, this equals
#' \eqn{(J - 1) \times p_{total}} where \eqn{p_{total}} includes
#' frontier coefficients, \eqn{\sigma_v}, and \eqn{\sigma_u}
#' (and \eqn{\mu} for truncated-normal). This test requires a
#' likelihood and is therefore only available for SFA-based
#' metafrontiers.
#'
#' \strong{Permutation test (DEA).} DEA has no likelihood, so the
#' poolability hypothesis is assessed by a permutation test. Under the
#' null of a single pooled technology, group labels are exchangeable:
#' reassigning observations to groups at random should not
#' systematically change the distance between the group frontiers and
#' the metafrontier. The observed statistic is the mean technology gap,
#' \eqn{S_{obs} = \mathrm{mean}(1 - TGR_i)}, and its null distribution
#' is approximated by refitting the metafrontier on \code{B} random
#' permutations of the group labels. The p-value is
#' \eqn{(1 + \#\{S_b \ge S_{obs}\}) / (B + 1)}, following the
#' aggregate-efficiency inference logic of Simar and Zelenyuk (2007).
#' The smoothed subsampling approach of Kneip, Simar, and Wilson (2016)
#' is the asymptotically rigorous alternative for testing hypotheses in
#' nonparametric production models; the permutation test offered here
#' is a computationally simple approximation. The default \code{B = 199}
#' is a pragmatic choice; p-values have resolution \eqn{1/(B + 1)}, so
#' increase \code{B} for finer resolution.
#'
#' @references
#' Simar, L. and Zelenyuk, V. (2007). Statistical inference for
#' aggregates of Farrell-type efficiencies. \emph{Journal of Applied
#' Econometrics}, 22(7), 1367--1394. \doi{10.1002/jae.991}
#'
#' Kneip, A., Simar, L. and Wilson, P.W. (2016). Testing hypotheses
#' in nonparametric models of production. \emph{Journal of Business &
#' Economic Statistics}, 34(3), 435--456.
#' \doi{10.1080/07350015.2015.1049747}
#'
#' @examples
#' set.seed(42)
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 200,
#'                              tech_gap = c(0, 0.5))
#' fit <- metafrontier(log_y ~ log_x1 + log_x2,
#'                     data = sim$data, group = "group")
#' poolability_test(fit)
#'
#' \donttest{
#' # DEA permutation test
#' fit_dea <- metafrontier(log_y ~ log_x1 + log_x2,
#'                         data = sim$data, group = "group",
#'                         method = "dea")
#' poolability_test(fit_dea, B = 99, seed = 1)
#' }
#'
#' @export
poolability_test <- function(object, B = 199, seed = NULL, ...) {
  data_name <- paste(deparse(substitute(object)), collapse = "")

  if (inherits(object, "metafrontier_dea")) {
    return(.poolability_permutation_dea(object, B = B, seed = seed,
                                        data_name = data_name))
  }

  if (!inherits(object, "metafrontier_sfa")) {
    stop("Poolability test requires an SFA- or DEA-based metafrontier.",
         call. = FALSE)
  }

  # Fit a pooled model (single frontier, no groups)
  pooled <- .fit_sfa_group(
    formula = object$formula,
    data = object$data,
    dist = object$group_models[[1]]$dist,
    control = list()
  )

  ll_pooled <- pooled$logLik
  ll_groups <- sum(object$logLik_groups)

  # LR statistic
  lr_stat <- -2 * (ll_pooled - ll_groups)

  # Degrees of freedom: difference in number of parameters
  k_pooled <- length(pooled$all_params)
  k_groups <- sum(sapply(object$group_models, function(m) {
    length(m$all_params)
  }))
  df <- k_groups - k_pooled

  p_value <- pchisq(lr_stat, df = df, lower.tail = FALSE)

  result <- list(
    statistic = c(`LR` = lr_stat),
    parameter = c(df = df),
    p.value = p_value,
    method = "Likelihood Ratio Test for Poolability of Group Frontiers",
    data.name = data_name,
    ll_pooled = ll_pooled,
    ll_groups = ll_groups
  )
  class(result) <- "htest"
  result
}


#' Permutation test for poolability of DEA-based metafrontiers
#' @noRd
.poolability_permutation_dea <- function(object, B, seed, data_name) {
  if (!is.null(seed)) set.seed(seed)

  s_obs <- mean(1 - object$tgr, na.rm = TRUE)

  # type/direction were not stored on objects fitted by older versions
  type <- if (!is.null(object$type)) object$type else "radial"
  direction <- if (!is.null(object$direction)) {
    object$direction
  } else {
    "proportional"
  }

  group_vec <- object$group_vec
  s_perm <- rep(NA_real_, B)

  for (b in seq_len(B)) {
    perm <- sample(group_vec)
    fit_b <- tryCatch(
      suppressWarnings(
        metafrontier(formula = object$formula,
                     data = object$data,
                     group = perm,
                     method = "dea",
                     orientation = object$orientation,
                     rts = object$rts,
                     type = type,
                     direction = direction)
      ),
      error = function(e) NULL
    )
    if (!is.null(fit_b)) {
      s_perm[b] <- mean(1 - fit_b$tgr, na.rm = TRUE)
    }
  }

  n_fail <- sum(is.na(s_perm))
  if (n_fail > 0.1 * B) {
    warning(n_fail, " of ", B, " permutation refits failed; the ",
            "p-value is based on the ", B - n_fail,
            " successful replicates only.", call. = FALSE)
  }

  s_perm <- s_perm[!is.na(s_perm)]
  b_eff <- length(s_perm)
  if (b_eff == 0L) {
    stop("All permutation refits failed; cannot compute a p-value.",
         call. = FALSE)
  }

  p_value <- (1 + sum(s_perm >= s_obs)) / (b_eff + 1)

  result <- list(
    statistic = c(`mean technology gap` = s_obs),
    parameter = c(B = b_eff),
    p.value = p_value,
    method = "Permutation test for poolability of group frontiers (DEA)",
    data.name = data_name
  )
  class(result) <- "htest"
  result
}
