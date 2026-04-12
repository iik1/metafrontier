#' Test Poolability of Group Frontiers
#'
#' Tests the null hypothesis that all groups share a common frontier
#' (i.e., the metafrontier coincides with all group frontiers) against
#' the alternative that group-specific frontiers differ. Uses a
#' likelihood ratio test.
#'
#' @param object a fitted \code{"metafrontier"} object with
#'   \code{method = "sfa"}.
#' @param ... additional arguments (currently unused).
#'
#' @return A list of class \code{"htest"} with components:
#'   \describe{
#'     \item{statistic}{the LR test statistic}
#'     \item{parameter}{degrees of freedom}
#'     \item{p.value}{p-value of the test}
#'     \item{method}{description of the test}
#'   }
#'
#' @details
#' The LR statistic is:
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
#' (and \eqn{\mu} for truncated-normal).
#'
#' @examples
#' set.seed(42)
#' sim <- simulate_metafrontier(n_groups = 2, n_per_group = 200,
#'                              tech_gap = c(0, 0.5))
#' fit <- metafrontier(log_y ~ log_x1 + log_x2,
#'                     data = sim$data, group = "group")
#' poolability_test(fit)
#'
#' @export
poolability_test <- function(object, ...) {
  if (!inherits(object, "metafrontier_sfa")) {
    stop("Poolability test requires an SFA-based metafrontier.",
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
    data.name = deparse(object$call),
    ll_pooled = ll_pooled,
    ll_groups = ll_groups
  )
  class(result) <- "htest"
  result
}
