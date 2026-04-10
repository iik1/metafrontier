#' Simulate Metafrontier Data
#'
#' Generates synthetic data from a known metafrontier data-generating
#' process. Useful for Monte Carlo simulations, package testing, and
#' teaching.
#'
#' @param n_groups integer. Number of technology groups (default 2).
#' @param n_per_group integer or integer vector. Number of observations
#'   per group. If a single value, the same number is used for all
#'   groups. If a vector, must be of length \code{n_groups}.
#' @param n_inputs integer. Number of input variables (default 2).
#' @param beta_meta numeric vector. Metafrontier coefficients
#'   (including intercept). Length must be \code{n_inputs + 1}.
#'   Default: \code{c(1.0, 0.5, 0.3)}.
#' @param tech_gap numeric vector of length \code{n_groups}. The
#'   technology gap for each group, defined as the reduction in the
#'   intercept relative to the metafrontier. Default: evenly spaced
#'   from 0 to 0.5.
#' @param sigma_u numeric vector of length \code{n_groups}. Standard
#'   deviation of the half-normal inefficiency term for each group.
#'   Default: \code{rep(0.3, n_groups)}.
#' @param sigma_v numeric. Standard deviation of the symmetric noise
#'   term. Default: 0.2.
#' @param seed integer or \code{NULL}. Random seed for
#'   reproducibility.
#'
#' @return A list with components:
#'   \describe{
#'     \item{data}{a data frame with columns \code{log_y}, \code{log_x1},
#'       \code{log_x2}, ..., \code{group}, and the true underlying values}
#'     \item{params}{a list of the true parameters used for generation}
#'   }
#'
#' @examples
#' sim <- simulate_metafrontier(n_groups = 3, n_per_group = 200,
#'                              sigma_u = c(0.2, 0.4, 0.3))
#' str(sim$data)
#' table(sim$data$group)
#'
#' # The true metafrontier coefficients
#' sim$params$beta_meta
#'
#' @export
simulate_metafrontier <- function(n_groups = 2L,
                                  n_per_group = 100L,
                                  n_inputs = 2L,
                                  beta_meta = NULL,
                                  tech_gap = NULL,
                                  sigma_u = NULL,
                                  sigma_v = 0.2,
                                  seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  # Default metafrontier coefficients
  if (is.null(beta_meta)) {
    beta_meta <- c(1.0, seq(0.5, 0.2, length.out = n_inputs))
  }
  if (length(beta_meta) != n_inputs + 1L) {
    stop("'beta_meta' must have length n_inputs + 1.", call. = FALSE)
  }

  # Default technology gaps
  if (is.null(tech_gap)) {
    tech_gap <- seq(0, 0.5, length.out = n_groups)
  }
  if (length(tech_gap) != n_groups) {
    stop("'tech_gap' must have length n_groups.", call. = FALSE)
  }

  # Default inefficiency SDs
  if (is.null(sigma_u)) {
    sigma_u <- rep(0.3, n_groups)
  }
  if (length(sigma_u) != n_groups) {
    stop("'sigma_u' must have length n_groups.", call. = FALSE)
  }

  # Handle per-group sample sizes
  if (length(n_per_group) == 1L) {
    n_per_group <- rep(n_per_group, n_groups)
  }
  if (length(n_per_group) != n_groups) {
    stop("'n_per_group' must be scalar or length n_groups.", call. = FALSE)
  }

  n_total <- sum(n_per_group)

  # Generate data
  frames <- vector("list", n_groups)
  group_names <- paste0("G", seq_len(n_groups))

  # Group-specific betas: shift intercept by tech_gap
  beta_groups <- vector("list", n_groups)

  for (g in seq_len(n_groups)) {
    n_g <- n_per_group[g]

    # Generate inputs (log scale, from uniform)
    X_inputs <- matrix(stats::runif(n_g * n_inputs, 0, 5),
                       nrow = n_g, ncol = n_inputs)
    colnames(X_inputs) <- paste0("log_x", seq_len(n_inputs))

    # Design matrix with intercept
    X <- cbind(1, X_inputs)

    # Group-specific frontier = metafrontier - tech_gap (intercept only)
    beta_g <- beta_meta
    beta_g[1] <- beta_meta[1] - tech_gap[g]
    beta_groups[[g]] <- beta_g

    # Frontier output
    frontier_y <- X %*% beta_g

    # Noise and inefficiency
    v <- stats::rnorm(n_g, 0, sigma_v)
    u <- abs(stats::rnorm(n_g, 0, sigma_u[g]))

    # Observed output
    log_y <- as.numeric(frontier_y + v - u)

    # True values
    true_te <- exp(-u)
    true_tgr <- exp(-tech_gap[g])  # constant for intercept-shift DGP

    df_g <- data.frame(
      X_inputs,
      log_y = log_y,
      group = group_names[g],
      true_te = true_te,
      true_tgr = true_tgr,
      true_te_star = true_te * true_tgr,
      true_u = u,
      true_v = v,
      stringsAsFactors = FALSE
    )

    frames[[g]] <- df_g
  }

  data <- do.call(rbind, frames)
  rownames(data) <- NULL
  data$group <- factor(data$group, levels = group_names)

  params <- list(
    beta_meta = beta_meta,
    beta_groups = setNames(beta_groups, group_names),
    tech_gap = setNames(tech_gap, group_names),
    sigma_u = setNames(sigma_u, group_names),
    sigma_v = sigma_v,
    n_groups = n_groups,
    n_per_group = setNames(n_per_group, group_names)
  )

  list(data = data, params = params)
}
