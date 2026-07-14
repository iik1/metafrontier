#' Simulate Metafrontier Data
#'
#' Generates synthetic data from a known metafrontier data-generating
#' process. Each group frontier lies weakly below the metafrontier,
#' consistent with groups facing different restrictions of a common
#' metatechnology (Battese, Rao and O'Donnell, 2004). Useful for Monte
#' Carlo simulations, package testing, and teaching.
#'
#' @param n_groups integer. Number of technology groups (default 2).
#' @param n_per_group integer or integer vector. Number of observations
#'   per group. If a single value, the same number is used for all
#'   groups. If a vector, must be of length \code{n_groups}.
#' @param n_inputs integer. Number of input variables (default 2).
#' @param beta_meta numeric vector. Metafrontier coefficients
#'   (including intercept). Length must be \code{n_inputs + 1}.
#'   Default: \code{c(1.0, seq(0.5, 0.2, length.out = n_inputs))},
#'   i.e. \code{c(1.0, 0.5, 0.2)} for the default two inputs. Ignored
#'   when \code{beta_groups} is supplied.
#' @param tech_gap numeric vector of length \code{n_groups}. The
#'   technology gap for each group, defined as the reduction in the
#'   intercept relative to the metafrontier. Default: evenly spaced
#'   from 0 to 0.5. Ignored (with a warning) when \code{beta_groups}
#'   is supplied.
#' @param sigma_u numeric vector of length \code{n_groups}. Standard
#'   deviation of the half-normal inefficiency term for each group.
#'   Default: \code{rep(0.3, n_groups)}.
#' @param sigma_v numeric. Standard deviation of the symmetric noise
#'   term. Default: 0.2.
#' @param seed integer or \code{NULL}. Random seed for
#'   reproducibility.
#' @param beta_groups optional group-specific frontier coefficients,
#'   including slopes: either an \code{n_groups} x \code{(n_inputs + 1)}
#'   numeric matrix (one row per group) or a list of \code{n_groups}
#'   numeric vectors of length \code{n_inputs + 1}. When supplied, it
#'   replaces the intercept-shift construction based on
#'   \code{tech_gap}; see Details. Default \code{NULL} (intercept-shift
#'   design).
#' @param input_means optional \code{n_groups} x \code{n_inputs} numeric
#'   matrix of per-group mean log-input levels. When supplied, the
#'   log-inputs for group \code{g} are drawn from a normal distribution
#'   centred at \code{input_means[g, ]}; see Details. Default
#'   \code{NULL} (identical uniform inputs across groups).
#' @param input_corr optional \code{n_inputs} x \code{n_inputs}
#'   correlation matrix for the log-inputs. When supplied, the
#'   log-inputs are drawn from a multivariate normal distribution with
#'   this correlation structure; see Details. Default \code{NULL}
#'   (independent inputs).
#'
#' @details
#' By default the group frontiers share the metafrontier slopes and
#' differ only in their intercepts, so the true technology gap ratio
#' (TGR) is constant within each group and equals
#' \code{exp(-tech_gap[g])}. When \code{beta_groups} is supplied the
#' group frontiers may differ in their slopes, in which case no single
#' log-linear metafrontier envelops all groups: the tightest log-linear
#' envelope is then a pseudo-true quantity. The returned
#' \code{true_tgr} is instead computed observation by observation
#' against the pointwise maximum over the group frontiers,
#' \eqn{TGR_i = \exp(x_i^	op eta_g - \max_j x_i^	op eta_j)}, which is
#' guaranteed to lie in (0, 1]. The true group frontier for each firm
#' is \eqn{x_i^	op eta_g}, \code{true_te} is generated exactly as in the
#' default design, and \code{true_te_star = true_te * true_tgr}. In
#' this case \code{params$beta_meta} is \code{NULL} and
#' \code{params$beta_groups} holds the supplied coefficients.
#'
#' By default the log-inputs are drawn i.i.d. from a uniform
#' distribution on \code{[0, 5]}, identically across groups. Supplying
#' \code{input_means} and/or \code{input_corr} switches to normal
#' log-inputs with standard deviation \code{5 / sqrt(12)} (matching the
#' spread of the uniform draws), centred at \code{input_means[g, ]}
#' (2.5 for every group and input when \code{input_means} is
#' \code{NULL}). When \code{input_corr} is supplied the draws are
#' multivariate normal with that correlation matrix; when it is
#' \code{NULL} but \code{input_means} is given, the inputs are drawn
#' independently.
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
#' # Group-specific slopes: per-observation true TGR
#' sim2 <- simulate_metafrontier(
#'   beta_groups = rbind(c(1.0, 0.5, 0.2), c(0.9, 0.6, 0.1))
#' )
#' range(sim2$data$true_tgr)
#'
#' @export
simulate_metafrontier <- function(n_groups = 2L,
                                  n_per_group = 100L,
                                  n_inputs = 2L,
                                  beta_meta = NULL,
                                  tech_gap = NULL,
                                  sigma_u = NULL,
                                  sigma_v = 0.2,
                                  seed = NULL,
                                  beta_groups = NULL,
                                  input_means = NULL,
                                  input_corr = NULL) {

  if (!is.null(seed)) set.seed(seed)

  custom_betas <- !is.null(beta_groups)

  if (custom_betas) {
    if (!is.null(tech_gap)) {
      warning("'tech_gap' is ignored when 'beta_groups' is supplied.",
              call. = FALSE)
      tech_gap <- NULL
    }
    if (is.matrix(beta_groups)) {
      if (nrow(beta_groups) != n_groups ||
          ncol(beta_groups) != n_inputs + 1L) {
        stop("'beta_groups' must be an n_groups x (n_inputs + 1) matrix.",
             call. = FALSE)
      }
      beta_groups <- lapply(seq_len(n_groups),
                            function(g) as.numeric(beta_groups[g, ]))
    } else if (is.list(beta_groups)) {
      if (length(beta_groups) != n_groups ||
          !all(lengths(beta_groups) == n_inputs + 1L)) {
        stop("'beta_groups' must be a list of n_groups vectors of length n_inputs + 1.",
             call. = FALSE)
      }
      beta_groups <- lapply(beta_groups, as.numeric)
    } else {
      stop("'beta_groups' must be a matrix or a list of numeric vectors.",
           call. = FALSE)
    }
    beta_meta <- NULL
  } else {
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

  # Input distribution: legacy i.i.d. uniform, or (correlated) normal
  gaussian_inputs <- !is.null(input_means) || !is.null(input_corr)
  sd_x <- 5 / sqrt(12)  # matches the spread of runif(0, 5)
  chol_corr <- NULL
  if (gaussian_inputs) {
    if (is.null(input_means)) {
      input_means <- matrix(2.5, nrow = n_groups, ncol = n_inputs)
    }
    input_means <- as.matrix(input_means)
    if (nrow(input_means) != n_groups || ncol(input_means) != n_inputs ||
        !is.numeric(input_means)) {
      stop("'input_means' must be an n_groups x n_inputs numeric matrix.",
           call. = FALSE)
    }
    if (!is.null(input_corr)) {
      input_corr <- as.matrix(input_corr)
      if (nrow(input_corr) != n_inputs || ncol(input_corr) != n_inputs ||
          !isSymmetric(unname(input_corr)) ||
          any(abs(diag(input_corr) - 1) > 1e-8)) {
        stop("'input_corr' must be a symmetric n_inputs x n_inputs correlation matrix with unit diagonal.",
             call. = FALSE)
      }
      chol_corr <- tryCatch(
        chol(input_corr),
        error = function(e) {
          stop("'input_corr' must be positive definite: Cholesky factorisation failed.",
               call. = FALSE)
        }
      )
    }
  }

  # Generate data
  frames <- vector("list", n_groups)
  group_names <- paste0("G", seq_len(n_groups))

  # Group-specific betas: shift intercept by tech_gap (unless supplied)
  if (!custom_betas) {
    beta_groups <- vector("list", n_groups)
  }

  for (g in seq_len(n_groups)) {
    n_g <- n_per_group[g]

    # Generate inputs (log scale)
    if (gaussian_inputs) {
      Z <- matrix(stats::rnorm(n_g * n_inputs),
                  nrow = n_g, ncol = n_inputs)
      if (!is.null(chol_corr)) Z <- Z %*% chol_corr
      X_inputs <- sweep(Z * sd_x, 2, input_means[g, ], "+")
    } else {
      X_inputs <- matrix(stats::runif(n_g * n_inputs, 0, 5),
                         nrow = n_g, ncol = n_inputs)
    }
    colnames(X_inputs) <- paste0("log_x", seq_len(n_inputs))

    # Design matrix with intercept
    X <- cbind(1, X_inputs)

    if (custom_betas) {
      beta_g <- beta_groups[[g]]
      # Frontier of every group at these inputs, for the pointwise envelope
      XB <- matrix(vapply(beta_groups,
                          function(b) as.numeric(X %*% b),
                          numeric(n_g)),
                   nrow = n_g)
      frontier_y <- XB[, g]
    } else {
      # Group-specific frontier = metafrontier - tech_gap (intercept only)
      beta_g <- beta_meta
      beta_g[1] <- beta_meta[1] - tech_gap[g]
      beta_groups[[g]] <- beta_g

      # Frontier output
      frontier_y <- X %*% beta_g
    }

    # Noise and inefficiency
    v <- stats::rnorm(n_g, 0, sigma_v)
    u <- abs(stats::rnorm(n_g, 0, sigma_u[g]))

    # Observed output
    log_y <- as.numeric(frontier_y + v - u)

    # True values
    true_te <- exp(-u)
    if (custom_betas) {
      # TGR against the pointwise maximum over group frontiers, in (0, 1]
      true_tgr <- exp(frontier_y - apply(XB, 1, max))
    } else {
      true_tgr <- exp(-tech_gap[g])  # constant for intercept-shift DGP
    }

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
    tech_gap = if (custom_betas) NULL else setNames(tech_gap, group_names),
    sigma_u = setNames(sigma_u, group_names),
    sigma_v = sigma_v,
    n_groups = n_groups,
    n_per_group = setNames(n_per_group, group_names)
  )

  list(data = data, params = params)
}


#' Simulate Panel Metafrontier Data
#'
#' @description
#' Generates a simulated panel dataset with known metafrontier
#' parameters for Monte Carlo studies and testing. Implements
#' the Battese-Coelli 1992 DGP with time-varying inefficiency.
#'
#' @param n_groups integer. Number of technology groups.
#' @param n_firms_per_group integer. Number of firms per group.
#' @param n_periods integer. Number of time periods.
#' @param beta_meta numeric vector. Metafrontier coefficients.
#' @param tech_gap numeric vector of group technology gaps.
#' @param sigma_u numeric. Standard deviation of the firm effect.
#' @param sigma_v numeric. Standard deviation of noise.
#' @param eta numeric. Time-decay parameter for BC92.
#' @param seed integer or NULL. Random seed.
#' @param attrition numeric in [0, 0.5]. Probability that each
#'   firm-period observation after a firm's first period is dropped
#'   independently, producing an unbalanced panel. Every firm's first
#'   period is always kept, so all firms remain in the data. The
#'   attrition draws are made after all other random numbers, so
#'   \code{attrition = 0} (the default) reproduces legacy balanced
#'   datasets exactly for the same seed. The realised share of at-risk
#'   observations dropped is stored in \code{params$attrition_share}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{data}{data frame with columns: \code{firm}, \code{year},
#'       \code{group}, \code{log_y}, \code{log_x1}, \code{log_x2},
#'       \code{true_te}, \code{true_u}, \code{true_v}}
#'     \item{params}{list of true parameter values used in generation,
#'       including \code{attrition} and the realised
#'       \code{attrition_share}}
#'   }
#'
#' @examples
#' sim <- simulate_panel_metafrontier(
#'   n_groups = 2, n_firms_per_group = 20,
#'   n_periods = 5, seed = 42
#' )
#' head(sim$data)
#' str(sim$params)
#'
#' # An unbalanced panel with roughly 20% attrition
#' sim_unbal <- simulate_panel_metafrontier(seed = 42, attrition = 0.2)
#' table(table(sim_unbal$data$firm))
#'
#' @export
simulate_panel_metafrontier <- function(n_groups = 2,
                                        n_firms_per_group = 30,
                                        n_periods = 5,
                                        beta_meta = c(1.0, 0.5, 0.3),
                                        tech_gap = NULL,
                                        sigma_u = 0.3,
                                        sigma_v = 0.2,
                                        eta = 0.05,
                                        seed = NULL,
                                        attrition = 0) {

  if (!is.numeric(attrition) || length(attrition) != 1L ||
      is.na(attrition) || attrition < 0 || attrition > 0.5) {
    stop("'attrition' must be a single number in [0, 0.5].", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  n_inputs <- length(beta_meta) - 1
  if (is.null(tech_gap)) {
    tech_gap <- seq(0, 0.5, length.out = n_groups)
  }
  if (length(sigma_u) == 1) sigma_u <- rep(sigma_u, n_groups)

  group_names <- paste0("G", seq_len(n_groups))
  frames <- list()

  for (j in seq_len(n_groups)) {
    beta_j <- beta_meta
    beta_j[1] <- beta_j[1] - tech_gap[j]

    # Firm-level effects
    u_i <- abs(rnorm(n_firms_per_group, 0, sigma_u[j]))

    for (f in seq_len(n_firms_per_group)) {
      for (t in seq_len(n_periods)) {
        # Time-varying inefficiency: u_it = u_i * exp(-eta*(t-T))
        u_it <- u_i[f] * exp(-eta * (t - n_periods))

        # Input draws
        x_vals <- rnorm(n_inputs, mean = 1, sd = 0.5)
        x_row <- c(1, x_vals)

        # Frontier value
        frontier_val <- sum(beta_j * x_row)
        v_it <- rnorm(1, 0, sigma_v)
        y_it <- frontier_val + v_it - u_it

        row <- data.frame(
          log_y = y_it,
          firm = paste0(group_names[j], "_F", f),
          year = t,
          group = group_names[j],
          true_u = u_it,
          true_te = exp(-u_it),
          stringsAsFactors = FALSE
        )
        for (inp in seq_len(n_inputs)) {
          row[[paste0("log_x", inp)]] <- x_vals[inp]
        }
        frames <- c(frames, list(row))
      }
    }
  }

  data <- do.call(rbind, frames)
  rownames(data) <- NULL
  data$group <- factor(data$group, levels = group_names)

  # True TGR
  X_mat <- cbind(1, as.matrix(data[, paste0("log_x", seq_len(n_inputs))]))
  meta_frontier <- X_mat %*% beta_meta
  group_frontier <- numeric(nrow(data))
  for (j in seq_len(n_groups)) {
    idx <- data$group == group_names[j]
    beta_j <- beta_meta
    beta_j[1] <- beta_j[1] - tech_gap[j]
    group_frontier[idx] <- X_mat[idx, ] %*% beta_j
  }
  data$true_tgr <- exp(group_frontier - meta_frontier)

  # Attrition: drop post-first-period observations independently. Drawn
  # after all other random numbers so attrition = 0 reproduces legacy
  # datasets exactly for the same seed.
  attrition_share <- 0
  if (attrition > 0 && n_periods > 1) {
    at_risk <- which(data$year > 1)
    dropped <- at_risk[stats::runif(length(at_risk)) < attrition]
    attrition_share <- length(dropped) / length(at_risk)
    if (length(dropped) > 0) {
      data <- data[-dropped, , drop = FALSE]
      rownames(data) <- NULL
    }
  }

  params <- list(
    beta_meta = beta_meta,
    tech_gap = setNames(tech_gap, group_names),
    sigma_u = setNames(sigma_u, group_names),
    sigma_v = sigma_v,
    eta = eta,
    n_groups = n_groups,
    n_firms_per_group = n_firms_per_group,
    n_periods = n_periods,
    attrition = attrition,
    attrition_share = attrition_share
  )

  list(data = data, params = params)
}
