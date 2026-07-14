#' Convert a Fitted Frontier Model to Metafrontier Format
#'
#' Generic function that extracts the components needed by
#' \code{\link{metafrontier}} from a pre-fitted frontier model.
#' Methods are provided for \pkg{sfaR} (\code{"sfacross"}),
#' \pkg{frontier} (\code{"frontier"}), and \pkg{Benchmarking}
#' (\code{"Farrell"}) objects, as well as plain lists with the
#' required fields.
#'
#' @param x a fitted frontier model object.
#' @param ... additional arguments passed to methods.
#'
#' @return A list of class \code{"metafrontier_model"} with components:
#'   \code{beta}, \code{te}, \code{X}, \code{y}, \code{sigma_v},
#'   \code{sigma_u}, \code{logLik}, \code{hessian}, \code{n},
#'   \code{dist}.
#'
#' @details
#' \code{metafrontier(models = ...)} calls this function internally on
#' each supplied model, so fitted \pkg{sfaR} or \pkg{frontier} objects
#' can be passed to \code{metafrontier()} directly. Manual conversion
#' is only needed for hand-built list models. Converting an object that
#' has already been converted is a no-op, so it is safe to pass
#' converted objects to \code{metafrontier()} as well.
#'
#' Note that \code{Benchmarking::dea()} (\code{"Farrell"}) objects do
#' not store the inputs, outputs, or frontier coefficients, so the
#' converted model carries only efficiency scores and cannot be used
#' with \code{metafrontier(models = ...)}; use the formula interface
#' with \code{method = "dea"} instead. Converting a Farrell object
#' therefore raises a warning.
#'
#' @examples
#' # Using a named list:
#' mod <- as_metafrontier_model(list(
#'   coefficients = c("(Intercept)" = 2, log_x1 = 0.5, log_x2 = 0.3),
#'   efficiency = runif(50, 0.7, 1),
#'   X = matrix(rnorm(150), 50, 3),
#'   y = rnorm(50, 5)
#' ))
#' str(mod)
#'
#' @export
as_metafrontier_model <- function(x, ...) {
  UseMethod("as_metafrontier_model")
}


#' @export
as_metafrontier_model.metafrontier_model <- function(x, ...) {
  x
}


#' @export
as_metafrontier_model.sfacross <- function(x, ...) {
  structure(.extract_sfacross(x), class = "metafrontier_model")
}


#' @export
as_metafrontier_model.frontier <- function(x, ...) {
  structure(.extract_frontier(x), class = "metafrontier_model")
}


#' @export
as_metafrontier_model.Farrell <- function(x, ...) {
  warning("Farrell (Benchmarking) objects do not store inputs, outputs, ",
          "or frontier coefficients; the converted model cannot be used ",
          "with metafrontier(models = ...). Use the formula interface ",
          "with method = \"dea\" instead.", call. = FALSE)
  structure(.extract_benchmarking(x), class = "metafrontier_model")
}


#' @export
as_metafrontier_model.list <- function(x, ...) {
  required <- c("coefficients", "efficiency", "X", "y")
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0) {
    stop("List is missing required fields: ",
         paste(missing_fields, collapse = ", "), ".", call. = FALSE)
  }
  structure(list(
    beta = x$coefficients,
    te = x$efficiency,
    X = x$X,
    y = as.numeric(x$y),
    sigma_v = if (is.null(x$sigma_v)) NA_real_ else x$sigma_v,
    sigma_u = if (is.null(x$sigma_u)) NA_real_ else x$sigma_u,
    logLik = if (is.null(x$logLik)) NA_real_ else x$logLik,
    hessian = x$hessian,
    n = length(x$y),
    dist = if (is.null(x$dist)) "hnormal" else x$dist
  ), class = "metafrontier_model")
}


#' @export
as_metafrontier_model.default <- function(x, ...) {
  stop("Unsupported model class: ", paste(class(x), collapse = ", "),
       ". Supported: 'sfacross' (sfaR::sfacross), 'frontier' ",
       "(frontier::sfa), 'Farrell' (Benchmarking::dea), or a named ",
       "list with 'coefficients', 'efficiency', 'X', and 'y'.",
       call. = FALSE)
}


#' Extract from frontier::sfa object
#'
#' frontier::sfa() returns an object of class "frontier" whose
#' dataTable matrix holds id, t, and y in the first three columns
#' followed by the model matrix, and whose mleParam vector holds the
#' nb frontier coefficients followed by sigmaSq and gamma.
#' @keywords internal
#' @noRd
.extract_frontier <- function(model) {
  if (!requireNamespace("frontier", quietly = TRUE)) {
    stop("Package 'frontier' is required to extract from frontier objects.",
         call. = FALSE)
  }

  all_coef <- model$mleParam
  n_beta <- model$nb
  beta <- all_coef[seq_len(n_beta)]

  te <- as.numeric(frontier::efficiencies(model))

  dat <- model$dataTable
  y <- as.numeric(dat[, 3])
  X <- dat[, 3 + seq_len(n_beta), drop = FALSE]

  sigma_sq <- all_coef["sigmaSq"]
  gamma <- all_coef["gamma"]
  sigma_v <- sqrt(sigma_sq * (1 - gamma))
  sigma_u <- sqrt(sigma_sq * gamma)

  list(
    beta = beta,
    te = te,
    X = X,
    y = y,
    sigma_v = as.numeric(sigma_v),
    sigma_u = as.numeric(sigma_u),
    logLik = model$mleLogl,
    hessian = NULL,
    n = length(y),
    dist = "hnormal"
  )
}


#' Extract from Benchmarking::dea object
#' @keywords internal
#' @noRd
.extract_benchmarking <- function(model) {
  if (!requireNamespace("Benchmarking", quietly = TRUE)) {
    stop("Package 'Benchmarking' is required to extract from Farrell objects.",
         call. = FALSE)
  }

  eff <- as.numeric(model$eff)

  # Benchmarking DEA objects store peers and lambdas
  peers <- if (!is.null(model$PEERS)) model$PEERS else NULL
  lambda <- if (!is.null(model$LAMBDA)) model$LAMBDA else NULL

  # Attempt to retrieve X and Y from the Farrell object.
  # Benchmarking::dea() does not store X/Y in its return value,

  # so these will generally be NULL.  The efficiency scores are
  # still usable for direct comparison, but the object cannot be
  # passed through metafrontier(models = ...) which requires X,
  # y, and beta for metafrontier envelope construction.
  X_ref <- model$XREF
  Y_ref <- model$YREF

  list(
    beta = NULL,  # DEA is nonparametric
    te = eff,
    X = X_ref,
    y = if (!is.null(Y_ref)) as.numeric(Y_ref[, 1]) else NULL,
    sigma_v = NA_real_,
    sigma_u = NA_real_,
    logLik = NA_real_,
    hessian = NULL,
    n = length(eff),
    dist = NA_character_,
    peers = peers,
    lambda = lambda
  )
}
