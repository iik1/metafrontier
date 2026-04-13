#' Convert a Fitted Frontier Model to Metafrontier Format
#'
#' Generic function that extracts the components needed by
#' \code{\link{metafrontier}} from a pre-fitted frontier model.
#' Methods are provided for \pkg{sfaR}, \pkg{frontier}, and
#' \pkg{Benchmarking} objects, as well as plain lists with the
#' required fields.
#'
#' @param x a fitted frontier model object.
#' @param ... additional arguments passed to methods.
#'
#' @return A list with components: \code{coefficients}, \code{efficiency},
#'   \code{X}, \code{y}, \code{sigma_v}, \code{sigma_u}, \code{logLik},
#'   \code{hessian}, \code{n}, \code{dist}.
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
as_metafrontier_model.sfacross <- function(x, ...) {
  .extract_sfacross(x)
}


#' @export
as_metafrontier_model.sfa <- function(x, ...) {
  .extract_frontier_sfa(x)
}


#' @export
as_metafrontier_model.Farrell <- function(x, ...) {
  .extract_benchmarking(x)
}


#' @export
as_metafrontier_model.list <- function(x, ...) {
  required <- c("coefficients", "efficiency", "X", "y")
  missing_fields <- setdiff(required, names(x))
  if (length(missing_fields) > 0) {
    stop("List is missing required fields: ",
         paste(missing_fields, collapse = ", "), ".", call. = FALSE)
  }
  list(
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
  )
}


#' @export
as_metafrontier_model.default <- function(x, ...) {
  stop("Unsupported model class: ", paste(class(x), collapse = ", "),
       ". Supported: sfaR::sfacross, frontier::sfa, Benchmarking::dea, ",
       "or a named list with 'coefficients', 'efficiency', 'X', and 'y'.",
       call. = FALSE)
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
