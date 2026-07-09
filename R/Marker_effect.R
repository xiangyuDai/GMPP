#' @title Estimate marker effects from a linear mixed model
#'
#' @description
#' Estimates additive (and optionally dominant) SNP marker effects
#' using GBLUP with \pkg{gaston}.
#'
#' @param Y           Phenotype vector.
#' @param Ka          Additive genomic relationship matrix.
#' @param Kd          Optional dominant genomic relationship matrix.
#'                    If \code{NULL}, only additive effects are estimated.
#' @param geno        Numeric matrix of genotypes coded as \code{{0, 1, 2}}
#'                    (individuals in rows, SNPs in columns).
#' @param LMM.control Named list of additional arguments passed to
#'                    \code{\link[gaston]{lmm.aireml}}, such as
#'                    \code{X} (design matrix, defaults to intercept),
#'                    \code{theta}, \code{EMsteps}, \code{max_iter},
#'                    \code{eps}, \code{verbose}, etc.
#'                    Parameters not set here use \pkg{gaston}'s defaults.
#'
#' @return A named list with components:
#'   \describe{
#'     \item{a}{Additive marker effects (vector).}
#'     \item{d}{Dominant marker effects (vector, if Kd is provided).}
#'   }
#'
#' @importFrom gaston lmm.aireml
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n_y <- 10; n_snp <- 100
#' geno <- matrix(sample(0:2, n_y * n_snp, replace = TRUE), n_y, n_snp)
#' beta <- rnorm(n_snp, 0, 1)
#' g <- as.vector(geno %*% beta)
#' y <- g + rnorm(n_y, 0, 1)
#' G <- GRM(geno, c("additive", "dominant"))
#' eff <- Marker_effect(Y = y, Ka = G$Ga, Kd = G$Gd, geno = geno)
#' }
Marker_effect <- function(Y, Ka, Kd = NULL, geno,
                          LMM.control = list()) {

  # ---- input validation ----

  # Y: phenotype
  if (!is.vector(Y) && !is.matrix(Y))
    stop("'Y' must be a matrix or numeric vector")
  if (any(is.na(Y)))
    stop("'Y' contains NA values")
  if (is.matrix(Y)) {
    if (ncol(Y) != 1)
      stop("'Y' should have only one column (", ncol(Y), " found) or a vector")
  }
  n <- length(Y)

  # Ka: additive GRM
  if (!is.matrix(Ka))
    stop("'Ka' must be a matrix")
  if (!is.numeric(Ka))
    stop("'Ka' must be numeric")
  if (nrow(Ka) != ncol(Ka))
    stop("'Ka' must be a square matrix")
  if (nrow(Ka) != n)
    stop("'Ka' dimensions (", nrow(Ka), " x ", ncol(Ka),
         ") do not match Y length (", n, ")")

  # Kd: dominant GRM (optional)
  if (!is.null(Kd)) {
    if (!is.matrix(Kd))
      stop("'Kd' must be a matrix")
    if (!is.numeric(Kd))
      stop("'Kd' must be numeric")
    if (nrow(Kd) != ncol(Kd))
      stop("'Kd' must be a square matrix")
    if (nrow(Kd) != n)
      stop("'Kd' dimensions (", nrow(Kd), " x ", ncol(Kd),
           ") do not match Y length (", n, ")")
  }

  # geno: genotype matrix
  if (!is.matrix(geno))
    stop("'geno' must be a matrix")
  if (!is.numeric(geno))
    stop("'geno' must be numeric")
  if (ncol(geno) < 1L)
    stop("'geno' must have at least 1 column")
  if (nrow(geno) != n)
    stop("'geno' rows (", nrow(geno),
         ") do not match Y length (", n, ")")
  if (any(is.na(geno)))
    stop("'geno' contains NA values")
  if (!all(geno %in% c(0, 1, 2)))
    stop("'geno' contains values outside {0, 1, 2}")

  # LMM.control
  if (!is.list(LMM.control))
    stop("'LMM.control' must be a list")
  if (length(LMM.control) > 0 && is.null(names(LMM.control)))
    stop("'LMM.control' must be a named list")

  # ---- fit linear mixed model ----
  K_list <- if (is.null(Kd)) list(Ka) else list(Ka, Kd)

  fit <- do.call(gaston::lmm.aireml, c(
    list(Y = Y, K = K_list),
    LMM.control
  ))

  # ---- helper ----
  Py     <- fit$Py
  Adjust <- function(x) x * (1 - x)
  effect <- list()

  # ---- allele frequency ----
  P <- colMeans(geno) / 2

  # ---- additive marker effects ----
  Za       <- sweep(geno, 2, 2 * P)
  Adj      <- 2 * sum(Adjust(P))
  sigma_a  <- fit$tau[1]
  effect$a <- (crossprod(Za, Py) * sigma_a) / Adj

  # ---- dominant marker effects ----
  if (!is.null(Kd)) {
    Zd        <- sweep(geno == 2, 2, -2 * (1 - P)^2, '*') +
                 sweep(geno == 1, 2,  2 * P * (1 - P), '*') +
                 sweep(geno == 0, 2, -2 * P^2, '*')
    Adj_d     <- sum((2 * Adjust(P))^2)
    sigma_d   <- fit$tau[2]
    effect$d  <- (crossprod(Zd, Py) * sigma_d) / Adj_d
  }

  return(effect)
}
