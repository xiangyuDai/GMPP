#' @title Calculate additive and dominant genomic relationship matrix
#'
#' @param geno   Numeric matrix of genotypes coded as \code{{0, 1, 2}}
#'               (individuals in rows, SNPs in columns).
#' @param method Character vector specifying which GRM(s) to compute.
#'               One or both of \code{"additive"} and \code{"dominant"}.
#' @return A named list with components:
#'   \describe{
#'     \item{Ga}{Additive GRM (if \code{"additive"} \%in\% method).}
#'     \item{Gd}{Dominant GRM (if \code{"dominant"} \%in\% method).}
#'   }
#'
#' @export
#'
#' @examples
#' geno <- matrix(sample(0:2, 200 * 50, replace = TRUE), 200, 50)
#' G <- GRM(geno)
#' G <- GRM(geno, "additive")
#' G <- GRM(geno, c("additive", "dominant"))
#' 

GRM <- function(geno, method = c("additive", "dominant")) {

  # ---- input validation ----
  if (!is.matrix(geno))
    stop("'geno' must be a matrix")
  if (ncol(geno) < 1L)
    stop("'geno' must have at least 1 column")
  if (!is.numeric(geno))
    stop("'geno' must be numeric")
  if (any(is.na(geno)))
    stop("'geno' contains NA values")
  if (!all(geno %in% c(0, 1, 2)))
    stop("'geno' contains values outside {0, 1, 2}")

  method <- match.arg(method, several.ok = TRUE)

  # ---- helper ----
  G      <- list()
  p      <- colMeans(geno) / 2
  Adjust <- function(x) x * (1 - x)

  # ---- additive GRM ----
  if ("additive" %in% method) {
    A    <- sweep(geno, 2, 2 * p)
    Adj  <- 2 * sum(Adjust(p))
    G$Ga <- tcrossprod(A) / Adj
  }

  # ---- dominant GRM ----
  if ("dominant" %in% method) {
    D     <- sweep(geno == 2, 2, -2 * (1 - p)^2, '*') +
             sweep(geno == 1, 2,  2 * p * (1 - p), '*') +
             sweep(geno == 0, 2, -2 * p^2, '*')
    Adj_d <- sum((2 * Adjust(p))^2)
    G$Gd  <- tcrossprod(D) / Adj_d
  }

  return(G)
}
