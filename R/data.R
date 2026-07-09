#' Simulated sire data
#'
#' Simulated genotypes and phenotypes for 10 sires.
#' Loaded via \code{data(simulated_mate)}.
#'
#' @format A numeric matrix with 10 rows and 101 columns.
#'   Column 1 is phenotype, columns 2--101 are SNP genotypes coded as \code{{0, 1, 2}}.
#'   Row names: \code{"S1"} to \code{"S10"}.
#' @source Simulated data.
#' @examples
#' data(simulated_mate)
#' dim(sire_data)
#' @name sire_data
NULL

#' Simulated dam data
#'
#' Simulated genotypes and phenotypes for 20 dams.
#' Loaded via \code{data(simulated_mate)}.
#'
#' @format A numeric matrix with 20 rows and 101 columns.
#'   Column 1 is phenotype, columns 2--101 are SNP genotypes coded as \code{{0, 1, 2}}.
#'   Row names: \code{"D1"} to \code{"D20"}.
#' @source Simulated data.
#' @examples
#' data(simulated_mate)
#' dim(dam_data)
#' @name dam_data
NULL
