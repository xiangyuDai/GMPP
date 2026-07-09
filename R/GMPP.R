#' @title Genomic mating to predict expected progeny genetic values
#'
#' @description
#' Predicts expected progeny genetic values for all possible crosses
#' between sires and dams under purebred, crossbred, or multi-population
#' mating designs.
#'
#' @param mating_sire    Numeric matrix. Row names = sire IDs
#'   (if \code{NULL}, \code{"1", "2", ...} are assigned).
#'   Column 1 = phenotype, remaining columns = SNP genotypes coded
#'   as \code{{0, 1, 2}}.
#' @param mating_dam     Numeric matrix in the same format as
#'   \code{mating_sire}, for dams.
#' @param add_sire       Optional. Additional sire population, same
#'   format as \code{mating_sire}. Default \code{NULL}. When provided
#'   without \code{add_dam}, produces (\code{mating_sire} \eqn{\times} \code{mating_dam})
#'   \eqn{\times} \code{add_sire} EPVs; when both are provided, produces
#'   (\code{mating_sire} \eqn{\times} \code{mating_dam}) \eqn{\times} (\code{add_sire} \eqn{\times} \code{add_dam})
#'   EPVs.
#' @param add_dam        Optional. Additional dam population, same
#'   format as \code{mating_sire}. Default \code{NULL}. When provided
#'   without \code{add_sire}, produces (\code{mating_sire} \eqn{\times} \code{mating_dam})
#'   \eqn{\times} \code{add_dam} EPVs; when both are provided, produces
#'   (\code{mating_sire} \eqn{\times} \code{mating_dam}) \eqn{\times} (\code{add_sire} \eqn{\times} \code{add_dam})
#'   EPVs.
#' @param ref_population Optional reference population. Numeric
#'   matrix. Column 1 = phenotype, remaining columns = SNP genotypes
#'   \code{{0, 1, 2}}. When provided,
#'   marker effects are estimated from this population rather than
#'   from the mating parents. Default \code{NULL}.
#' @param method         Genetic model. Must include \code{"additive"}.
#'   \code{c("additive", "dominant")} adds dominant effects. Default
#'   \code{"additive"}.
#' @param lmm.sire       Named list of arguments passed to
#'   \code{\link[gaston]{lmm.aireml}} for \code{mating_sire}.
#'   Default \code{list()}.
#' @param lmm.dam        Named list for \code{mating_dam}.
#' @param lmm.add_sire   Named list for \code{add_sire}.
#' @param lmm.add_dam    Named list for \code{add_dam}.
#' @param lmm.ref        Named list for \code{ref_population}.
#' @param lmm.purebred   Named list for the combined purebred
#'   population (only used when \code{crossbred = FALSE}).
#' @param crossbred      Logical. If \code{TRUE} (default), sires and
#'   dams are treated as distinct populations (crossbred). If
#'   \code{FALSE}, they are pooled as a single purebred population.
#' @param verbose        Logical. If \code{TRUE} (default), progress
#'   messages are printed to the console.
#' @param block_size     Integer. Number of SNPs processed per
#'   iteration.  Larger values are faster but use more memory.
#'   Default \code{100}.
#'
#' @return Expected progeny values. Format depends on the mating design:
#'   \describe{
#'     \item{Case 1 (no add_*):}{Numeric matrix, rows = sires,
#'       columns = dams.}
#'     \item{Case 2 (one add_*):}{Numeric matrix, rows = sire\eqn{\times}dam
#'       (\code{"sire|dam"}), columns = additional individuals.}
#'     \item{Case 3 (both add_*):}{Numeric matrix, rows = sire\eqn{\times}dam,
#'       columns = add_sire\eqn{\times}add_dam
#'       (\code{"add_sire|add_dam"}).}
#'   }
#'
#' @references
#' Aliloo, H., Pryce, J.E., Gonzalez-Recio, O., Cocks, B.G.,
#' Goddard, M.E. and Hayes, B.J. (2017). Including nonadditive
#' genetic effects in mating programs to maximize dairy farm
#' profitability. \emph{Journal of Dairy Science}, \bold{100}(2),
#' 1203-1222.
#'
#' @importFrom gaston lmm.aireml
#' @export
#'
#' @seealso \code{\link{GRM}} for genomic relationship matrix
#'   computation, \code{\link{Marker_effect}} for marker effect
#'   estimation.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' n_sire <- 10; n_dam <- 20; n_snp <- 100
#'
#' # Simulate sire population
#' sire_geno <- matrix(sample(0:2, n_sire * n_snp, replace = TRUE),
#'                     n_sire, n_snp)
#' rownames(sire_geno) <- paste0("S", 1:n_sire)
#' beta_sire <- rnorm(n_snp, 0, 1)
#' g_sire <- as.vector(sire_geno %*% beta_sire)
#' sire_data <- cbind(pheno = g_sire + rnorm(n_sire, 0, 1), sire_geno)
#'
#' # Simulate dam population
#' dam_geno <- matrix(sample(0:2, n_dam * n_snp, replace = TRUE),
#'                    n_dam, n_snp)
#' rownames(dam_geno) <- paste0("D", 1:n_dam)
#' beta_dam <- rnorm(n_snp, 0, 1)
#' g_dam <- as.vector(dam_geno %*% beta_dam)
#' dam_data <- cbind(pheno = g_dam + rnorm(n_dam, 0, 1), dam_geno)
#'
#' result <- GMPP(mating_sire = sire_data, mating_dam = dam_data,
#'                    method = "additive", crossbred = TRUE)
#' head(result)
#' }
GMPP <- function(mating_sire, mating_dam, add_sire = NULL, add_dam = NULL,
               ref_population = NULL, method = "additive",
               lmm.sire = list(), lmm.dam = list(),
               lmm.add_sire = list(), lmm.add_dam = list(),
               lmm.ref = list(), lmm.purebred = list(),
               crossbred = TRUE, verbose = TRUE, block_size = 100L) {

  # ---- welcome ----
  if (verbose) {
    line <- strrep("=", 64L)
    art <- c(
      "  ____    __  __    ____     ____  ",
      " / ___|  |  \\/  |  |  _ \\   |  _ \\ ",
      "| |  _   | |\\/| |  | |_) |  | |_) |",
      "| |_| |  | |  | |  |  __/   |  __/ ",
      " \\____|  |_|  |_|  |_|      |_|    "
    )
    max_w <- max(nchar(art))
    art <- vapply(art, function(l) paste0(l, strrep(" ", max_w - nchar(l))), "",
                  USE.NAMES = FALSE)
    art <- paste0(strrep(" ", 14L), art)
    art[5] <- paste0(art[5], strrep(" ", 64L - nchar(art[5]) - 6L), "v1.0.0")
    message(line, "\n",
            paste0(art, collapse = "\n"), "\n",
            line, "\n")
  }

  # ---- input validation ----
  method <- match.arg(method, several.ok = TRUE)
  if (!"additive" %in% method)
    stop("'method' must include \"additive\"")

  if (!is.logical(crossbred) || length(crossbred) != 1L || is.na(crossbred))
    stop("'crossbred' must be a single logical value (TRUE or FALSE)")
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("'verbose' must be a single logical value (TRUE or FALSE)")
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      is.na(block_size) || block_size < 1L)
    stop("'block_size' must be a single positive integer")
  if (block_size != as.integer(block_size))
    stop("'block_size' must be a whole number")

  validate_lmm <- function(x, name) {
    if (!is.list(x))
      stop("'", name, "' must be a list")
    if (length(x) > 0 && is.null(names(x)))
      stop("'", name, "' must be a named list")
    x$verbose <- verbose
    return(x)
  }
  lmm.sire     <- validate_lmm(lmm.sire,     "lmm.sire")
  lmm.dam      <- validate_lmm(lmm.dam,      "lmm.dam")
  lmm.add_sire <- validate_lmm(lmm.add_sire, "lmm.add_sire")
  lmm.add_dam  <- validate_lmm(lmm.add_dam,  "lmm.add_dam")
  lmm.ref      <- validate_lmm(lmm.ref,      "lmm.ref")
  lmm.purebred <- validate_lmm(lmm.purebred, "lmm.purebred")

  validate_pop <- function(x, name) {
    if (!is.matrix(x))
      stop("'", name, "' must be a numeric matrix")
    if (ncol(x) < 2L)
      stop("'", name, "' must have at least 2 columns (phenotype + genotypes), ",
           ncol(x), " found")
    if (!is.numeric(x[, 1]))
      stop("'", name, "' column 1 (phenotype) must be numeric")
    if (any(is.na(x[, 1])))
      stop("'", name, "' column 1 (phenotype) contains NA values")
    geno <- x[, -1, drop = FALSE]
    if (!is.numeric(geno))
      stop("'", name, "' genotype columns must be numeric")
    if (any(is.na(geno)))
      stop("'", name, "' genotype columns contain NA values")
    if (!all(geno %in% c(0, 1, 2)))
      stop("'", name, "' genotype columns contain values outside {0, 1, 2}")
    rn <- rownames(x)
    if (!is.null(rn) && anyDuplicated(rn))
      stop("'", name, "' has duplicate rownames")
  }

  validate_pop(mating_sire, "mating_sire")
  validate_pop(mating_dam,  "mating_dam")

  nsnp <- ncol(mating_sire) - 1L
  if (ncol(mating_dam) - 1L != nsnp)
    stop("'mating_sire' and 'mating_dam' must have the same number of SNPs (",
         nsnp, " vs ", ncol(mating_dam) - 1L, ")")

  if (!is.null(add_sire)) {
    validate_pop(add_sire, "add_sire")
    if (ncol(add_sire) - 1L != nsnp)
      stop("'add_sire' must have the same number of SNPs as 'mating_sire' (",
           nsnp, " expected, ", ncol(add_sire) - 1L, " found)")
  }
  if (!is.null(add_dam)) {
    validate_pop(add_dam, "add_dam")
    if (ncol(add_dam) - 1L != nsnp)
      stop("'add_dam' must have the same number of SNPs as 'mating_sire' (",
           nsnp, " expected, ", ncol(add_dam) - 1L, " found)")
  }
  if (!is.null(ref_population)) {
    validate_pop(ref_population, "ref_population")
    if (ncol(ref_population) - 1L != nsnp)
      stop("'ref_population' must have the same number of SNPs as 'mating_sire' (",
           nsnp, " expected, ", ncol(ref_population) - 1L, " found)")
  }

  ns <- nrow(mating_sire); nd <- nrow(mating_dam)

  # ---- header ----
  if (verbose) {
    design <- if (crossbred) "crossbred" else "purebred"
    mod <- paste(method, collapse = " + ")
    n_cross  <- ns * nd
    n_total  <- if (!is.null(add_sire)) n_cross * nrow(add_sire) else n_cross
    n_total  <- if (!is.null(add_dam))  n_total * nrow(add_dam) else n_total

    message(paste0(strrep("-", 24), " Genomic Mating ", strrep("-", 24)))
    message(sprintf(
      "  Design: %s  |  Model: %s  |  Markers: %d",
      design, mod, nsnp))
    if (!is.null(ref_population)) message("  reference population: ", nrow(ref_population))
    message("  mating_sire: ", ns)
    message("  mating_dam: ", nd)
    if (!is.null(add_sire)) message("  add_sire: ", nrow(add_sire))
    if (!is.null(add_dam))  message("  add_dam: ", nrow(add_dam))
    message("  EPVs: ", n_total)
  }

  # ---- helper: fit marker effects with clean progress ----
  fit_marker <- function(Y, Ka, Kd = NULL, geno, lmm, label) {
    n <- length(Y)
    if (verbose) message(sprintf("  - %-20s (n=%d) ", label, n), appendLF = FALSE)
    gaston_out <- utils::capture.output({
      eff <- Marker_effect(Y = Y, Ka = Ka, Kd = Kd, geno = geno,
                           LMM.control = lmm)
    })
    if (verbose) {
      iter_lines <- grep("^\\[Iteration \\d+\\]", gaston_out, value = TRUE)
      if (length(iter_lines) == 0L) {
        message("done (0 iter)")
        warning("'", label, "' model fitting produced 0 iterations; ",
                "marker effect estimates may be unreliable")
        return(eff)
      }
      iter_nums  <- as.integer(gsub("^\\[Iteration (\\d+)\\].*", "\\1", iter_lines))
      last_iter  <- max(iter_nums)
      n_iter     <- length(unique(iter_nums))
      message(sprintf("done (%d iter)", n_iter))
      last_lines <- gaston_out[grep(sprintf("^\\[Iteration %d\\]", last_iter), gaston_out)]
      key_lines  <- grep("theta =|log L =|\\|\\|gradient", last_lines, value = TRUE)
      for (l in key_lines) message("    ", l)
    }
    return(eff)
  }

  # ---- step 1: estimation of genotypic values of progeny ----
  if (verbose) message("[1/2] Genotypic value estimation")

  geno_value <- list()

  if (!is.null(ref_population)) {
    p         <- colMeans(ref_population[, -1]) / 2
    G_Matrix  <- GRM(ref_population[, -1], method = method)
    snp_effect <- fit_marker(
      Y = ref_population[, 1], Ka = G_Matrix$Ga, Kd = G_Matrix$Gd,
      geno = ref_population[, -1], lmm = lmm.ref, label = "ref_population"
    )

    if (!is.null(snp_effect$d)) {
      geno_value$AA <- (2 - 2 * p) * snp_effect$a + (-2 * (1 - p)^2) * snp_effect$d
      geno_value$Aa <- (1 - 2 * p) * snp_effect$a + (2 * p * (1 - p)) * snp_effect$d
      geno_value$aa <- (0 - 2 * p) * snp_effect$a + (-2 * p^2) * snp_effect$d
    } else {
      geno_value$AA <- (2 - 2 * p) * snp_effect$a
      geno_value$Aa <- (1 - 2 * p) * snp_effect$a
      geno_value$aa <- (0 - 2 * p) * snp_effect$a
    }

  } else if (crossbred) {
    if ("dominant" %in% method)
      warning("Reference population = NULL; Under crossbred mating design, ",
              "genotypic values do not include dominance deviation effects.")

    p_sire <- colMeans(mating_sire[, -1]) / 2
    p_dam  <- colMeans(mating_dam[, -1]) / 2

    G_sire <- GRM(mating_sire[, -1], method = "additive")
    G_dam  <- GRM(mating_dam[, -1], method = "additive")

    sire_effect <- fit_marker(
      Y = mating_sire[, 1], Ka = G_sire$Ga,
      geno = mating_sire[, -1], lmm = lmm.sire, label = "mating_sire"
    )
    dam_effect <- fit_marker(
      Y = mating_dam[, 1], Ka = G_dam$Ga,
      geno = mating_dam[, -1], lmm = lmm.dam, label = "mating_dam"
    )

    mating_AA <- ((2 - 2 * p_sire) * sire_effect$a + (2 - 2 * p_dam) * dam_effect$a) / 2
    mating_Aa <- ((1 - 2 * p_sire) * sire_effect$a + (1 - 2 * p_dam) * dam_effect$a) / 2
    mating_aa <- ((0 - 2 * p_sire) * sire_effect$a + (0 - 2 * p_dam) * dam_effect$a) / 2

    if (is.null(add_sire) && is.null(add_dam)) {
      geno_value$AA <- mating_AA
      geno_value$Aa <- mating_Aa
      geno_value$aa <- mating_aa
    }

    if (xor(is.null(add_sire), is.null(add_dam))) {
      add_data  <- if (!is.null(add_sire)) add_sire else add_dam
      lmm_add   <- if (!is.null(add_sire)) lmm.add_sire else lmm.add_dam
      add_label  <- if (!is.null(add_sire)) "add_sire" else "add_dam"
      p_add     <- colMeans(add_data[, -1]) / 2
      G_add     <- GRM(add_data[, -1], method = "additive")
      add_effect <- fit_marker(
        Y = add_data[, 1], Ka = G_add$Ga,
        geno = add_data[, -1], lmm = lmm_add, label = add_label
      )
      add_AA <- (2 - 2 * p_add) * add_effect$a
      add_Aa <- (1 - 2 * p_add) * add_effect$a
      add_aa <- (0 - 2 * p_add) * add_effect$a
      geno_value$AA <- (mating_AA + add_AA) / 2
      geno_value$Aa <- (mating_Aa + add_Aa) / 2
      geno_value$aa <- (mating_aa + add_aa) / 2
    }

    if (!is.null(add_sire) && !is.null(add_dam)) {
      p_add_sire <- colMeans(add_sire[, -1]) / 2
      p_add_dam  <- colMeans(add_dam[, -1]) / 2
      G_add_sire <- GRM(add_sire[, -1], method = "additive")
      G_add_dam  <- GRM(add_dam[, -1], method = "additive")
      add_sire_effect <- fit_marker(
        Y = add_sire[, 1], Ka = G_add_sire$Ga,
        geno = add_sire[, -1], lmm = lmm.add_sire, label = "add_sire"
      )
      add_dam_effect <- fit_marker(
        Y = add_dam[, 1], Ka = G_add_dam$Ga,
        geno = add_dam[, -1], lmm = lmm.add_dam, label = "add_dam"
      )
      add_AA <- ((2 - 2 * p_add_sire) * add_sire_effect$a +
                 (2 - 2 * p_add_dam) * add_dam_effect$a) / 2
      add_Aa <- ((1 - 2 * p_add_sire) * add_sire_effect$a +
                 (1 - 2 * p_add_dam) * add_dam_effect$a) / 2
      add_aa <- ((0 - 2 * p_add_sire) * add_sire_effect$a +
                 (0 - 2 * p_add_dam) * add_dam_effect$a) / 2
      geno_value$AA <- (mating_AA + add_AA) / 2
      geno_value$Aa <- (mating_Aa + add_Aa) / 2
      geno_value$aa <- (mating_aa + add_aa) / 2
    }

  } else {
    # Purebred (crossbred = FALSE)
    if (!is.null(add_sire) || !is.null(add_dam))
      warning("'add_sire' and 'add_dam' are ignored when 'crossbred = FALSE'; ",
              "all individuals should be included in 'mating_sire' and 'mating_dam'")
    purebred <- rbind(mating_sire, mating_dam)
    p <- colMeans(purebred[, -1]) / 2
    G_purebred <- GRM(purebred[, -1], method = method)
    purebred_effect <- fit_marker(
      Y = purebred[, 1], Ka = G_purebred$Ga, Kd = G_purebred$Gd,
      geno = purebred[, -1], lmm = lmm.purebred, label = "mating_sire + mating_dam"
    )

    if (!is.null(purebred_effect$d)) {
      geno_value$AA <- (2 - 2 * p) * purebred_effect$a + (-2 * (1 - p)^2) * purebred_effect$d
      geno_value$Aa <- (1 - 2 * p) * purebred_effect$a + (2 * p * (1 - p)) * purebred_effect$d
      geno_value$aa <- (0 - 2 * p) * purebred_effect$a + (-2 * p^2) * purebred_effect$d
    } else {
      geno_value$AA <- (2 - 2 * p) * purebred_effect$a
      geno_value$Aa <- (1 - 2 * p) * purebred_effect$a
      geno_value$aa <- (0 - 2 * p) * purebred_effect$a
    }
  }

  # ---- step 2: calculation of expected progeny genetic values ----
  # Strategy: gamete A-allele probability = a = geno/2  (0, 0.5, 1)
  #   EPV_k = a1*a2*w2 + a1*w1 + a2*w1 + gaa_k
  #   where w1 = gAa - gaa,  w2 = gAA - 2*gAa + gaa
  # Accumulate per-individual statistics across SNP blocks, then assemble.
  if (verbose) message("[2/2] Progeny value computation")

  get_id <- function(x) {
    id <- rownames(x)
    if (is.null(id)) id <- as.character(seq_len(nrow(x)))
    return(id)
  }

  sire_id <- get_id(mating_sire)
  dam_id  <- get_id(mating_dam)
  ns <- nrow(mating_sire); nd <- nrow(mating_dam); nc <- ns * nd
  si <- rep(1:ns, times = nd)
  di <- rep(1:nd, each  = ns)
  sire_id_extract <- sire_id[si]
  dam_id_extract  <- dam_id[di]

  gAA <- geno_value$AA; gAa <- geno_value$Aa; gaa <- geno_value$aa
  if (block_size > nsnp) {
    warning("'block_size' (", block_size, ") exceeds number of SNPs (", nsnp,
            "), setting block_size = ", nsnp)
  }
  block_size <- min(block_size, nsnp)

  # ---- helper: extract gamete probabilities for a block of loci ----
  # gamete A-allele prob = geno/2
  gam_blk <- function(x, idx) x[, idx + 1L, drop = FALSE] / 2

  # ============================================================
  # Case 1: no additional populations
  # ============================================================
  if (is.null(add_sire) && is.null(add_dam)) {

    B     <- matrix(0, ns, nd)
    b_s   <- numeric(ns)
    b_d   <- numeric(nd)
    const <- 0

    for (start in seq(1L, nsnp, by = block_size)) {
      end <- min(start + block_size - 1L, nsnp)
      idx <- start:end
      w1 <- gAa[idx] - gaa[idx]
      w2 <- gAA[idx] - 2 * gAa[idx] + gaa[idx]

      a_s <- gam_blk(mating_sire, idx)
      a_d <- gam_blk(mating_dam,  idx)

      b_s   <- b_s   + as.vector(a_s %*% w1)
      b_d   <- b_d   + as.vector(a_d %*% w1)
      const <- const + sum(gaa[idx])
      B     <- B     + sweep(a_s, 2, w2, "*") %*% t(a_d)
    }

    EPV <- B + b_s + rep(b_d, each = ns) + const
    rownames(EPV) <- sire_id
    colnames(EPV) <- dam_id
    if (verbose) message(paste0(strrep("-", 29), " Done ", strrep("-", 29)))
    return(EPV)
  }

  # ============================================================
  # Case 2: exactly one additional population
  # ============================================================
  if (xor(is.null(add_sire), is.null(add_dam))) {

    add_data <- if (!is.null(add_sire)) add_sire else add_dam
    add_id   <- get_id(add_data)
    na <- nrow(add_data)

    B_sa  <- matrix(0, ns, na)
    B_da  <- matrix(0, nd, na)
    b_s   <- numeric(ns)
    b_d   <- numeric(nd)
    b_a   <- numeric(na)
    const <- 0

    for (start in seq(1L, nsnp, by = block_size)) {
      end <- min(start + block_size - 1L, nsnp)
      idx <- start:end
      w1 <- gAa[idx] - gaa[idx]
      w2 <- gAA[idx] - 2 * gAa[idx] + gaa[idx]

      a_s <- gam_blk(mating_sire, idx)
      a_d <- gam_blk(mating_dam,  idx)
      a_a <- gam_blk(add_data,    idx)

      b_s   <- b_s   + as.vector(a_s %*% w1)
      b_d   <- b_d   + as.vector(a_d %*% w1)
      b_a   <- b_a   + as.vector(a_a %*% w1)
      const <- const + sum(gaa[idx])
      B_sa  <- B_sa  + sweep(a_s, 2, w2, "*") %*% t(a_a)
      B_da  <- B_da  + sweep(a_d, 2, w2, "*") %*% t(a_a)
    }

    # Assemble from per-individual statistics
    #   B_sa[si,] expands ns*na -> nc*na by repeating rows (same size as output)
    EPV <- 0.5 * (B_sa[si, , drop = FALSE] +
                  B_da[di, , drop = FALSE] +
                  b_s[si] + b_d[di]) +
           rep(b_a, each = nc) + const
    rownames(EPV) <- paste(sire_id_extract, dam_id_extract, sep = "|")
    colnames(EPV) <- add_id
    if (verbose) message(paste0(strrep("-", 29), " Done ", strrep("-", 29)))
    return(EPV)
  }

  # ============================================================
  # Case 3: both additional populations
  # ============================================================
  if (!is.null(add_sire) && !is.null(add_dam)) {

    add_sire_id <- get_id(add_sire)
    add_dam_id  <- get_id(add_dam)
    ns_add <- nrow(add_sire); nd_add <- nrow(add_dam)
    nc_add <- ns_add * nd_add

    si_add <- rep(1:ns_add, times = nd_add)
    di_add <- rep(1:nd_add, each  = ns_add)
    add_sire_id_extract <- add_sire_id[si_add]
    add_dam_id_extract  <- add_dam_id[di_add]

    B_ss <- matrix(0, ns, ns_add)
    B_sd <- matrix(0, ns, nd_add)
    B_ds <- matrix(0, nd, ns_add)
    B_dd <- matrix(0, nd, nd_add)
    b_s   <- numeric(ns)
    b_d   <- numeric(nd)
    b_as  <- numeric(ns_add)
    b_ad  <- numeric(nd_add)
    const <- 0

    for (start in seq(1L, nsnp, by = block_size)) {
      end <- min(start + block_size - 1L, nsnp)
      idx <- start:end
      w1 <- gAa[idx] - gaa[idx]
      w2 <- gAA[idx] - 2 * gAa[idx] + gaa[idx]

      a_s  <- gam_blk(mating_sire, idx)
      a_d  <- gam_blk(mating_dam,  idx)
      a_as <- gam_blk(add_sire,    idx)
      a_ad <- gam_blk(add_dam,     idx)

      b_s   <- b_s   + as.vector(a_s  %*% w1)
      b_d   <- b_d   + as.vector(a_d  %*% w1)
      b_as  <- b_as  + as.vector(a_as %*% w1)
      b_ad  <- b_ad  + as.vector(a_ad %*% w1)
      const <- const + sum(gaa[idx])
      B_ss  <- B_ss  + sweep(a_s,  2, w2, "*") %*% t(a_as)
      B_sd  <- B_sd  + sweep(a_s,  2, w2, "*") %*% t(a_ad)
      B_ds  <- B_ds  + sweep(a_d,  2, w2, "*") %*% t(a_as)
      B_dd  <- B_dd  + sweep(a_d,  2, w2, "*") %*% t(a_ad)
    }

    # Precompute terms constant across add_sire x add_dam pairs
    sd_linear <- rep(b_s, nd) + rep(b_d, each = ns)

    # Assembly: column by column to avoid large temporaries
    EPV <- matrix(0, nc, nc_add)
    for (pq in seq_len(nc_add)) {
      p <- si_add[pq]
      q <- di_add[pq]

      col_pq <- 0.25 * (rep(B_ss[, p], nd) +
                        rep(B_sd[, q], nd) +
                        rep(B_ds[, p], each = ns) +
                        rep(B_dd[, q], each = ns)) +
                0.5 * (sd_linear + b_as[p] + b_ad[q]) + const
      EPV[, pq] <- col_pq
    }

    rownames(EPV) <- paste(sire_id_extract, dam_id_extract, sep = "|")
    colnames(EPV) <- paste(add_sire_id_extract, add_dam_id_extract, sep = "|")
    if (verbose) message(paste0(strrep("-", 29), " Done ", strrep("-", 29)))
    return(EPV)
  }
}
