#!/usr/bin/env Rscript

# Generate pls::plsr(method = "simpls") reference fixtures for multivar tests.
#
# Uses the public plsr API (not only the internal simpls.fit helper) so
# coefficients, intercepts, scores, loadings, training predictions, and
# newdata predictions stay aligned with R package pls.
#
# Prints a Scala object to stdout. Regenerate with:
#   Rscript tools/r-parity/generate_multivar_pls_r_parity_fixtures.R \
#     > modules/core/shared/src/test/scala/multivar/family/paired/PlsRegressionRReferenceFixtures.scala

suppressPackageStartupMessages(library(pls))

fmt_num <- function(x) {
  if (is.nan(x)) return("Double.NaN")
  if (is.infinite(x)) {
    return(if (x > 0) "Double.PositiveInfinity" else "Double.NegativeInfinity")
  }
  if (identical(x, 0) || identical(x, 0L)) return("0.0")
  formatC(as.numeric(x), digits = 17, format = "fg", flag = "#")
}

fmt_vec <- function(x) {
  paste0("Vector(", paste(vapply(as.numeric(x), fmt_num, character(1)), collapse = ", "), ")")
}

fmt_matrix <- function(x, indent = 4) {
  x <- as.matrix(x)
  if (ncol(x) == 0L || nrow(x) == 0L) {
    stop("cannot format empty matrix")
  }
  pad <- paste(rep(" ", indent), collapse = "")
  row_pad <- paste(rep(" ", indent + 4), collapse = "")
  rows <- apply(x, 1, function(row) paste0(row_pad, fmt_vec(row)))
  paste0(
    "GaleNumerics.matrixFromRows(\n",
    pad, "Vector(\n",
    paste(rows, collapse = ",\n"),
    "\n",
    pad, ")\n",
    paste(rep(" ", indent - 2), collapse = ""), ")"
  )
}

fmt_dvec <- function(x, indent = 4) {
  paste0("GaleNumerics.vectorFromArray(Array(", paste(vapply(as.numeric(x), fmt_num, character(1)), collapse = ", "), "))")
}

canonicalize_latent <- function(projection, scores, xloadings, yloadings) {
  projection <- as.matrix(projection)
  scores <- as.matrix(scores)
  xloadings <- as.matrix(xloadings)
  yloadings <- as.matrix(yloadings)
  for (col in seq_len(ncol(projection))) {
    best <- which.max(abs(projection[, col]))
    if (projection[best, col] < 0) {
      projection[, col] <- -projection[, col]
      scores[, col] <- -scores[, col]
      xloadings[, col] <- -xloadings[, col]
      yloadings[, col] <- -yloadings[, col]
    }
  }
  list(
    projection = projection,
    scores = scores,
    xloadings = xloadings,
    yloadings = yloadings
  )
}

plsr_case <- function(name, x, y, ncomp, new_x, seed_note) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  new_x <- as.matrix(new_x)
  storage.mode(x) <- "double"
  storage.mode(y) <- "double"
  storage.mode(new_x) <- "double"

  df <- data.frame(Y = I(y), X = I(x))
  fit <- plsr(Y ~ X, ncomp = ncomp, data = df, method = "simpls", validation = "none")

  coef_mat <- drop(coef(fit, ncomp = ncomp, intercept = FALSE))
  if (is.null(dim(coef_mat))) {
    coef_mat <- matrix(coef_mat, ncol = 1L)
  }
  intercept <- as.numeric(fit$Ymeans - as.numeric(fit$Xmeans %*% coef_mat))
  train_pred <- drop(predict(fit, ncomp = ncomp))
  if (is.null(dim(train_pred))) {
    train_pred <- matrix(train_pred, ncol = 1L)
  }
  new_pred <- drop(predict(fit, newdata = data.frame(X = I(new_x)), ncomp = ncomp))
  if (is.null(dim(new_pred))) {
    new_pred <- matrix(new_pred, ncol = 1L)
  }

  canon <- canonicalize_latent(
    fit$projection,
    scores(fit),
    loadings(fit),
    Yloadings(fit)
  )

  list(
    name = name,
    seed_note = seed_note,
    x = x,
    y = y,
    new_x = new_x,
    ncomp = ncomp,
    coefficients = coef_mat,
    intercept = intercept,
    x_rotations = canon$projection,
    x_scores = canon$scores,
    x_loadings = canon$xloadings,
    y_loadings = canon$yloadings,
    train_predictions = train_pred,
    new_predictions = new_pred,
    x_means = as.numeric(fit$Xmeans),
    y_means = as.numeric(fit$Ymeans)
  )
}

emit_case <- function(case) {
  c(
    paste0("  /** ", case$seed_note, " */"),
    paste0("  val ", case$name, ": PlsPlsrReference ="),
    "    PlsPlsrReference(",
    paste0("      x = ", fmt_matrix(case$x, 6), ","),
    paste0("      y = ", fmt_matrix(case$y, 6), ","),
    paste0("      newX = ", fmt_matrix(case$new_x, 6), ","),
    paste0("      components = ", case$ncomp, ","),
    paste0("      coefficients = ", fmt_matrix(case$coefficients, 6), ","),
    paste0("      intercept = ", fmt_dvec(case$intercept), ","),
    paste0("      xRotations = ", fmt_matrix(case$x_rotations, 6), ","),
    paste0("      xScores = ", fmt_matrix(case$x_scores, 6), ","),
    paste0("      xLoadings = ", fmt_matrix(case$x_loadings, 6), ","),
    paste0("      yLoadings = ", fmt_matrix(case$y_loadings, 6), ","),
    paste0("      trainingPredictions = ", fmt_matrix(case$train_predictions, 6), ","),
    paste0("      newPredictions = ", fmt_matrix(case$new_predictions, 6)),
    "    )"
  )
}

set.seed(42)
n <- 12
p <- 4
m <- 2
ncomp <- 2
X <- matrix(rnorm(n * p), n, p)
beta <- cbind(c(1.0, -0.5, 0.25, 0.1), c(0.2, 0.8, -0.3, 0.4))
Y <- X %*% beta + matrix(rnorm(n * m, sd = 0.1), n, m)
newX <- matrix(rnorm(5 * p), 5, p)

multivariate <- plsr_case(
  "multivariate",
  X,
  Y,
  ncomp,
  newX,
  "Synthetic multivariate response with m < p (eigen of S'S branch). set.seed(42)."
)

univariate <- plsr_case(
  "univariate",
  X,
  Y[, 1, drop = FALSE],
  ncomp,
  newX,
  "Synthetic univariate response with two components. set.seed(42), same X as multivariate."
)

set.seed(7)
nw <- 15
pw <- 3
mw <- 6
Xw <- matrix(rnorm(nw * pw), nw, pw)
Yw <- matrix(rnorm(nw * mw), nw, mw)
newXw <- matrix(rnorm(4 * pw), 4, pw)
wide <- plsr_case(
  "wideResponse",
  Xw,
  Yw,
  2,
  newXw,
  "Synthetic wide response with m > p (eigen of SS' branch). set.seed(7)."
)

# Compact real-data check: yarn density on a NIR column subset.
data(yarn)
yarn_x <- yarn$NIR[, c(10L, 40L, 80L, 120L, 160L, 200L), drop = FALSE]
yarn_y <- matrix(as.numeric(yarn$density), ncol = 1L)
yarn_new <- yarn_x[c(1L, 8L, 15L, 22L), , drop = FALSE]
yarn_case <- plsr_case(
  "yarnSubset",
  yarn_x,
  yarn_y,
  2,
  yarn_new,
  "pls::yarn density ~ six NIR wavelengths, method='simpls', ncomp=2."
)

cat("package multivar\n")
cat("package family.paired\n\n")
cat("import multivar.core.*\n\n")
cat("import gale.linalg.DMat\n")
cat("import gale.linalg.DVec\n\n")
cat("/** Golden PLS fixtures from tools/r-parity/generate_multivar_pls_r_parity_fixtures.R.\n")
cat("  *\n")
cat("  * Provenance: R package pls ", as.character(packageVersion("pls")), ";\n", sep = "")
cat("  * public API `plsr(..., method = \"simpls\", validation = \"none\")`.\n")
cat("  * Latent signs are canonicalized so the largest-magnitude x-rotation entry\n")
cat("  * is positive; coefficients and predictions are sign-invariant.\n")
cat("  */\n")
cat("object PlsRegressionRReferenceFixtures:\n")
cat("  final case class PlsPlsrReference(\n")
cat("      x: DMat,\n")
cat("      y: DMat,\n")
cat("      newX: DMat,\n")
cat("      components: Int,\n")
cat("      coefficients: DMat,\n")
cat("      intercept: DVec,\n")
cat("      xRotations: DMat,\n")
cat("      xScores: DMat,\n")
cat("      xLoadings: DMat,\n")
cat("      yLoadings: DMat,\n")
cat("      trainingPredictions: DMat,\n")
cat("      newPredictions: DMat\n")
cat("  )\n\n")
writeLines(emit_case(multivariate))
cat("\n")
writeLines(emit_case(univariate))
cat("\n")
writeLines(emit_case(wide))
cat("\n")
writeLines(emit_case(yarn_case))
cat("\n")
