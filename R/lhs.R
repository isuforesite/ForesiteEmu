#' Generate a LHS from a given matrix.
#'
#' @param X
#' @param n
#' @param seed
#'
#' @return
#' @export
#'
#' @import lhs
lhs <- function(X, n, seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)

  lh <- lhs::randomLHS(n, ncol(X))
  lx <- do.call(
    cbind,
    lapply(1:ncol(X), function(i) {
      quantile(X[, i], probs = lh[, i])
    })
  )

  rownames(lx) <- NULL
  colnames(lx) <- colnames(X)

  lx
}

# lhs ---------------------------------------------------------------------
# If family = NULL, use empirical copula
copula_lhs <- function(x, n, family = NULL, seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)

  # Fit copula
  if (is.null(family)) {
    cop  <- copula::empCopula(copula::pobs(x))
  } else {
    cop  <- copula::archmCopula(family, dim = ncol(x))
    fit  <- copula::fitCopula(cop, copula::pobs(x), start = 1)
    cop  <- copula::archmCopula(family, dim = ncol(x), param = fit@estimate)
  }

  # Sample from fitted copula
  u    <- copula::rCopula(n, cop)

  # Transform iid uniform to LH uniform
  uHL  <- copula::rLatinHypercube(u)

  # Transform LH uniform to original scale
  xUL  <- do.call(
    cbind,
    lapply(1:ncol(x), function(i) { quantile(x[, i], probs = uHL[, i]) })
  )

  colnames(xUL) <- NULL
  rownames(xUL) <- NULL

  xUL
}

# LHS sampling paired with SVD --------------------------------------------
fsvd_lhs <- function(x, n, seed = NULL) {
  dec <- svd(x) # u contains the values in the uncorrelated space
  lx  <- lhs(dec$u, n, seed)
  lX  <- as.matrix(lx) %*% diag(dec$d) %*% t(dec$v)

  rownames(lX) <- NULL
  colnames(lX) <- colnames(x)

  lX
}

fdiff_lhs <- function(x, n, seed = NULL) {
  di <- x[, 2:ncol(x)] - x[, 1:(ncol(x) - 1)]

  lx  <- lhs(cbind(x[, 1, drop = FALSE], di), n, seed)
  lX  <- t(apply(lx, 1, cumsum))

  rownames(lX) <- NULL
  colnames(lX) <- colnames(x)

  lX
}
