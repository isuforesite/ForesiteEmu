#' Title
#'
#' @param x
#' @param nLHS
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
lhs_weather <- function(x, nLHS, seed) {
  xWide    <- reshape_curves_to_wide(x, "yday", "year")
  xLHSWide <- copula_lhs(x = xWide[, -1], n = nLHS, seed = seed)
  colnames(xLHSWide) <- colnames(xWide[, -1])

  xLHS     <- reshape_lhs_to_long(xLHSWide, xWide)
  colnames(xLHS)[9] <- "weather_sample_id"

  xLHS
}
