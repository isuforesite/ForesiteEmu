#' Title
#'
#' @param x
#' @param unitKey
#' @param soilLayers
#' @param nLHS
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
lhs_soil <- function(x, unitKey, soilLayers, nLHS, seed) {
  xFilled  <- fill_curves(x, unitKey, soilLayers)

  xWide    <- reshape_curves_to_wide(xFilled, "bottom", unitKey)
  xLHSWide <- copula_lhs(x = xWide[, -1], n = nLHS, seed = seed)
  colnames(xLHSWide) <- colnames(xWide[, -1])

  xLHSLong <- reshape_lhs_to_long(xLHSWide, xWide)
  colnames(xLHSLong)[1]  <- "layer"
  colnames(xLHSLong)[16] <- "soil_sample_id"

  xLHSLong
}

# Soil profiles -----------------------------------------------------------
fill_curves <- function(x, key, thresholds) {
  mTh     <- max(thresholds)
  kSplits <- split(x, x[, key])
  l       <- lapply(kSplits, function(df) {
    do.call(
      rbind,
      lapply(thresholds, function(th) {
        if (min(df$hzdept_r) != 0) {
          dfNew          <- df[which.min(df$hzdept_r), ]
          dfNew$hzdept_r <- 0
          dfNew$hzdepb_r <- min(df$hzdept_r)
          df             <- rbind(dfNew, df)
        }

        if (max(df$hzdepb_r) < mTh) {
          dfNew          <- df[which.max(df$hzdepb_r), ]
          dfNew$hzdept_r <- max(df$hzdept_r)
          dfNew$hzdepb_r <- mTh
          df             <- rbind(dfNew, df)
        }

        ind <- df$hzdept_r <= th
        w   <- which.max(df[ind, "hzdepb_r"])
        cbind(df[ind, ][w, ], bottom = th)
      })
    )
  })

  out <- do.call(rbind, l)
  rownames(out) <- NULL

  out
}

# Data handling -----------------------------------------------------------
reshape_curves_to_wide <- function(x, time, key) {
  reshape(x, timevar = time, idvar = key, direction = "wide")
}

reshape_lhs_to_long <- function(x, wideDF) {
  l     <- attr(wideDF, "reshapeWide")
  dfOut <- reshape(
    as.data.frame(x),
    timevar   = l$timevar,
    idvar     = l$idvar,
    times     = l$times,
    varying   = l$varying,
    direction = "long"
  )

  for (t in l$times)
    colnames(dfOut) <- gsub(paste(".", t, sep = ""), "", colnames(dfOut))

  rownames(dfOut) <- NULL

  dfOut
}