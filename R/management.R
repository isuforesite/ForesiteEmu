#' Title
#'
#' @param x
#' @param nlhs
#' @param seed
#'
#' @return
#' @export
#'
#' @import data.table
lhs_management <- function(x, nlhs, seed) {
  # Set up data -------------------------------------------------------------
  if (!("data.table" %in% class(x)))
    x <- data.table::data.table(x)

  # Constants:
  #   tillage.implement, tillage.date, fertilizer.formulation
  # Discrete:
  #   tillage.depth, fertilizer.crop, fertilizer.rotation,
  #   fertilizer.application_method, fertilizer.app1_date, fertilizer.app2_date,
  #   crops.cultivar, crops.sow_crop, crops.sowing_density, crops.sowing_depth,
  #   crops.row_spacing, crops.harvest_crop.
  # Continuous:
  #   tillage.residue, fertilizer.app1_kg_n_ha, fertilizer.app2_kg_n_ha,
  #   crops.planting_date

  discreteCols <- c(
    "tillage.implement", "tillage.date", "fertilizer.formulation",
    "tillage.depth", "fertilizer.crop", "fertilizer.rotation",
    "fertilizer.application_method", "fertilizer.app1_date",
    "fertilizer.app2_date", "crops.cultivar", "crops.sow_crop",
    "crops.sowing_density", "crops.sowing_depth", "crops.row_spacing",
    "crops.harvest_crop"
  )

  continuousCols <- c(
    "tillage.residue", "fertilizer.app1_kg_n_ha", "fertilizer.app2_kg_n_ha",
    "crops.planting_date"
  )

  data.table::setkeyv(x, cols = discreteCols)

  # Crop planting date: from Dates to days since starting date.
  dateStart <- min(as.Date(x$crops.planting_date))
  x$crops.planting_date <- as.numeric(
    as.Date(x$crops.planting_date) - dateStart
  )

  # Sampling setup ----------------------------------------------------------
  nTotal               <- nrow(x)
  discreteCombinations <- unique(x[, ..discreteCols])
  nUniqueCombinations  <- nrow(discreteCombinations)
  nCombinations        <- x[, .N, by = discreteCols]$N

  if (nlhs < nUniqueCombinations)
    stop(
      sprintf(
        "Abort: LHS sample size < combinations of discrete variables (%d).",
        nUniqueCombinations
      )
    )

  # Draw samples ------------------------------------------------------------
  # 1. Draw a sample corresponding to each unique comb'n of discrete variables
  l <- lapply(1:nUniqueCombinations, function(i) {
    subdf    <- x[discreteCombinations[i, ..discreteCols], ..continuousCols]
    dtout    <- if (nrow(subdf) > 1) {
      nDraws   <- round(nCombinations[i] / nTotal * nlhs, 0)
      subdfLHS <- copula_lhs(
        x    = as.data.frame(subdf),
        n    = max(min(nDraws, nrow(subdf)), 0),
        seed = seed
      )
      colnames(subdfLHS) <- continuousCols
      data.table::data.table(subdfLHS)
    } else {
      subdf
    }

    cbind(discreteCombinations[i, ..discreteCols], dtout)
  })

  # 2. Bind and reconstruct crop planting dates
  out <- data.table::rbindlist(l)

  out$crops.planting_date <- as.character(dateStart + out$crops.planting_date)
  out$management_sample_id = 1:nrow(out)

  out
}

#' Title
#'
#' @param NResidue
#' @param NFerlitizer
#' @param NPlantingDate
#' @param m
#'
#' @return
#' @export
#'
#' @examples
management_population <- function(
  NResidue = NULL, NFerlitizer = NULL, NPlantingDate = NULL, m = 1) {

  if (is.null(NResidue))
    NResidue      <- length(seq(from = 0.4, to = 1.0, by = 0.1)) * m

  if (is.null(NFerlitizer))
    NFerlitizer   <- length(seq(from = 110, to = 230, by =  10)) * m

  if (is.null(NPlantingDate))
    NPlantingDate <- length(
      seq(from = as.Date("2020-04-20"), to = as.Date("2020-06-15"), by =  10)
    ) * m

  # Create non-random features ----------------------------------------------
  # 1. Make crop-rotation matrix
  rotCropScheme <- data.frame(
    rotation = c("maize-soybean", "maize-soybean", "cont-maize"),
    crop     = c("maize",         "soybean",       "maize"),
    stringsAsFactors = TRUE
  )

  # 2. Make tillage matrix
  tillageMatrix <- expand.grid(
    management = "tillage",
    implement  = "user_defined",
    rotation   = rotCropScheme$rotation,
    depth      = c(0, 127, 152, 179, 203),
    residue    = 1:NResidue,
    date       = as.Date("2020-11-10")
  )

  colnames(tillageMatrix) <- paste0("tillage.", colnames(tillageMatrix))

  # 3. Make fertilizer matrix
  applicationMatrix <- expand.grid(
    rotation    = rotCropScheme$rotation,
    crop        = rotCropScheme$crop,
    application = c("single", "split")
  )

  applicationMatrix <- unique( # Void impossibles
    applicationMatrix[!(
      applicationMatrix$crop == "soybean" &
        applicationMatrix$application == "split"
    ), ]
  )

  fertilizerMatrix  <- expand.grid(
    management  = "fertilizer",
    formulation = "uan_n",
    crop        = applicationMatrix$crop,
    rotation    = applicationMatrix$rotation,
    application = applicationMatrix$application,
    kg_n_na     = 1:NFerlitizer,
    app1_date   = 1, # Placeholder to be adjusted after drawing a planting date
    app2_date   = 1  # Idem
  )

  colnames(fertilizerMatrix) <- paste0("fertilizer.", colnames(fertilizerMatrix))

  # 4. Make crop matrix
  soyCultivarMatrix  <- expand.grid(
    cultivar       = c("MG_2", "MG_3"),
    sow_crop       = "soybean",
    planting_date  = 1:NPlantingDate,
    sowing_density = c(35, 37, 40, 42),
    sowing_depth   = 38,
    row_spacing    = c(381, 762),
    harvest_crop   = "soybean"
  )

  cornCultivarMatrix <- expand.grid(
    cultivar       = c("B_105", "B_110", "B_115"),
    sow_crop       = "maize",
    planting_date  = 1:NPlantingDate,
    sowing_density = 8,
    sowing_depth   = 51,
    row_spacing    = c(381, 508, 762),
    harvest_crop   = "maize"
  )

  cultivarMatrix    <- rbind(
    soyCultivarMatrix, cornCultivarMatrix
  )

  colnames(cultivarMatrix) <- paste0("crops.", colnames(cultivarMatrix))

  # 5. Make management matrix (combine all the previous)
  nrow(tillageMatrix)
  nrow(fertilizerMatrix)
  nrow(cultivarMatrix)

  expandedMatrix <- clean_expansion(expand.grid.df(
    clean_expansion(expand.grid.df(tillageMatrix, fertilizerMatrix)),
    cultivarMatrix
  ))

  # 6 Final tweaks
  expandedMatrix$fertilizer.application_method <-
    expandedMatrix$fertilizer.application


  # Add random features -----------------------------------------------------
  # 1. Residue
  #    0 if depth = 0,
  #    0.4-1.0 uniformly otherwise.
  expandedMatrix$tillage.residue <- 0

  ind <- expandedMatrix$tillage.depth > 0
  expandedMatrix$tillage.residue[ind] <-
    runif(n = sum(ind), min = 0.4, max = 1.0)

  # 2. Fertilizer application
  #    0 if crop = soybean
  #    110-230 (more weight on 140-200) if crop = maize, rot = maize-soybean
  #    170-280 (more weight on 200-250) if crop = maize, rot = cont-maize
  expandedMatrix$fertilizer.kg_n_ha <- 0

  ind <- expandedMatrix$fertilizer.crop == "maize" &
    expandedMatrix$fertilizer.rotation == "maize-soybean"
  expandedMatrix$fertilizer.kg_n_ha[ind] <-
    msm::rtnorm(
      n = sum(ind), mean = mean(c(140, 200)), sd = 20, lower = 110, upper = 230
    )

  ind <- expandedMatrix$fertilizer.crop == "maize" &
    expandedMatrix$fertilizer.rotation == "cont-maize"
  expandedMatrix$fertilizer.kg_n_ha[ind] <-
    msm::rtnorm(
      n = sum(ind), mean = mean(c(200, 250)), sd = 20, lower = 170, upper = 280
    )

  # 3. Planting date
  #    20-apr to 15-jun uniformly.
  dateFrom <- as.Date("2020-04-20")
  dateTo   <- as.Date("2020-06-15")
  expandedMatrix$crops.planting_date <-
    sample(
      x       = seq(dateFrom, dateTo, by = 1),
      size    = nrow(expandedMatrix),
      replace = TRUE
    )

  # 4. Final tweaks: fertilizers
  #    If application = single, app1 = 100% (planting date), app2 =  0% (+25 ds)
  #    If application = split,  app1 =  40% (planting date), app2 = 60% (+25 ds)
  expandedMatrix$fertilizer.app1_kg_n_ha <- expandedMatrix$fertilizer.kg_n_ha
  expandedMatrix$fertilizer.app2_kg_n_ha <- 0
  expandedMatrix$fertilizer.app1_date    <- expandedMatrix$crops.planting_date
  expandedMatrix$fertilizer.app2_date    <- expandedMatrix$crops.planting_date

  ind <- expandedMatrix$fertilizer.application_method == "split"
  expandedMatrix$fertilizer.app1_kg_n_ha[ind] <-
    0.4 * expandedMatrix$fertilizer.kg_n_ha[ind]
  expandedMatrix$fertilizer.app2_kg_n_ha[ind] <-
    0.6 * expandedMatrix$fertilizer.kg_n_ha[ind]
  expandedMatrix$fertilizer.app2_date[ind] <-
    expandedMatrix$crops.planting_date[ind] + 25

  # Prepare output ----------------------------------------------------------
  # column inclusion and names as designed by Matt
  colsKeep <- c(
    "tillage.implement", "tillage.depth", "tillage.residue", "tillage.date",
    "fertilizer.formulation", "fertilizer.crop", "fertilizer.rotation",
    "fertilizer.application_method", "fertilizer.app1_kg_n_ha",
    "fertilizer.app1_date", "fertilizer.app2_kg_n_ha", "fertilizer.app2_date",
    "crops.cultivar", "crops.sow_crop", "crops.planting_date",
    "crops.sowing_density", "crops.sowing_depth", "crops.row_spacing",
    "crops.harvest_crop"
  )

  colsKeep[!(colsKeep %in% colnames(expandedMatrix))]

  unique(expandedMatrix[, colsKeep])
}

clean_expansion <- function(x) {
  if (all(c("fertilizer.crop", "crops.sow_crop") %in% colnames(x)))
    x <- x[x$fertilizer.crop == x$crops.sow_crop, ]

  if (all(c("tillage.rotation", "fertilizer.crop") %in% colnames(x)))
    x <- x[!(x$tillage.rotation == "cont-maize" & x$fertilizer.crop == "soybean"), ]

  if (all(c("tillage.rotation", "fertilizer.rotation") %in% colnames(x)))
    x <- x[x$tillage.rotation == x$fertilizer.rotation, ]

  rownames(x) <- NULL

  unique(x)
}

# https://stackoverflow.com/a/21911221/2860744
expand.grid.df <- function(...) Reduce(function(...) merge(..., by = NULL), list(...))