#' @importFrom dplyr %>%
NULL

## simple atching --> using just one adduct

## function to match features. columns must be provided by position

#' Match features based on m/z using only one adduct
#'
#' Function intended to match features to metabolites based on m/z. In this
#' case, only one adduct is taken into account. For a vectorized version, see
#' the \code{matchFeaturesComplex} function.
#'
#' @param data \code{data.frame} with features to be matched. It must have the
#'   following columns: \itemize{\item \code{col.data.mz}: experimental m/z
#'   (with ions). \item \code{col.data.name}: column with unique feature names.
#'   \item Extra columns can be included and will be ignored.} column with
#'   experimental m/z (with adduct) and a column with unique feature names.
#' @param reference \code{data.frame} with metabolites to be used as reference.
#'   It must have the following columns: \itemize{ \item \code{col.ref.mz}:
#'   charged mass (m/z) with the corresponding ion to be used \item
#'   \code{col.ref.common.name}: common name of metabolites (optional). \item
#'   \code{col.ref.systematic.name}: systematic name of metabolites (optional).
#'   \item \code{col.ref.neutral.mass}: neutral mass of metabolites (optional).
#'   \item \code{col.ref.type}: class of metabolite (optional).}
#' @param error Error to be used. It is calculated as p.p.m
#' @param col.data.mz Column name in \code{data} with experimental m/z (with
#'   ions) of features.
#' @param col.data.name Column name in \code{data} with unique feature names.
#' @param col.ref.mz Column name in \code{reference} with charged mass (m/z)
#'   with the corresponding ion to be used.
#' @param col.ref.common.name Column name in \code{reference} with common name
#'   of metabolites (optional).
#' @param col.ref.systematic.name Column name in \code{reference} with
#'   systematic name of metabolites (optional).
#' @param col.ref.neutral.mass Column name in \code{reference} with neutral
#'   class of metabolites (optional).
#' @param col.ref.type Column name in \code{reference} with class of of
#'   metabolites (optional).
#'
#' @return A \code{list} with two elements: a filtered version of the original
#'   \code{data.frame} with only found features and a full version of the
#'   original \code{data.frame} with all features. The information associated
#'   with matching is aggregated as new columns in both cases.
#'
#' @export
#'
#' @seealso \code{matchFeaturesComplex}
#'
#'
matchFeaturesSimple <- function(
  data,
  reference,
  col.data.mz,
  col.data.name,
  col.ref.mz,
  col.ref.common.name = "COMMON",
  col.ref.systematic.name = "SYSTEMATIC_NAME",
  col.ref.neutral.mass = NULL,
  col.ref.type = "TYPE",
  error = 5
) {
  if (is.character(col.data.mz)) {
    col.data.mz <- which(colnames(data) == col.data.mz)
  }
  if (is.character(col.data.name)) {
    col.data.name <- which(colnames(data) == col.data.name)
  }
  # to deal with NAs in the database, we remove all the NAs for the selected col
  reference <- reference[!is.na(reference[col.ref.mz]), ]
  # matching features with metabolites using m/z
  results <- apply(
    X = data,
    MARGIN = 1,
    FUN = function(x, vec, err) {
      ppmErr <- abs((as.numeric(x[col.data.mz]) - vec) / vec * 10^6)[, 1]
      filf <- ppmErr <= err
      filtered <- reference[filf, ]
      if (nrow(filtered) != 0) {
        filtered$ppm_error <- ppmErr[filf]
        l <- filtered # list(filtered)
        # names(l) <- x[col.data.name]
        return(l)
      }
    },
    vec = reference[col.ref.mz],
    err = error
  )
  # removing NULL
  results[sapply(X = results, FUN = is.null)] <- NULL
  # results[sapply(X = results, FUN = function(x) any(is.na(x)))] <- NULL
  resList <- .buildFinalMatching(
    data = data,
    results = results,
    col.data.name = col.data.name,
    col.ref.common.name = col.ref.common.name,
    col.ref.systematic.name = col.ref.systematic.name,
    col.ref.type = col.ref.type,
    col.ref.neutral.mass = col.ref.neutral.mass
  )
  return(resList)
}


## function to summarized results into a list of two elements:
# full matrix with found features marked
# filtered matrix with only found features
.buildFinalMatching <- function(
  data,
  results,
  col.data.name,
  col.ref.common.name,
  col.ref.systematic.name,
  col.ref.type,
  col.ref.neutral.mass
) {
  # filtered names
  filfNames <- names(results)
  filfComp <- lapply(
    X = results,
    FUN = function(x) {
      x[order(x[["ppm_error"]]),
        c(col.ref.common.name, col.ref.type,
          col.ref.systematic.name,
          col.ref.neutral.mass, "ppm_error")]
    }
  )
  filfCompMin <- lapply(X = filfComp, FUN = function(x) x[1, ])
  ## concatenate name, mass and error for multiple matches
  typeComp <- sapply(
    X = filfComp[filfNames],
    FUN = function(x) {
      paste(
        sapply(
          X = seq_along(x[, col.ref.common.name]),
          FUN = function(y) {
            paste0(
              x[, col.ref.common.name][y], " (",
              x[y, col.ref.type], " - ",
              paste(
                round(x[y, -which(names(x[y, ]) %in% c(col.ref.common.name,
                                                       col.ref.systematic.name,
                                                       col.ref.type))], 3),
                collapse = " - "
              ), ")", collapse = " - ")
          }
        ), collapse = " | "
      )
    }
  )
  ## create data filtered
  dataFilf <- data %>% filter(data[, col.data.name] %in% filfNames)
  dataFilf[dataFilf == ""] <- NA
  # filtered compounds
  dataFilf <- dataFilf %>% mutate(
    Compound.MinErr = unlist(
      sapply(
        X = filfCompMin[filfNames],
        FUN = function(x) {
          if (any(dim(x[, col.ref.common.name]) == 0)) NULL
          else x[, col.ref.common.name]
        }
      )
    ),
    Type.MinErr = unlist(
      sapply(
        X = filfCompMin[filfNames],
        FUN = function(x) {
          if (any(dim(x[, col.ref.type]) == 0)) NULL
          else x[, col.ref.type]
        }
      )
    ),
    PPM.MinErr = unlist(
      sapply(
        X = filfCompMin[filfNames],
        FUN = function(x) {
          if (any(dim(x[, "ppm_error"]) == 0)) NULL
          else x[, "ppm_error"]
        }
      )
    ),
    NeutralMass.MinErr = unlist(
      sapply(
        X = filfCompMin[filfNames],
        FUN = function(x) {
          if (any(dim(x[, col.ref.neutral.mass]) == 0)) NULL
          else x[, col.ref.neutral.mass]
        }
      )
    ),
    SystematicName.MinErr = unlist(
      sapply(
        X = filfCompMin[filfNames],
        FUN = function(x) {
          if (any(dim(x[, col.ref.systematic.name]) == 0)) NULL
          else x[, col.ref.systematic.name]
        }
      )
    ),
    AllCompounds = typeComp,
    Type = sapply(
      X = filfComp[filfNames],
      FUN = function(x) {
        paste0(x[, col.ref.type], collapse = " | ")
      }
    )
  )
  ## create data completed with found features marked
  dataMatrixMark <- data
  dataMatrixMark$Mark <- data[, col.data.name] %in% filfNames
  dataMatrixMark[dataMatrixMark == ""] <- NA
  ## type
  dataMatrixMark$Type <- NA
  dataMatrixMark$Type[data[, col.data.name] %in% typeComp] <- typeComp
  ## compounds
  dataMatrixMark$Compound <- "-"
  dataMatrixMark$Compound[data[, col.data.name] %in% filfNames] <- filfComp

  return(list(dataFilf, dataMatrixMark))
}

.ppmErrCalc <- function(x, ref, col.data.mz) {
  return(abs((as.numeric(x[col.data.mz]) - ref) / ref * 10^6)[, 1])
}
