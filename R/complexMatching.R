## function to match features. columns must be provided by position

#' Match features based on m/z using only more than one adduct (vectorized
#' version of \code{matchFeaturesSimple})
#'
#' Function intended to match features to metabolites based on m/z. This is a
#' vectorized version of \code{matchFeaturesSimple} intented to consider more
#' than one adduct.
#'
#' @param data \code{data.frame} with features to be matched. It must have the
#'   following columns: \itemize{\item \code{col.data.mz}: experimental m/z
#'   (with ions). \item \code{col.data.name}: column with unique feature names.
#'   \item Extra columns can be included and will be ignored.} column with
#'   experimental m/z (with adduct) and a column with unique feature names.
#' @param reference \code{data.frame} with metabolites to be used as reference.
#'   It must have the following columns: \itemize{ \item \code{cols.ref.mz}:
#'   charged masses (m/z) with the corresponding ions to be used \item
#'   \code{col.ref.common.name}: common name of metabolites (optional). \item
#'   \code{col.ref.systematic.name}: systematic name of metabolites (optional).
#'   \item \code{col.ref.neutral.mass}: neutral mass of metabolites (optional).
#'   \item \code{col.ref.type}: class of metabolite (optional).}
#' @param error Error to be used. It is calculated as p.p.m
#' @param col.data.mz Column name in \code{data} with experimental m/z (with
#'   ions) of features.
#' @param col.data.name Column name in \code{data} with unique feature names.
#' @param cols.ref.mz Column names in \code{reference} with charged masses (m/z)
#'   with the corresponding ions to be used.
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
matchFeaturesComplex <- function(
  data,
  reference,
  col.data.mz,
  col.data.name,
  cols.ref.mz,
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
  reference <- reference[!is.na(reference[cols.ref.mz]), ]
  # matching features with metabolites using m/z (vectorized)
  listRes <- lapply(
    X = cols.ref.mz,
    FUN = function(x) {
      .matchErrPpms(
        col.ref.mz = x,
        data = data,
        reference = reference,
        col.data.mz = col.data.mz,
        error = error
      )
    }
  )
  filfListRes <- sapply(X = listRes, FUN = is.null)
  listRes <- listRes[!filfListRes]
  namesListRes <- cols.ref.mz[!filfListRes]
  ## build final matching vectorized
  resFinalList <- lapply(
    X = listRes,
    FUN = function(results) {
      finalMatching <- .buildFinalMatching(
        data = data,
        results = results,
        col.data.name = col.data.name,
        col.ref.common.name = col.ref.common.name,
        col.ref.systematic.name = col.ref.systematic.name,
        col.ref.type = col.ref.type,
        col.ref.neutral.mass = col.ref.neutral.mass
      )
      return(finalMatching)
    }
  )
  names(resFinalList) <- namesListRes
  return(resFinalList)
}


.matchErrPpms <- function(
  data,
  reference,
  col.ref.mz,
  col.data.mz,
  error
) {
  # matching features with metabolites using m/z
  results <- apply(
    X = data,
    MARGIN = 1,
    FUN = function(x, vec, err) {
      ppmErr <- abs((as.numeric(x[col.data.mz]) - vec) / vec * 10^6)[, 1]
      # .ppmErrCalc(x = x, ref = vec, col.data.mz = col.data.mz)
      # abs((as.numeric(x[col.data.mz]) - vec) / vec * 10^6)[, 1]
      filfNA <- !is.na(ppmErr)
      filfErr <- ppmErr <= err
      filtered <- reference[filfNA & filfErr, ]
      if (nrow(filtered) != 0) {
        filtered$ppm_error <- ppmErr[filfNA & filfErr]
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
  return(results)
}

