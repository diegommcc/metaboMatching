

## calculate adducts for reference

#' Calculate adducts for reference metabolites
#'
#' It calculates adducts for reference metabolites and aggregates information to
#' the original \code{data.frame}. Importantly, it uses prior information about
#' which adducts are more likely to be present according metabolite subclass.
#' !!!!!!!This documentation is not completed.!!!!!!!!!
#'
#' @param reference \code{data.frame} with metabolites to be used as reference.
#'   It must have the following columns: \itemize{ \item \code{cols.ref.mz}:
#'   charged masses (m/z) with the corresponding ions to be used \item
#'   \code{col.ref.common.name}: common name of metabolites (optional). \item
#'   \code{col.ref.systematic.name}: systematic name of metabolites (optional).
#'   \item \code{col.ref.neutral.mass}: neutral mass of metabolites (optional).
#'   \item \code{col.ref.type}: class of metabolite (optional).}
#' @param guide.adds Error to be used. It is calculated as p.p.m
#' @param col.ionization.guide Column name in \code{data} with experimental m/z (with
#'   ions) of features.
#' @param col.mz.ref Column name in \code{data} with unique feature names.
#' @param col.subclass.ref Column names in \code{reference} with charged masses (m/z)
#'   with the corresponding ions to be used.
#' @param mass.ions Column name in \code{reference} with common name
#'   of metabolites (optional).
#'
#' @return A \code{list} with two elements: a filtered version of the original
#'   \code{data.frame} with only found features and a full version of the
#'   original \code{data.frame} with all features. The information associated
#'   with matching is aggregated as new columns in both cases.
#'
#' @export
#'
calculateAdducts <- function(
  reference,
  guide.adds,
  col.ionization.guide,
  col.mz.ref,
  col.subclass.ref,
  mass.ions
) {
  for (ioni.sel in unique(guide.adds[[col.ionization.guide]])) {
    guide.filf <- guide.adds %>% filter(.data[[col.ionization.guide]] == ioni.sel)
    guide.filf <- guide.filf[, -which(colnames(guide.filf) == col.ionization.guide)]
    reference <- .addAdductsRef(
      guide.adds = guide.filf,
      reference = reference,
      sub.class.ref = col.subclass.ref,
      mz.ref = col.mz.ref,
      masses.ions = mass.ions
    )
  }
  return(reference)
}

.formatNames <- function(x) {
  x1 <- substring(x, first = 4)
  namesIons <- gsub(pattern = "\\][-|+]", replacement = "", x = x1)
  ## operation to do
  op <- substring(x, first = 3, last = 3)
  return(list(namesIons, op))
}

.addAdductsRef <- function(
  guide.adds,
  reference,
  sub.class.ref,
  mz.ref,
  masses.ions
) {
  operation <- data.frame(
    Ion = .formatNames(rownames(guide.adds))[[1]],
    Operation = .formatNames(rownames(guide.adds))[[2]]
  )
  selRef <- sapply(
    X = colnames(guide.adds),
    FUN = function(x) reference[[sub.class.ref]] == x
  )
  ## order columns just in case
  selRef <- selRef[, colnames(guide.adds)]
  ress <- lapply(
    X = seq(ncol(selRef)),
    FUN = function(ncol) {
      sel <- reference[selRef[, ncol], ]
      toDo <- guide.adds[, ncol]
      ## if there are not metabolites for this ionization mode
      if (!any(toDo)) return(sel)
      toDoG <- operation[toDo, ]
      for (i in seq(nrow(toDoG))) {
        if (toDoG$Operation[i] == "-") {
          sel[[paste0("M", toDoG$Operation[i], toDoG$Ion[i])]] <-
            sel[[mz.ref]] - masses.ions[[toDoG$Ion[i]]]
        } else {
          sel[[paste0("M", toDoG$Operation[i], toDoG$Ion[i])]] <-
            sel[[mz.ref]] + masses.ions[[toDoG$Ion[i]]]
        }
      }
      return(sel)
    }
  )
  return(do.call(what = plyr::rbind.fill, args = ress))
}
