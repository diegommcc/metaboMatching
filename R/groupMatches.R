#' @importFrom openxlsx loadWorkbook createWorkbook addWorksheet writeData saveWorkbook
NULL


#' Group matches based on m/z
#'
#' Similar matches are grouped by introducing a row of lines. This function is
#' just with exploratory proposes.
#'
#' @param data \code{data.frame} with m/z values to be used.
#' @param col.data.mz Name or position of m/z values in \code{data}.
#' @param error Number of considered decimals to group similar features.
#'
#' @return \code{data.frame} provided in \code{data} with extra rows separating
#'   similar matches.
#'
#' @export
#'
groupByMz <- function(
  data,
  col.data.mz,
  error = 2
) {
  dataFilfSorted <- data[order(data[, col.data.mz], decreasing = TRUE), ]
  vec <- round(as.numeric(dataFilfSorted[, col.data.mz]), error)
  listRes <- lapply(
    X = as.list(unique(vec)),
    FUN = function(x) {
      true <- vec == x
      if (length(vec) != max(which(true)) || max(which(true)) != 1)
        return(rbind(dataFilfSorted[true, ], "-"))
      else
        return(dataFilfSorted[true, ])
    }
  )
  return(do.call(rbind, listRes))
}

#' Group matches based on m/z and explort results as a xlsx file
#'
#' Alternative version of \code{\link{groupByMz}}. It generates a xlsx file with
#' similar features grouped by squares. This function uses the \pkg{openxlsx}
#' package, and it is not able to handle big matrices. Don't try to use it with
#' matrices larger than 500 rows.
#'
#' @param data \code{data.frame} with m/z values to be used.
#' @param col.data.mz Name or position of m/z values in \code{data}.
#' @param path.xlsx Path in which xlsx file will be created.
#' @param error Number of considered decimals to group similar features.
#' @param return.wb Boolean indicating if Workbook is returned. \code{FALSE} by
#'   default.
#' @param border.colour Color of borders (\code{'black'} by default).
#' @param border.style Style of borders (\code{'medium'} by default). See
#'   \code{openxlsx::writeData} to check available options.
#'
#' @return If \code{return.wb} is FALSE, nothing will be returned. Otherwise, a
#'   Workbook object is returned.
#'
#' @export
#'
groupByMzToXLSX <- function(
  data,
  col.data.mz,
  path.xlsx,
  add.sheet = FALSE,
  name.sheet = NULL,
  error = 2,
  return.wb = FALSE,
  border.colour = "black",
  border.style = "medium"
) {
  dataFilfSorted <- data[order(data[, col.data.mz], decreasing = TRUE), ]
  vec <- round(as.numeric(dataFilfSorted[, col.data.mz]), error)
  ## add new sheets to a pre-existing file
  if (add.sheet) {
    wb <- loadWorkbook(path.xlsx)
  } else {
    wb <- createWorkbook()
  }
  if (is.null(name.sheet)) name.sheet <- "Grouped"
  addWorksheet(wb, name.sheet)

  for (i in unique(vec)) {
    true <- vec == i
    if (1 %in% which(true)) {
      writeData(
        wb = wb, sheet = name.sheet, x = dataFilfSorted[true, ],
        startCol = 1, startRow = which(true), array = TRUE,
        borders = "surrounding", colNames = TRUE,
        borderStyle = border.style, borderColour = border.colour
      )
    } else {
      writeData(
        wb = wb, sheet = name.sheet, x = dataFilfSorted[true, ],
        startCol = 1, startRow = which(true)[1], array = TRUE,
        borders = "surrounding", colNames = FALSE,
        borderStyle = border.style, borderColour = border.colour
      )
    }
  }
  saveWorkbook(
    wb, file = path.xlsx,
    overwrite = TRUE, returnValue = FALSE
  )
  ## if return wb
  if (return.wb) return(wb)
}


## the problem with this function is when borders are created. I think it is related to colNames
