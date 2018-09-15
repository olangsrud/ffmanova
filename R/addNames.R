
# New function in 2018 used to add names to the output from manova5050(), rotationtests() and unitests() 
# and thus all output from ffmanova().
addNames <- function(x, rowNames, colNames = NULL, residNames = "Residuals") {
  allTermNames <- c(rowNames, residNames)
  for (i in seq_len(length(x))) {
    if (is.matrix(x[[i]])) {
      rownames(x[[i]]) <- rep_len(allTermNames, NROW(x[[i]]))
      if (!is.null(colNames)) {
        colnames(x[[i]]) <- rep_len(colNames, NCOL(x[[i]]))
      }
    } else {
      if (!is.logical(x[[i]])) 
        names(x[[i]]) <- rep_len(allTermNames, length(x[[i]]))
    }
  }
  x
}