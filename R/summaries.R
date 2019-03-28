### summaries.R: print and summary functions


#' Print method for ffmanova
#' 
#' Print method for objects of class \code{"ffmanova"}.  It prints an ANOVA
#' table.
#' 
#' The function constructs an anova table, and prints it using
#' \code{\link{printCoefmat}} with tailored arguments.
#' 
#' @param x \code{"ffmanova"} object.  Typically created by
#' \code{\link{ffmanova}}.
#' @param digits positive integer.  Minimum number of significant digits to be
#' used for printing most numbers.
#' @param \dots further arguments sent to the underlying
#' \code{\link{printCoefmat}}.
#' @return Invisibly returns the original object.
#' @author Bjørn-Helge Mevik
#' @seealso \code{\link{ffmanova}}, \code{\link{printCoefmat}}
#' @keywords print internal
#' @importFrom stats printCoefmat
#' @export
print.ffmanova <- function(x, digits = max(getOption("digits") - 3, 3), ...) {
    cat("--- 50-50 MANOVA ",
        sum(x$df), " objects -- ",
        ncol(x$pRaw), " responses",
        if (x$stand) " (standardised)",
        ":\n", sep = "")
    tab <- with(x, data.frame(df, exVarSS, c(nPC, NA), c(nBU, NA),
                              c(exVarPC, NA), c(exVarBU, NA), c(pValues, NA)))
    dimnames(tab) <- list(c(x$termNames, "Residuals"),
                          c("Df", "exVarSS", "nPC", "nBU", "exVarPC",
                            "exVarBU", "p-Value"))
    if(x$termNames[1]=="(Intercept)")
      tab <- tab[-1,]                     # Drop the (Intercept) row
    printCoefmat(tab, digits = digits, cs.ind = 2, tst.ind = NULL,
                 zap.ind = c(1,3,4), has.Pvalue = TRUE, na.print = "", ...)
    invisible(x)
}
