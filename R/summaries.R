### summaries.R: print and summary functions
### $Id$

print.ffmanova <- function(x, digits = max(getOption("digits") - 3, 3), ...) {
    cat("--- 50-50 MANOVA Version R1.0 ---",
        sum(x$df), "(ok?) objects --",
        ncol(x$pRaw), "(ok?) responses:\n")
    tab <- with(x, data.frame(df, exVarSS, c(nPC, NA), c(nBU, NA),
                              c(exVarPC, NA), c(exVarBU, NA), c(pValues, NA)))
    dimnames(tab) <- list(c(x$termNames, "Residuals"),
                          c("Df", "exVarSS", "nPC", "nBU", "exVarPC",
                            "exVarBU", "p-Value"))
    tab <- tab[-1,]                     # Drop the (Intercept) row
    printCoefmat(tab, digits = digits, cs.ind = 2, tst.ind = NULL,
                 zap.ind = c(1,3,4), has.Pvalue = TRUE, na.print = "", ...)
    invisible(tab)
}
