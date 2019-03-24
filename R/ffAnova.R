#' Type II* Anova
#' 
#' Analysis of variance table for a linear model using type \code{II}* sums of squares,
#' which are described in Langsrud et al. (2007).
#' Type \code{II}* extends the type \code{II} philosophy to continuous variables. 
#' The results are invariant to scale changes and pitfalls are avoided. 
#' 
#' This function is a variant of \code{\link{ffmanova}} for the univariate special case.
#' The two input parameters will be interpreted by \code{\link{model.frame}}.
#'
#' @param formula A model \code{\link{formula}} or an \code{R} object 
#'     (preferably an \code{\link{lm}} object).
#' @param data An optional data frame, list or environment.
#'
#' @return An object of class \code{"anova"} (see \code{\link{anova}}).
#' @export
#' @author Øyvind Langsrud and Bjørn-Helge Mevik
#' @references 
#' Langsrud, Ø., Jørgensen, K., Ofstad, R. and Næs, T. (2007):
#'  \dQuote{Analyzing Designed Experiments with Multiple Responses},
#'  \emph{Journal of Applied Statistics}, \bold{34}, 1275-1296.
#'
#' @examples
#' # Generate example data
#' set.seed(123)
#' a <- c(0, 0, 0, 10, 10, 10, 1, 1, 1)
#' A <- as.character(a)  # A is categorical
#' b <- 1:9
#' y <- rnorm(9)/10 + a  # y depends strongly on a (and A)
#' a100 <- a + 100  # change of scale (origin)
#' b100 <- b + 100  # change of scale (origin)
#' 
#' # Four ways of obtaining the same results
#' ffAnova(y ~ A * b)
#' ffAnova(y ~ A * b100)
#' ffAnova(lm(y ~ A * b))
#' ffAnova(y ~ A * b, data.frame(A = A, y = y, b = 1:9))
#' 
#' # Second order continuous variable
#' ffAnova(y ~ a + I(a^2))
#' 
#' # Model equivalent to 'y ~ A * b'
#' ffAnova(y ~ (a + I(a^2)) * b)
#' 
#' # Demonstrating similarities and differences using package car
#' if (!require(car))        # Package car is loaded if available 
#'   Anova <- function(x) {  # Replacement function if car not available
#'     warning("No results since package car is not available")}
#' 
#' lm_Ab <- lm(y ~ A * b)
#' lm_Ab100 <- lm(y ~ A * b100)
#' 
#' # Type II same as type II* in this case
#' Anova(lm_Ab)      # Type II
#' Anova(lm_Ab100)   # Type II
#' ffAnova(lm_Ab)    # Type II*
#' ffAnova(lm_Ab100) # Type II*
#' 
#' # Type III depends on scale
#' Anova(lm_Ab, type = 3)
#' Anova(lm_Ab100, type = 3)
#' 
#' lm_a <- lm(y ~ a + I(a^2))
#' lm_a100 <- lm(y ~ a100 + I(a100^2))
#' 
#' # Now Type II depends on scale
#' Anova(lm_a)      # Type II
#' Anova(lm_a100)   # Type II
#' ffAnova(lm_a)    # Type II*
#' ffAnova(lm_a100) # Type II*
ffAnova <- function(formula, data = NULL) {
  if (class(formula)[1] == "ffmanova") {
    return(ffmanova2anova(formula))
  }
  # ffmanova(formula, data, outputClass = 'anova') # problem when lm imput
  sysCall <- sys.call()
  sysCall[[1]] <- as.name("ffmanova")
  sysCall$outputClass <- "anova"
  parentFrame <- parent.frame()
  return(eval(sysCall, envir = parentFrame))
}



ffmanova2anova <- function(res) {
  if (length(res$ffModel$msError) != 1) {
    stop("Only a single response variable allowed.")
  }
  ssTot <- res$ffModel$ssTot
  if (!is.null(res$ffModel$scaleY)) {
    ssTot <- ssTot * res$ffModel$scaleY^2
  }
  tab <- with(res, data.frame(df, exVarSS * ssTot, exVarSS * ssTot/df, c(stat, NA), c(pValues, NA)))
  dimnames(tab) <- list(c(res$termNames, "Residuals"), c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
  if (res$termNames[1] == "(Intercept)") 
    tab <- tab[-1, ]
  tstat <- tab[, 1] == 1 & !is.na(tab[, 4])
  tab[tstat, 4] <- tab[tstat, 4]^2
  return(structure(tab, class = c("anova", "data.frame"), 
                   heading = c("Anova Table (Type II* tests)\n", paste("Response:", res$ffModel$colnamesY))))
}



ffAnova_Version_1_0_1 <- function(formula, data = NULL) {
  
  # Code copied from ffmanova
  # Duplicate code here to avoid interfering with old/stable code

  ## Get the model frame.  META: This is unneccessary general for the
  ## moment, but perhaps subset and na.action will be added later.
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]                # Retain only the named arguments
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## Get the terms
  mt <- attr(mf, "terms")
  
  ## Get the data matrices:
  mm <- model.matrix(mt, mf)
  Y <- as.matrix(model.response(mf, "numeric"))
  
  #-- code in ffmanova but not ffAnova--  if (stand) Y <- stdize(Y, center = FALSE, avoid.zero.divisor = TRUE)
  
  ### New code for ffAnova
  if (ncol(Y) != 1) {
    stop("Only a single response variable allowed.")
  } else {                           # colnames(Y) <- deparse(formula[[2]]) not working when lm object input
    if (max(abs(Y - mf[, 1])) == 0)  #  unnecessary check? # much faster than identical(as.vector(Y),as.vector(mf[,1]))
      colnames(Y) <- colnames(mf)[1]
  }  ### End 'New code for ffAnova'
  
  ## Create a `fator/term index matrix':
  mOld <- attr(mt, "factors")
  ## Fix any I() terms:
  mNew <- fixModelMatrix(mOld)
  ## add constant term
  mNew <- cbind(`(Intercept)` = 0, mNew)
  ## transpose
  model <- t(mNew)
  
  ## Split the model matrix into matrices for each term:
  termNr <- attr(mm, "assign") + 1
  D <- vector("list", max(termNr))
  for (i in seq(along = D)) 
    D[[i]] <- mm[, termNr == i, drop = FALSE]
  
  xObj <- x_Obj(D, model)
  xyObj <- xy_Obj(xObj, Y)
  
  nTerms <- length(xyObj$xObj$df_D_test)
  
  
  ### New code for ffAnova
  res <- c(manova5050(xyObj, FALSE), unitests(xyObj))
  ssTot <- sum((Y - mean(Y))^2)
  
  tab <- with(res, data.frame(df, exVarSS * ssTot, exVarSS * ssTot/df, c(stat, NA), c(pValues, NA)))
  dimnames(tab) <- list(c(res$termNames, "Residuals"), c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
  tab <- tab[-1, ]
  tstat <- tab[, 1] == 1 & !is.na(tab[, 4])
  tab[tstat, 4] <- tab[tstat, 4]^2
  structure(tab, 
            class = c("anova", "data.frame"), 
            heading = c("Anova Table (Type II* tests)\n", paste("Response:", colnames(Y))))
}

