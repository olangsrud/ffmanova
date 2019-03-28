## ffmanova 1.1.0

* New function `predict.ffmanova`. Thus, the generic function, predict, can take output from ffmanova as input. 
  - Description: _The same predictions as lm can be obtained. With some variables missing in input, adjusted means or mean predictions are computed (Langsrud et al., 2007). Linear combinations of such predictions, with standard errors, can also be obtained._
* New function `ffAnova`. 
  - Description: _Analysis of variance table for a linear model using type II&ast; sums of squares, which are described in Langsrud et al. (2007). Type II&ast; extends the type II philosophy to continuous variables. The results are invariant to scale changes and pitfalls are avoided._
* Models without a constant term is supported.
* The function ffmanova has been extended with new input parameters and additional elements in output is possible. These changes are related to the above new functionality. By calling ffmanova with returnModel=FALSE, output should be the same as in version 1.0.0. 


## ffmanova 1.0.0

* In January 2019, ffmanova got version number 1.0.0.  Since the first realise on CRAN in 2006, ffmanova has been a stable working horse. Cosmetic changes and changes to meet new package standards have been made. The source code for the main computations is close to a direct translation of an implementation in matlab. 
  - Description of `ffmanova`: _General linear modeling of fixed-effects models with multiple responses is performed. The function calculates 50-50 MANOVA p-values, ordinary univariate p-values and adjusted p-values using rotation testing._
