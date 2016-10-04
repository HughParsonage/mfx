#' Odds ratios for a beta regression.
#' 
#' This function estimates a beta regression model and calculates the
#' corresponding odds ratios.
#' 
#' The underlying link function in the mean model (mu) is "logit". If both
#' \code{robust=TRUE} and \code{!is.null(clustervar1)} the function overrides
#' the \code{robust} command and computes clustered standard errors.
#' 
#' @aliases betaor print.betaor
#' @param formula an object of class ``formula'' (or one that can be coerced to
#' that class).
#' @param data the data frame containing these data. This argument must be
#' used.
#' @param robust if \code{TRUE} the function reports White/robust standard
#' errors.
#' @param clustervar1 a character value naming the first cluster on which to
#' adjust the standard errors.
#' @param clustervar2 a character value naming the second cluster on which to
#' adjust the standard errors for two-way clustering.
#' @param control a list of control arguments specified via
#' \code{\link[betareg]{betareg.control}}.
#' @param link.phi as in the \code{\link[betareg]{betareg}} function.
#' @param type as in the \code{\link[betareg]{betareg}} function.
#' @return \item{oddsratio}{a coefficient matrix with columns containing the
#' estimates, associated standard errors, test statistics and p-values.}
#' \item{fit}{the fitted \code{\link[betareg]{betareg}} object.}
#' \item{call}{the matched call.}
#' @seealso \code{\link{betamfx}}, \code{\link[betareg]{betareg}}
#' @references Francisco Cribari-Neto, Achim Zeileis (2010). Beta Regression in
#' R. Journal of Statistical Software 34(2), 1-24.
#' 
#' Bettina Gruen, Ioannis Kosmidis, Achim Zeileis (2012). Extended Beta
#' Regression in R: Shaken, Stirred, Mixed, and Partitioned. Journal of
#' Statistical Software, 48(11), 1-25.
#' @examples
#' 
#' # simulate some data
#' set.seed(12345)
#' n = 1000
#' x = rnorm(n)
#' 
#' # beta outcome
#' y = rbeta(n, shape1 = plogis(1 + 0.5 * x), shape2 = (abs(0.2*x)))
#' # use Smithson and Verkuilen correction
#' y = (y*(n-1)+0.5)/n
#' 
#' data = data.frame(y,x)
#' betaor(y~x|x, data=data)
#' 
#' @export 
betaor <-
function(formula, data, robust = FALSE, clustervar1 = NULL, 
                   clustervar2 = NULL, control = betareg.control(), 
                   link.phi = NULL, type = "ML"){
  
  res = betaorest(formula, data, robust, clustervar1, clustervar2, 
                  control, link.phi, type)
  
  est = NULL
  zstat = log(res$or$oddsratio)*res$or$oddsratio/res$or$se
  
  est$oddsratio = cbind(OddsRatio = res$or$oddsratio,
                        StdErr = res$or$se,
                        z.value = zstat,
                        p.value = 2*pt(-abs(zstat), df = Inf))
  
  colnames(est$oddsratio) = c("OddsRatio","Std. Err.","z","P>|z|")
  rownames(est$oddsratio) =  rownames(res$or)
  est$fit = res$fit
  est$call = match.call() 
  class(est) = "betaor"
  est
}
