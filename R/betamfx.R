#' Marginal effects for a beta regression.
#' 
#' This function estimates a beta regression model and calculates the
#' corresponding marginal effects.
#' 
#' The underlying link function in the mean model (mu) is ``logit''. If both
#' \code{robust=TRUE} and \code{!is.null(clustervar1)} the function overrides
#' the \code{robust} command and computes clustered standard errors.
#' 
#' @aliases betamfx print.betamfx
#' @param formula an object of class ``formula'' (or one that can be coerced to
#' that class).
#' @param data the data frame containing these data. This argument must be
#' used.
#' @param atmean default marginal effects represent the partial effects for the
#' average observation.  If \code{atmean = FALSE} the function calculates
#' average partial effects.
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
#' @return \item{mfxest}{a coefficient matrix with columns containing the
#' estimates, associated standard errors, test statistics and p-values.}
#' \item{fit}{the fitted \code{\link[betareg]{betareg}} object.} \item{dcvar}{a
#' character vector containing the variable names where the marginal effect
#' refers to the impact of a discrete change on the outcome. For example, a
#' factor variable.} \item{call}{the matched call.}
#' @param ... Arguments passed to betareg::betareg.
#' @seealso \code{\link{betaor}}, \code{\link[betareg]{betareg}}
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
#' betamfx(y~x|x, data=data)
#' @export 
#' 
#' @importFrom stats binomial
#' @importFrom stats coef
#' @importFrom stats complete.cases
#' @importFrom stats dlogis
#' @importFrom stats dnorm
#' @importFrom stats formula
#' @importFrom stats glm
#' @importFrom stats glm.control
#' @importFrom stats hatvalues
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats na.omit
#' @importFrom stats plogis
#' @importFrom stats pnorm
#' @importFrom stats printCoefmat
#' @importFrom stats pt
#' @importFrom stats vcov
#' 

betamfx <-
function(formula, data, atmean = TRUE, robust = FALSE, clustervar1 = NULL, 
                   clustervar2 = NULL, control = betareg.control(), 
                   link.phi = NULL, type = "ML", ...){
  
  res = betamfxest(formula = formula, data = data, atmean, robust, 
                   clustervar1 = clustervar1, clustervar2 = clustervar2,
                   control = control, link.phi = link.phi, type = type, ...)
  
  est = NULL
  est$mfxest = cbind(dFdx = res$mfx$mfx,
                     StdErr = res$mfx$se,
                     z.value = res$mfx$mfx/res$mfx$se,
                     p.value = 2*pt(-abs(res$mfx$mfx/res$mfx$se), df = Inf))
  colnames(est$mfxest) = c("dF/dx","Std. Err.","z","P>|z|")
  rownames(est$mfxest) =  rownames(res$mfx)
  
  est$fit = res$fit
  est$dcvar = rownames(res$mfx[res$mfx$discretechgvar==1,])  
  est$call = match.call() 
  class(est) = "betamfx"
  est
}
