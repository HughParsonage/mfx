#' Incidence rate ratios for a negative binomial regression.
#' 
#' This function estimates a negative binomial regression model and calculates
#' the corresponding incidence rate ratios.
#' 
#' If both \code{robust=TRUE} and \code{!is.null(clustervar1)} the function
#' overrides the \code{robust} command and computes clustered standard errors.
#' 
#' @aliases negbinirr print.negbinirr
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
#' @param start starting values for the parameters in the
#' \code{\link[MASS]{glm.nb}} model.
#' @param control see \code{\link[stats]{glm.control}}.
#' @return \item{irr}{a coefficient matrix with columns containing the
#' estimates, associated standard errors, test statistics and p-values.}
#' \item{fit}{the fitted \code{\link[MASS]{glm.nb}} object.} \item{call}{the
#' matched call.}
#' @seealso \code{\link{negbinmfx}}, \code{\link[MASS]{glm.nb}}
#' @examples
#' 
#' # simulate some data
#' set.seed(12345)
#' n = 1000
#' x = rnorm(n)
#' y = rnegbin(n, mu = exp(1 + 0.5 * x), theta = 0.5)
#' 
#' data = data.frame(y,x)
#' 
#' negbinirr(formula=y~x,data=data)
#' @export 
negbinirr <-
function(formula, data, robust = FALSE, clustervar1=NULL, 
                     clustervar2=NULL, start = NULL, control = glm.control()){

  res = negbinirrest(formula, data, robust, clustervar1, clustervar2, 
                     start = start, control = control)
  
  est = NULL
  zstat = log(res$irr$irr)*res$irr$irr/res$irr$se
  
  est$irr = cbind(IRR = res$irr$irr,
                        StdErr = res$irr$se,
                        z.value = zstat,
                        p.value = 2*pt(-abs(zstat), df = Inf))
  colnames(est$irr) = c("IRR","Std. Err.","z","P>|z|")
  rownames(est$irr) =  rownames(res$irr)
  
  est$fit = res$fit
  est$call = match.call() 
  class(est) = "negbinirr"
  est
}
