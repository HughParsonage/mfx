#' Odds ratios for a logit regression.
#' 
#' This function estimates a binary logistic regression model and calculates
#' the corresponding odds ratios.
#' 
#' If both \code{robust=TRUE} and \code{!is.null(clustervar1)} the function
#' overrides the \code{robust} command and computes clustered standard errors.
#' 
#' @aliases logitor print.logitor
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
#' \code{\link[stats]{glm}} model.
#' @param control see \code{\link[stats]{glm.control}}.
#' @return \item{oddsratio}{a coefficient matrix with columns containing the
#' estimates, associated standard errors, test statistics and p-values.}
#' \item{fit}{the fitted \code{\link[stats]{glm}} object.} \item{call}{the
#' matched call.}
#' @seealso \code{\link{logitmfx}}, \code{\link[stats]{glm}}
#' @examples
#' 
#' # simulate some data
#' set.seed(12345)
#' n = 1000
#' x = rnorm(n)
#' 
#' # binary outcome
#' y = ifelse(pnorm(1 + 0.5*x + rnorm(n))>0.5, 1, 0)
#' 
#' data = data.frame(y,x)
#' logitor(formula=y~x, data=data)
#' @export
logitor <-
function(formula, data, robust = FALSE, clustervar1 = NULL, 
                   clustervar2 = NULL, start = NULL, control = list()){
  res = logitorest(formula, data, robust, clustervar1, clustervar2, start, control)
  
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
  class(est) = "logitor"
  est
}
