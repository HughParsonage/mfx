#' Marginal effects for a logit regression.
#' 
#' This function estimates a binary logistic regression model and calculates
#' the corresponding marginal effects.
#' 
#' If both \code{robust=TRUE} and \code{!is.null(clustervar1)} the function
#' overrides the \code{robust} command and computes clustered standard errors.
#' 
#' @aliases logitmfx print.logitmfx
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
#' @param start starting values for the parameters in the
#' \code{\link[stats]{glm}} model.
#' @param control see \code{\link[stats]{glm.control}}.
#' @return \item{mfxest}{a coefficient matrix with columns containing the
#' estimates, associated standard errors, test statistics and p-values.}
#' \item{fit}{the fitted \code{\link[stats]{glm}} object.} \item{dcvar}{a
#' character vector containing the variable names where the marginal effect
#' refers to the impact of a discrete change on the outcome. For example, a
#' factor variable.} \item{call}{the matched call.}
#' @seealso \code{\link{logitor}}, \code{\link[stats]{glm}}
#' @references William H. Greene (2008). Econometric Analysis (6th ed.).
#' Prentice Hall, N.Y. pp 770-787.
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
#' logitmfx(formula=y~x, data=data)
#' @export
logitmfx <-
function(formula, data, atmean = TRUE, robust = FALSE, clustervar1 = NULL, 
                    clustervar2 = NULL, start = NULL, control = list()){
  res = logitmfxest(formula, data, atmean, robust, clustervar1, clustervar2, start, control)
  
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
  class(est) = "logitmfx"
  est
}
