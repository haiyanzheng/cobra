#' Random numbers from generalized Gamma distribution: f(x) =
#'
#' f(x) = f/s^a/Gamma(a/f) x x^(a-1) x exp(-(x/s)^f)
#' E(X^r) = s^r Gamma( (a+r)/f) / Gamma(a/f)
#'
#' Note: x ~ gamma(a/f,scale=s^f), y = x^(1/f) ~ genGamma(a,s,f),
#'
#' @param n number of simulations
#' @param a shape parameter
#' @param s scale parameter
#' @param f family parameter
#'
#' @return
#' @export
#'
#' @examples
#' a = 2; s = 3; f = 4
#' x = rgengamma(1e6,a,s,f)
#' mean(x)
#' # exact expectation
#' s*gamma( (a+1)/f) / gamma(a/f)
#' mean(x^3)
#' # Exact 3rd moment
#' s^3*gamma( (a+3)/f) / gamma(a/f)
#'
rgengamma = function(n=1,a,s,f) {
x = rgamma(n,a/f,s^(-f)  )
  return(x^(1/f))
}

a = 4; s = 5; f = 3
x = rgengamma(1e6,a,s,f)

r= 3
mean(x^3)
s^r*gamma((a+r)/f)/gamma(a/f)


