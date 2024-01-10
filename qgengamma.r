# exponential data, generalized Gamma prior for hazard
# f(x) = f/s^a/Gamma(a/f) x x^(a-1) x exp(-(x/s)^f)
# E(X^r) = s^r Gamma( (a+r)/f) / Gamma(a/f)
#
# Note: x ~ gamma(a/f,scale=s^f), y = x^(1/f) ~ genGamma(a,s,f),

qgengamma = function(p,a,s,f) {
  q = qgamma(p,a/f,1/s^f)
  return( q^(1/f))
}

#
# a = 2; s=10; f = 5;
# r = rgengamma(1e6,a,s,f)
# hist(r)
# q = quantile(r,c(0.025,0.5,0.975))
# q
# qgengamma(c(0.025,0.5,0.975),a,s,f)
#
