
###################################################
###
### Historical control (exponential data)
### ESS of the generalized Gamma prior (N=0)
### Expected posterior ESS - N (when N\neq 0)
###
###################################################

rm(list = ls())

source("qgengamma.r")
source("rgengamma.r")

set.seed(123)

# posterior simulation (size)
n.sim = 10000

# gen-gamma paramters, a & f, s = 1
afAll = rbind(c(13.7, 14.8),
              c(10, 11.1)
              )

af.labels = apply(afAll, 1, function(e)
  paste("a=", e[1], " f=", e[2]))


# planned N
Nall = c(0, 247)

outAll = array(NA, c(2, nrow(rbind(afAll)), length(Nall), n.sim))


dimnames(outAll) = list(c("MTM.P", "ELIR"),
                        af.labels,
                        Nall,
                        1:n.sim)

for (j.af in 1:nrow(afAll)) {
  for (j.N in 1:length(Nall)) {
    a = afAll[j.af, 1]
    f = afAll[j.af, 2]
    N = Nall[j.N]
    
    for (ss in 1:n.sim) {
      theta = rgengamma(1, a, s = 1, f)
      
      if (N > 0)
        y = rgamma(1, N, theta)
      else
        y = 0
      
      q.up = qgengamma(0.9999, a, 1, f)
      # grid of theta values
      x = seq(0, q.up, length = 10000)
      # log-posterior
      logp = (a - 1 + N) * log(x) - x ^ f - x * y
      mode = x[which.max(logp)]
      
      pp = exp(logp - max(logp))
      pp = pp / sum(pp)
      
      mn = sum(pp * x)
      vr = sum(pp * (x - mn) ^ 2)
      
      theta.post = sample(x,
                          size = 10000,
                          replace = TRUE,
                          p = pp)
      
      outAll[1, j.af, j.N, ss] = (a - 1 + N) + f * (f - 1) * mode ^ f -
        N
      outAll[2, j.af, j.N, ss] = mean(a - 1 + N + f * (f - 1) * theta.post ^ f) - N
      
    }
    cat("a=", a, " f=", f, " N=", N, "\n")
  }
}


for (j in 1:nrow(rbind(afAll))) {
  cat(af.labels[j], "\n")
  print(round(apply(outAll[, j, ,], c(1, 2), mean), 0))
}

