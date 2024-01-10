library("survival")

# lambda: hazard rate
# eta: hazard rate of identical exponential censoring distribution
# t_erl: duration of enrollment
# t_flw: duration of follow-up
varfun <- function(lambda, eta, t_erl, t_flw){
  expdiff <- exp(-(t_flw-t_erl)*(lambda+eta)) - exp(-t_flw*(lambda+eta))
  lambda^2*1/((lambda/(lambda+eta))*(1-expdiff/(t_erl*(lambda+eta))))
}

#-------------- sample size calculator -----------#
# theta0: NI margin
# theta1: effect size, log-odds ratio
# alpha: tolerable type I error rate
# beta: tolerable type II error rate
# lambdaC: hazard rate on the control
# trtassign: treatment assignment, i.e., c(n_E/(n_E + n_C), n_C/(n_E + n_C))
# eta: hazard rate of identical exponential censoring distribution
# t_erl: duration of enrollment
# t_flw: duration of follow-up
freqN <- function(alpha, pwr, theta0, theta1,
                  lambdaC, trtassign,
                  eta, t_erl, t_flw, print.out = TRUE){
  effsize <- log(theta0) - log(theta1)
  asymVar <- varfun(lambda = theta1*lambdaC, eta, t_erl, t_flw)/((theta1*lambdaC)^2*trtassign[1]) +
    varfun(lambda = lambdaC, eta, t_erl, t_flw)/(lambdaC^2*trtassign[2])
  
  N <- (qnorm(1-alpha, 0, 1) + qnorm(pwr, 0, 1))^2/effsize^2*asymVar
  
  if(print.out == T){
    cat("Power:", pwr*100, "%\nOne-sided significance level:",
        alpha*100, "%.\nNI margin =",
        theta0, "%\nThe sample size required is:")
    cat("\n", N, "in total, with ", N*trtassign[1], "in the treatment group and ",
        N*trtassign[2], "in the control.")
  }
  
}

freqN(alpha = 0.05, pw = 0.8, theta0 = 1.542, theta1 = 1,
      lambdaC = 0.1274, trtassign = c(0.5, 0.5), 
      eta = 0, t_erl = 4, t_flw = 6)  

freqN(alpha = 0.05, pw = 0.9, theta0 = 1.542, theta1 = 1,
      lambdaC = 0.1274, trtassign = c(0.5, 0.5), 
      eta = 0, t_erl = 4, t_flw = 6)

freqN(alpha = 0.05, pw = 0.8, theta0 = 1.427, theta1 = 1,
      lambdaC = 0.1274, trtassign = c(0.5, 0.5), 
      eta = 0, t_erl = 4, t_flw = 6)  

freqN(alpha = 0.05, pw = 0.8, theta0 = 1.450, theta1 = 1,
      lambdaC = 0.1274, trtassign = c(0.5, 0.5), 
      eta = 0, t_erl = 4, t_flw = 6)  


#################################################################
# Default configuration suggests no use of historical control data
sim.NItrials <- function(lambda.hc = 0.1274, lambda.trt = 0.1965,
                         eta = 0, t_erl = 4, t_flw = 6, theta0, alpha = 0.05,
                         n = c(168, 168)){

  time.cc.sim <- rexp(n[2], lambda.hc) 
  time.cc.censor <- ifelse(rep(eta, n[2])==0, rep(Inf, n[2]),
                           rexp(n[2], eta))  # censoring time in the control
  erl.cc <- runif(n[2])*t_erl  # enrolment time
  status.cc <- rep(1, n[2])  # 
  status.cc[time.cc.sim>time.cc.censor | time.cc.sim>t_flw-t_erl] <- 0
  obs.time.cc <- pmin(time.cc.sim, t_flw-t_erl, time.cc.censor)
  grp.cc <- rep(0, n[2]) # group label for ctr
  
  time.trt.sim <- rexp(n[1], lambda.trt)
  time.trt.censor <- ifelse(rep(eta, n[1])==0, rep(Inf, n[1]),
                            rexp(n[1], eta))
  erl.trt <- runif(n[1])*t_erl  # enrolment time
  status.trt <- rep(1, n[1])  # 
  status.trt[time.trt.sim>time.trt.censor | time.trt.sim>t_flw-t_erl] <- 0
  obs.time.trt <- pmin(time.trt.sim, t_flw-t_erl, time.trt.censor)
  grp.trt <- rep(1, n[1]) # group label for ctr
  
  NIdata <- data.frame(time = c(obs.time.cc, obs.time.trt),
                       status = c(status.cc, status.trt),
                       grp = c(grp.cc, grp.trt)
  )
  
  mod <- coxph(Surv(time, status)~grp, NIdata)
  phi.hat <- summary(mod)$coefficients[,"coef"]
  sigma.phi.hat1 <- summary(mod)$coefficients[,"se(coef)"]
  
  teststat <- (phi.hat - log(theta0))/sigma.phi.hat1
  
  return(
    c(teststat, pnorm(teststat, 0, 1), pnorm(teststat, 0, 1)<=alpha)
  )
}

###############################################
# The following incorporates historical data 
# which had been translated into ESS

#------------------- ESS = 202 -------------------#
# Equivalently, Generalised-Gamma(a=13.7, s=1, f=14.8)
# power

set.seed(123)
simRes <- sapply(1:10000, function(x){
  if(x%%1000 == 0) print(x)
  sim.NItrials(lambda.hc = 0.1819, lambda.trt = 0.1819,
               eta = 0, t_erl = 4, t_flw = 6, theta0 = 1.43, alpha = 0.05,
               n = c(247, 449))
})

mean(simRes[3,]==1)   # 0.8126


# type I error rate
set.seed(123)
simRes <- sapply(1:10000, function(x){
  if(x%%1000 == 0) print(x)
  sim.NItrials(lambda.hc = 0.1274, lambda.trt = 0.1819,
               eta = 0, t_erl = 4, t_flw = 6, theta0 = 1.43, alpha = 0.05,
               n = c(247, 449))
})

mean(simRes[3,]==1)  # 0.0481

#------------------- ESS = 110 -------------------#
# Equivalently, Generalised-Gamma(a=10, s=1, f=11.1)
# power

set.seed(123)
simRes <- sapply(1:10000, function(x){
  if(x%%1000 == 0) print(x)
  sim.NItrials(lambda.hc = 0.1819, lambda.trt = 0.1819,
               eta = 0, t_erl = 4, t_flw = 6, theta0 = 1.43, alpha = 0.05,
               n = c(247, 357))
})

mean(simRes[3,]==1)   # 0.7708


# type I error rate
set.seed(123)
simRes <- sapply(1:10000, function(x){
  if(x%%1000 == 0) print(x)
  sim.NItrials(lambda.hc = 0.1274, lambda.trt = 0.1819,
               eta = 0, t_erl = 4, t_flw = 6, theta0 = 1.43, alpha = 0.05,
               n = c(247, 357))
})

mean(simRes[3,]==1)  # 0.052

