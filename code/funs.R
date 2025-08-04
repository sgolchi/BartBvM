require(rstan)
require(rstanarm)
require(dplyr)
require(ggplot2)
require(lhs)
require(BART)
require(SoftBart)
require(mvtnorm)
require(rjags)
require(coda)
require(foreach)
require(doParallel)
require(doSNOW)

#' Logit transformation
#' 
#' @param x A scalar between 0 and 1
#' @returns A real scalar
#' @export
logit = function(x) log(x/(1-x))

#' Logit inverse transformation
#' 
#' @param x A real scalar
#' @returns A scalar between 0 and 1
#' @export
expit = function(x) exp(x)/(1+exp(x))

#' Sampling distribution of posterior summaries for a logistic regression model
#' 
#' This function generates the Monte Carlo distribution of a number of posterior
#' summaries for a logistic regression model that adjusts for a binary covariate
#' and a binary treatment assignment in a randomized trial
#' 
#' @param b0 A scalar, intercept in logistic regression
#' @param b1 A scalar, coefficient of a covariate
#' @param psi0 A scalar, main effect of treatment
#' @param psi1 scalar, interaction effect treatment*covariate
#' @param N Sample size
#' @param J Number of Monte Carlo iterations
#' @param prior prior for regression coefficients
#' @param prior_intercept prior for the intercept
#' @returns A list including the posterior mean and standard deviation for the 
#' parameter of interest (main effect of the treatment, psi0)
sdist = function(b0, b1, psi0, psi1, N, J = 1000, 
                 prior = normal(location = rep(0, 3), scale = rep(2, 3)), 
                 prior_intercept = normal(location = 0, scale = 2)) {
  mle.main = mle.int = pmean.main = psd.main = psd.int = pmean.int = 
    pos.main = pos.int =  NULL
  for (j in 1:J) {
    x = rbinom(N, 1, 0.5)
    A = rbinom(N, 1, 0.5)
    p = expit(b0 + b1*x + (psi0 + psi1*x)*A)
    y = NULL
    for (n in 1:N) y = c(y, rbinom(1, 1, p[n])) 
    data = data.frame(y = y, x = x, A = A) 
    m0 = glm(y~x*A, family = 'binomial', data = data)
    m1 = stan_glm(y~x*A, family = 'binomial', data = data, refresh = 0, 
                  chains = 1, prior = prior, prior_intercept = prior_intercept)
    pars = as.matrix(m1)
    mle.main[j] = m0$coefficients[3]
    pmean.main[j] = m1$coefficients[3]
    psd.main[j] = m1$ses[3]
    pos.main[j] = mean(pars[,3]>0)
    mle.int[j] = m0$coefficients[4]
    pmean.int[j] = m1$coefficients[4]
    psd.int[j] = m1$ses[4]
    pos.int[j] = mean(pars[,4]>0)
  }
  out = list(mle.main = mle.main, pmean.main = pmean.main, psd.main = psd.main, 
             pos.main = pos.main, mle.int = mle.int, pmean.int = pmean.int, 
             psd.int = psd.int, pos.int = pos.int)
  return(out)
}

#' Design prior weights
#' 
#' This function computes the density of a point in the parameter space under 
#' the design prior
#' 
#' @param par A vector of size 4 the coordinates of a point in 
#' the 4-d parameter space of the simple logistic regression
#' @param means A vector of size 4 means of normal design priors
#' @param sds A vector of size 4 means of normal design priors
#' @returns A scalar density
dpw = function(par = c(c0, c1, c2, c3), means = c(0, 0, 0.3, 0), 
               sds = c(0.6, 0.075, 0.15, 0.05)){
  w = c()
  for (i in 1:length(par)) w[i] = dnorm(par[i], means[i], sds[i])
  return(prod(w))
}

#' BvM estimate of power
#' 
#' @param y A sacalar, the variance parameter (lambda) of the  in log scale
#' @param N Sample size
#' @param psi The hypothesized value of the parameter of interest
#' @param u Decision threshold
#' @returns Power, a scalar between 0 and 1
#' @export
pow = function(y, N, psi, u) pnorm((sqrt(N)/exp(y))*psi - qnorm(u), 0, 1)


#' BvM estimate of power
#' 
#' This function is used to simulate data, implement censoring, approximate
#' posteriors, estimate "lambda" for the first interim analysis, and estimate 
#' posterior probabilities for the sequential design involving survival data
#' 
#' @param n the sample size for the first interim analysis
#' @param c_vec the vector of constants c to reflect how much more data are 
#'       collected in the subsequent stages, c_vec = 1 for non-sequential design
#' @param cutpoints the vector of interval boundaries for the piecewise 
#'       exponential function
#' @param beta a coefficient for the hazard function
#' @param prob: a probability of an observation being randomized to the 
#'       treatment (i.e., x1 = 1)
#' @param censor: a probability indicating what proportion of eligible patients
#'       are censored at the beginning of each interval 
#' @param deltaL: the lower interval endpoint for the hazard ratio hypothesis 
#'       (in terms of beta, not exp(beta))
#' @param R: the number of simulation repetitions used to estimate each sampling 
#'       distribution
#' @param jags_model: the name of the JAGS file 
#' @param prec_beta: the prior precision for the regression coefficient
#' @param prec_logh: the prior precisions for the log hazards
#' @param burnin: the number of burnin iterations for MCMC
#' @param nchains: the number of MCMC chains
#' @param nthin: the thinning parameter for MCMC
#' @param ndraws: the number of MCMC draws that should be retained
#' @param seed: a numeric seed used for random number generation
#' 
#' @returns vector of posterior SD (first column) followed by posterior 
#'       probabilities at each stage of the experiment
getPostSummaries <- function(n, c_vec, cutpoints, hazards, beta, prob, censor,
                             deltaL, R = 10000, jags_model, prec_beta, 
                             prec_logh, burnin = 1000, nchains = 1, nthin = 2, 
                             ndraws = 5000, seed = 1){
  ## this subfunction generates data before accounting for dropout
  ## inputs are a subset of those from the main function
  dataGen <- function(n, cutpoints, hazards, beta, prob) {
    
    ## Generate covariate x1 for treatment assignment
    x <- cbind(rbinom(n, 1, prob))
    
    ## Compute interval widths
    intervals <- c(0, cutpoints)
    widths <- diff(intervals)
    
    ## Identify unique covariate combinations (in this case, x1 = 0 or 1)
    uc <- unique(x)
    
    ## Precompute adjusted hazards for each unique combination
    adjusted_hazards_list <- lapply(1:nrow(uc), function(i) {
      linpred <- sum(beta * c(uc[i, ]))
      hazards * exp(linpred)
    })
    
    ## Initialize results and generate pseudorandom sequence for
    ## CDF inversion
    y <- rep(0, n)
    event <- rep(0, n)
    u <- runif(n)
    
    for (i in seq_along(adjusted_hazards_list)) {
      ## Select patients belonging to treatment group
      group_idx <- which(x[,1] == uc[i,1])
      
      ## Compute adjusted hazards for this group
      adjusted_hazards <- adjusted_hazards_list[[i]]
      
      ## Simulate survival times for this group
      for (j in group_idx) {
        cumulative_hazard <- 0
        for (k in seq_along(widths)) {
          ## Increment cumulative hazard
          incremental_hazard <- adjusted_hazards[k] * widths[k]
          cumulative_hazard <- cumulative_hazard + incremental_hazard
          
          ## Check if event occurs in this interval based on the CDF and u 
          ## realization
          if (u[j] > exp(-cumulative_hazard)) {
            y[j] <- intervals[k] - (log(u[j]) + cumulative_hazard - 
                                      incremental_hazard)/adjusted_hazards[k]
            event[j] <- 1
            break
          }
        }
        
        ## If no event experienced, survival time exceeds last interval 
        ## (administrative censoring)
        if (event[j] == 0) {
          y[j] <- max(intervals)
        }
      }
    }
    
    return(data.frame(y = y, event = event, x1 = x[, 1]))
  }
  
  ## this subfunction censors data to account for dropout
  ## besides dat that comes from the dataGen() function; 
  ## inputs are a subset of those from the main function
  dataDrop <- function(dat, censor, n, cutpoints){
    
    ## generate dropout indicators (constant proportion for each clinical visit) 
    drop.mat <- cbind(matrix(rbinom(n*length(cutpoints), 1, 
                                    rep(censor, each = length(cutpoints)*n)), 
                             byrow = FALSE, nrow = n), rep(1, n))
    
    ## determine first clinical visit that patient does not attend; they will be
    ## right censored at their previous visit
    drop.vec <- apply(drop.mat, 1, function(x){which(x == 1)[1]})
    
    ## update event time and event indicator depending on whether patient 
    ## is censored
    intervals <- c(0, cutpoints)
    for (k in 1:length(cutpoints)){
      dat$event <- ifelse(drop.vec != k, dat$event, 
                          ifelse(dat$y < intervals[k], dat$event, 0))
      dat$y <- ifelse(drop.vec != k, dat$y, 
                      ifelse(dat$y < intervals[k], dat$y, intervals[k]))
    }
    
    ## round survival times up to the nearest day
    dat$y <- ceiling(dat$y)
    
    return(dat)
  }
  
  ## this subfunction processes the data for the JAGS model
  ## dat is the output from dataDrop() function; cutpoints are an input 
  ## from the main function
  processData <- function(dat, cutpoints){
    
    ## obtain summary statistics for the data
    ## Compute interval widths
    intervals <- c(0, cutpoints)
    widths <- diff(intervals)
    
    ## get the total number of events by segment of baseline hazard function
    ev <- dat$event
    tm <- dat$y
    e_seg <- sum((tm<=cutpoints[1])*ev)
    for (k in 2:length(cutpoints)){
      e_seg <- c(e_seg, sum((tm > cutpoints[k-1] & tm<=cutpoints[k])*ev))
    }
    
    ## get the total events for the treatment group
    e_grp <- sum((dat$x1)*ev)
    
    ## get the total time at risk for x1 = 0
    tm0 <- subset(dat, dat$x1 == 0)$y
    e_0 <- sum(pmin(tm0, cutpoints[1]))
    for (k in 2:length(cutpoints)){
      e_0 <- c(e_0, sum(pmin(pmax(tm0 - cutpoints[k-1], 0), widths[k])))
    } 
    
    ## get the total time at risk for x1 = 1
    tm1 <- subset(dat, dat$x1 == 1)$y
    e_1 <- sum(pmin(tm1, cutpoints[1]))
    for (k in 2:length(cutpoints)){
      e_1 <- c(e_1, sum(pmin(pmax(tm1 - cutpoints[k-1], 0), widths[k])))
    } 
    
    ## return the summary statistics for the expedited JAGS model
    return(list(e_seg, e_grp, e_0, e_1))
  }
  
  ## this subfunction returns a posterior sample for beta (conditional and 
  ## marginal estimands coincide); everything except dat_list is an input for
  ## the main function, and dat_list is returned by processData()
  getPosterior <- function(jags_model, dat_list, prec_beta, prec_logh, burnin,
                           nchains, nthin, ndraws){
    
    ## initialize model using summary statistics from processData() function
    model.fit <- jags.model(file=jags_model,
                            data=list(e_seg = dat_list[[1]], 
                                      e_grp = dat_list[[2]], 
                                      e_0 = dat_list[[3]], e_1 = dat_list[[4]], 
                                      K = length(dat_list[[1]]), 
                                      prec_beta = prec_beta, 
                                      prec_logh = prec_logh, zero = 0), 
                            n.chains = nchains, quiet = TRUE)
    
    ## update models and extract posterior draws
    update(model.fit, burnin, progress.bar = "none")
    model.samples <- coda.samples(model.fit, c("beta"), 
                                  n.iter=ndraws, thin=nthin, 
                                  progress.bar = "none")
    
    as.matrix(model.samples)
  }
  
  ## function to get posterior "lambda" term for the beta parameter;
  ## not the "lambda" of hazard ratio = exp(beta); we need the sample size
  ## n corresponding to the analysis at which we estimate lambda; mod.samples
  ## is the output from getPosterior()
  getLambda <- function(mod.samples, n){
    sd(as.numeric(mod.samples))*sqrt(n)
  }
  
  ## this function gets posterior probabilities; interval endpoint deltaL is 
  ## specified in the main function in terms of beta and not hazard 
  ## ratio = exp(beta); mod.samples is the output from getPosterior()
  getPostProbs <- function(mod.samples, deltaL){
    mean(as.numeric(mod.samples) > deltaL)
  }
  
  ## we now implement this in parallel across all analyses
  sim.results <- foreach(i=1:R, .packages=c('rjags', 'coda'), .combine=rbind,
                         .options.snow=opts) %dopar% {
                           
                           ## use previously defined functions to get posterior 
                           ## summaries for the first analysis
                           set.seed(seed + i)
                           ## generate, censor, and process data
                           dat.acc <- dataGen(n = n, cutpoints = cutpoints, 
                                              hazards = hazards, beta = beta, 
                                              prob = prob)
                           dat.acc <- dataDrop(dat.acc, censor = censor, n = n, 
                                               cutpoints = cutpoints)
                           dat.proc <- processData(dat = dat.acc, 
                                                   cutpoints = cutpoints)
                           ## implement posterior approximation
                           jags.path <- paste(getwd(), '/', jags_model, sep='')
                           post.beta <- getPosterior(jags_model = jags.path, 
                                                     dat_list = dat.proc, 
                                                     prec_beta = prec_beta, 
                                                     prec_logh = prec_logh, 
                                                     burnin = burnin, 
                                                     nchains = nchains,
                                                     nthin = nthin, 
                                                     ndraws = ndraws*nthin)
                           
                           ## get "lambda" for first analysis only 
                           ## (can modify this to get more lambda estimates)
                           lambda <- getLambda(mod.samples = post.beta, n = n)
                           
                           ## get posterior probability for first analysis
                           probs <- getPostProbs(mod.samples = post.beta, 
                                                 deltaL = deltaL)
                           
                           ## augment the data set by repeating this process to 
                           ## add data from later stages this can be commented 
                           ## out if not considering the sequential aspect
                           for (j in 2:length(c_vec)){
                             dat.new <- dataGen(n = ceiling((c_vec[j]-c_vec[j-1])*n), 
                                                cutpoints = cutpoints, 
                                                hazards = hazards, beta = beta, 
                                                prob = prob)
                             dat.new <- dataDrop(dat.new, censor = censor, 
                                                 n = ceiling((c_vec[j]-c_vec[j-1])*n), 
                                                 cutpoints = cutpoints)
                             dat.acc <- rbind(dat.acc, dat.new)
                             dat.proc <- processData(dat = dat.acc, 
                                                     cutpoints = cutpoints)
                             post.beta <- getPosterior(jags_model = jags.path, 
                                                       dat_list = dat.proc, 
                                                       prec_beta = prec_beta, 
                                                       prec_logh = prec_logh, 
                                                       burnin = burnin, 
                                                       nchains = nchains,
                                                       nthin = nthin, 
                                                       ndraws = ndraws*nthin)
                             probs.temp <- getPostProbs(mod.samples = post.beta, 
                                                        deltaL = deltaL)
                             probs <- cbind(probs, probs.temp)
                           }
                           c(lambda, probs)
                         }
  
  return(sim.results)
}


#' Sequential power
#' 
#' this function estimates power for a sequential design with a particular set 
#' of design parameters
#' 
#' @param psi_star this is the hypothesized value for the parameter of interest
#' @param psi0 the null value for the parameter of interest
#' @param lambda the variance parameter "lambda" value estimated by BART
#' @param ns: a vector of sample sizes for the sequential design 
#'     (i.e., c(n1, n2, n3) for a three-stage design)
#' @param us a vector of decision thresholds; this vector should have the same 
#'     dimension as ns; the vector should be on the probability scale 
#'     (e.g., c(0.99, 0.95, 0.9))
#' @returns The "cumulative" stopping probabilities for each analysis; 
#'     the final probability is power
#' @export
pwr_seq <- function(psi_star, psi0, lambda, ns, us){
  pwrs <- NULL
  ## use formula 8 from the text to compute power
  pwrs[1] <- pnorm(sqrt(ns[1])*(psi_star-psi0)/lambda - qnorm(us[1]), 0, 1)
  
  ## if we indeed have a sequential design
  if (length(ns) > 2){
    ## construct the correlation between for the "error" epsilon terms
    ## across all analyses; this correlation should be sqrt(n_a/n_b),
    ## where n_a < n_b are the sample sizes for the two relevant analyses 
    cor.mat <- outer(ns, ns, 
                     FUN = function(x, y) sqrt(pmin(x, y) / pmax(x, y)))
    
    ## for each subsequent analysis, calculate the probability that the
    ## experiment (i) did not end in an earlier stage and (ii) stops for
    ## success at analysis i
    
    ## to do this, we need to calculate the probability that a normal
    ## variable with mean sqrt(ns[j])*(psi_star-psi0)/lambda and variance 
    ## 1 is less than qnorm(u[j]) for all analyses j < i and the 
    ## a normal variable with mean sqrt(ns[i])*(psi_star-psi0)/lambda and
    ## variance 1 is greater than qnorm(us[i])
    for (i in 2:length(ns)){
      pwrs[i] <- pmvnorm(lower = c(rep(-Inf, i-1), qnorm(us[i])), 
                         upper = c(qnorm(us[1:(i-1)]), Inf), 
                         mean = sqrt(ns[1:i])*(psi_star-psi0)/lambda, 
                         corr = cor.mat[1:i, 1:i])
    }
  }
  return(cumsum(pwrs))
}

#' Integrated risk
#' 
#' @param c0 A scalar penalty for type I error
#' @param c1 A scalar penalty for type II error
#' @param c2 A scalar cost per sample
#' @param aT1E integrated type I error
#' @param aT2E integrated type II error
#' @param En expected sample size
#' @returns Integrated cost
irisk = function(c0, c1, c2, aT1E, aT2E, En) {
  return(round(c0*aT1E + c1*aT2E + c2*En, 0))
}
