source('code/funs.R')
###############################################################################
# Monte Carlo simulation with embedded Gibbs sampling for obtaining posterior  
# summaries, including variance estimates over the training set               
###############################################################################
n.train = 60                 # size of the training set
N.train = 350                # training sample size
R = 10                     # number of the Monte Carlo simulations
a = randomLHS(n.train, 6)    # training set: Latin Hypercube Sample 
###############################################################################
# mapping the LHS onto the parameter space                                    
###############################################################################
range = cbind(round(seq(0.047, 0.064, length.out = 100), 3),
              round(seq(0.08, 0.11, length.out = 100), 3),
              round(seq(0.034, 0.046, length.out = 100), 3),
              round(seq(0.017, 0.023, length.out = 100), 3),
              round(seq(0.01, 0.05, length.out = 100), 3),
              round(seq(log(1.2), log(1.4), length.out = 100), 3))
samp = a
for (i in 1:n.train) {
  for (j in 1:6) {
    samp[i,j] = range[which.min(abs(a[i,j] - seq(0, 1, length.out = 100))),j]
  }
}

## set up parallelization
cores = detectCores()
cl = makeSOCKcluster(cores[1]-1)

registerDoSNOW(cl)
pb = txtProgressBar(max = R, style = 3)
progress = function(n) setTxtProgressBar(pb, n)
opts = list(progress = progress)         
###############################################################################
# generating posterior summaries including variance estimates at the training 
# set                                                                         
###############################################################################
df = NULL
for (j in 1:nrow(samp)){
  probs.temp = getPostSummaries(n = N.train, c_vec = c(1, 1.5, 2), 
                                cutpoints = c(7, 14, 21, 28),
                                hazards = samp[j, 1:4], 
                                beta = samp[j, 6], 
                                prob = 0.5, censor = samp[j, 5],
                                deltaL = 0, R = R, jags_model = "code/jags_surv.txt",
                                prec_beta = 0.01, prec_logh = rep(0.01, 4),
                                burnin = 1000, nchains = 2, nthin = 1, 
                                ndraws = 3000, seed = i*100)
  df = rbind(df, data.frame(cbind(do.call("rbind", replicate(R, samp[j,], 
                                                             simplify = FALSE)),
                                  probs.temp)))
}
names(df) = c('h1', 'h2', 'h3', 'h4', 'cp', 'beta', 'lambda', 'p1', 
              'p2', 'p3')
#write.csv(df, paste0("post_sum_train_", N.train, ".csv"), row.names = FALSE)