source('code/funs.R')

set.seed(12345)

### Training set and sample size ###
n.train = 40
N.train = 500
a = randomLHS(n.train, 4)

range = cbind(round(seq(-1.2, 1.2, length.out = 100), 2),
              round(seq(-0.15, 0.15, length.out = 100), 2),
              round(seq(0, 0.6, length.out = 100), 2),
              round(seq(-0.1, 0.1, length.out = 100), 2))
samp = a
for (i in 1:n.train) {
  for (j in 1:4) {
    samp[i,j] = range[which.min(abs(a[i,j] - seq(0, 1, length.out = 100))),j]
  }
}

### Simulating lambda over the training set ###
pmean.m  = psd.m = pos.m = mle.m = pmean.i = psd.i = pos.i = mle.i = 
  c0 = c1 = c2 = c3 = ss = NULL
for (i in 1:nrow(samp)) {
  out = sdist(samp[i, 1], samp[i, 2], samp[i, 3], samp[i, 4], N=N.train, 
              J = 1000)
  pmean.m = c(pmean.m, out$pmean.main)
  psd.m = c(psd.m, out$psd.main)
  pos.m = c(pos.m, out$pos.main)
  mle.m = c(mle.m, out$mle.main)
  pmean.i = c(pmean.i, out$pmean.int)
  psd.i = c(psd.i, out$psd.int)
  pos.i = c(pos.i, out$pos.int)
  mle.i = c(mle.i, out$mle.int)
  m = length(out$mle.main)
  c0 = c(c0, rep(samp[i, 1], m))
  c1 = c(c1, rep(samp[i, 2], m))
  c2 = c(c2, rep(samp[i, 3], m))
  c3 = c(c3, rep(samp[i, 4], m))
}
df2 = data.frame(c0 = c0, c1 = c1, c2 = c2, c3 = c3,
                 mle.main = mle.m, pmean.main = pmean.m, psd.main = psd.m, 
                 pos.main = pos.m,
                 mle.int = mle.i, pmean.int = pmean.i, psd.int = psd.i, 
                 pos.int = pos.i)
#save(df2, file = 'data/dfsampn40N500sims1000iter.Rdata')

######################################################
# get leave-one-out CV BART predictions over the sample
load('data/dfsampn40N500sims1000iter.Rdata')
dfs2 = df2 %>% group_by(c0, c1, c2, c3) %>% summarize(vpmean = sd(pmean.main), 
                                                      vmle = sd(mle.main),
                                                      vpost = mean(psd.main),
                                                      pow = mean(pos.main>0.975))
dfs2$lambda = dfs2$vpmean * sqrt(N.train)
dfs2$lambda2 = dfs2$vpost * sqrt(N.train)
mpred = low = upp = c()
for (i in 1:nrow(dfs2)) {
  X2 = cbind(dfs2$c0[i], dfs2$c1[i], dfs2$c2[i], dfs2$c3[i])
  X1 = cbind(dfs2$c0[-i], dfs2$c1[-i], dfs2$c2[-i], dfs2$c3[-i])
  bartfit = SoftBart::softbart(X = X1, 
                               Y = log(dfs2$lambda2[-i]),
                               X_test = X2)
  mpred[i] = bartfit$y_hat_test_mean
  low[i] = quantile(bartfit$y_hat_test, p = 0.025)
  upp[i] = quantile(bartfit$y_hat_test, p = 0.975)
}
lamdf = data.frame(lambda = log(dfs2$lambda2), lambdahat = mpred, low = low, 
                   upp = upp)
ggplot(lamdf, aes(x = lambda, y = lambdahat)) + 
  geom_point(colour = 'grey40', size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'grey50') +
  geom_errorbar(aes(ymin=low, ymax=upp), colour = 'grey40') +
  theme_bw() + xlab('simulated') + ylab('LOOCV BART prediction') +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14))
#ggsave(file = 'plots/lambdaBART.png')


lamdf$pow = pow(lamdf$lambdahat, N.train, psi = dfs2$c2, u = 0.975)
lamdf$plow = pow(lamdf$low, N.train, psi = dfs2$c2, u = 0.975)
lamdf$pupp = pow(lamdf$upp, N.train, psi = dfs2$c2, u = 0.975)
lamdf$powobs = dfs2$pow
ggplot(lamdf, aes(x = powobs, y = pow)) + 
  geom_point(colour = 'grey40', size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'grey50') +
  geom_errorbar(aes(ymin=plow, ymax=pupp), colour = 'grey40') +
  theme_bw() + xlab('simulated') + ylab('LOOCV BART prediction') +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14))
#ggsave(file = 'plots/powerBART.png')


## Now we generate a grid over \Theta_d ##
grid = expand.grid(round(seq(-1.2, 1.2, length.out = 20), 2),
              round(seq(-0.15, 0.15, length.out = 20), 2),
              round(seq(0, 0.6, length.out = 20), 2),
              round(seq(-0.1, 0.1, length.out = 20), 2))
X2 = grid
X1 = data.frame(cbind(dfs2$c0, dfs2$c1, dfs2$c2, dfs2$c3))
Y = log(dfs2$lambda2)
names(X2) = names(X1) = c('c0', 'c1', 'c2', 'c3')
bartfit = SoftBart::softbart(X = X1, 
                             Y = Y,
                             X_test = X2)
mpreds = bartfit$y_hat_test_mean
lows = apply(bartfit$y_hat_test, 2, quantile, p = 0.025)
upps = apply(bartfit$y_hat_test, 2, quantile, p = 0.975)
dfall = data.frame(lambdahat = mpreds, lamlow = lows, lamupp = upps)
dfall = cbind(X2, dfall)

w = apply(dfall[, 1:4], 1, dpw)
dfall$w = w

dflarge = do.call("rbind", replicate(5, dfall, simplify = FALSE))
dflarge$pow = c(pow(dfall$lambdahat, N = 400, psi = dfall$c2, u = 0.975),
              pow(dfall$lambdahat, N = 600, psi = dfall$c2, u = 0.975),
              pow(dfall$lambdahat, N = 800, psi = dfall$c2, u = 0.975),
              pow(dfall$lambdahat, N = 1000, psi = dfall$c2, u = 0.975),
              pow(dfall$lambdahat, N = 1200, psi = dfall$c2, u = 0.975))
dflarge$N = as.factor(c(rep(c(400, 600, 800, 1000, 1200), each = nrow(dfall))))
dfsum = dflarge %>% group_by(c2, N) %>% summarize(int.pow = mean(pow), 
                                                  int.poww = sum(pow*w/sum(w)))
dfsim = data.frame(c2 = dfs2$c2, pow = dfs2$pow, N = as.factor(rep(500, 40)))
ggplot(dfsum, aes(x = c2, y = int.pow, color = N, linetype = N)) + 
  geom_line(linewidth = 1) +
  geom_point(data = dfsim, aes(x = c2, y = pow), size = 2) +
  xlab('Main treatment effect') + ylab('Integrated power') +
  scale_color_brewer(palette="Dark2") +
  theme_bw() + 
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'))
#ggsave(file = 'plots/IntPowerVsN.png')
