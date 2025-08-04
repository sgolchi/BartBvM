source('code/funs.R')

################################################################################
# PlatinumCAN design:                                                          
# Bayesian operating characteristics for design comparison                     
################################################################################
df = read.csv("data/post_sum_train_350.csv")
N.train = 350
dfs = df %>% group_by(h1, h2, h3, h4, cp, beta) %>% 
  summarize(lambda = mean(lambda), 
            p1 = mean(p1>0.99), p2 = mean(p2>0.98),
            p3 = mean(p3>0.975))

################################################################################
# Generate the test set from the design prior                                  
################################################################################
n.test = 100000
X2 = data.frame(h1 = exp(rnorm(n.test, log(0.055), 0.15)),
                h2 = exp(rnorm(n.test, log(0.095), 0.15)),
                h3 = exp(rnorm(n.test, log(0.040), 0.15)),
                h4 = exp(rnorm(n.test, log(0.020), 0.15)),
                cp = runif(n.test, 0.01, 0.05),
                beta = rnorm(n.test, log(1.3), 0.075))
X1 = data.frame(cbind(dfs$h1, dfs$h2, dfs$h3, dfs$h4, dfs$cp, dfs$beta))
Y = log(dfs$lambda)
names(X1) = c('h1', 'h2', 'h3', 'h4', 'cp', 'beta')

################################################################################
# Obtain predictions of the variance of the sampling distribution via          
# Soft (probabilistic) BART over the test set                                  
################################################################################
bartfit = SoftBart::softbart(X = X1, 
                             Y = Y, X_test = X2)
mpred = bartfit$y_hat_test_mean
low = apply(bartfit$y_hat_test, 2, quantile, p = 0.025)
upp = apply(bartfit$y_hat_test, 2, quantile, p = 0.975)
dfall = data.frame(lambdahat = mpred, lamlow = low, lamupp = upp)
dfall = cbind(X2, dfall)


################################################################################
# Estimate cumulative integrated power, assurance and utility based operating 
# characteristics for two designs:
################################################################################
# Design 1: three analysis with decision thresholds 0.99, 0.98, and 0.975.     
################################################################################
powerSeq = NULL
for (i in 1:nrow(dfall)) {
  powerSeq = rbind(powerSeq, 
                   pwr_seq(dfall$beta[i], 0, exp(dfall$lambdahat[i]), 
                           ns = c(350, 500, 700), us  = c(0.99, 0.98, 0.975)))
}
powerSeq1 = as.data.frame(powerSeq)
assur1 = apply(powerSeq1, 2, mean)
EN1 = round(sum(c(assur1[1], diff(assur1))*c(350, 500, 700)), 0)
names(powerSeq1) = paste(c('p'), 1:ncol(powerSeq), sep = '')
dfdes1 = cbind(dfall, powerSeq1)
bgrid = seq(-0.05, 0.6, 0.01)
dfdes1$beta2 = NULL
for (i in 1:nrow(dfdes1)) {
  dfdes1$beta2[i] = bgrid[which.min(abs(dfdes1$beta[i] - bgrid))]
}
dfsum1 = dfdes1 %>% group_by(beta2) %>% 
  summarize(ip1 = mean(p1), ip2 = mean(p2), ip3 = mean(p3))
dfsum01= reshape(data = dfsum1,
                 idvar = 'beta2',
               varying = 2:4, 
               sep= "",
               timevar= "interim",
               times = c(1, 2, 3),
               new.row.names= 1:189,
               direction = "long")
dfsum01$analysis = as.factor(dfsum01$interim)
ggplot(dfsum01, aes(x = beta2, y = ip, color = analysis)) + 
  geom_line(linewidth = 1) +
  ylab('integrated power') +
  xlab('log(HR)') +
  scale_color_grey(start=0.8, end=0.3) + 
  geom_vline(xintercept = log(1.3), linetype = 3, color = 'grey40') +
  geom_vline(xintercept = log(1), linetype = 3, color = 'grey40') +
  geom_hline(yintercept = 0.425, linetype = 3, color = 'grey40') +
  geom_hline(yintercept = 0.7, linetype = 3, color = 'grey40') +
  geom_hline(yintercept = 0.87, linetype = 3, color = 'grey40') +
  geom_hline(yintercept = 0.04, linetype = 3, color = 'grey40') +
  scale_y_continuous(breaks=c(0.04,0.425,0.7, 0.87)) +
  scale_x_continuous(breaks=round(c(0,log(1.3)),2)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'))
#ggsave(file = 'plots/design1_ns350_500_700_us0.99_0.98_0.975.png')

powerSeq = NULL
dfnull = dfall[dfall$beta<0,]
for (i in 1:nrow(dfnull)) {
  powerSeq = rbind(powerSeq, 
                   pwr_seq(dfnull$beta[i], 0, exp(dfnull$lambdahat[i]), 
                           ns = c(350, 500, 700), us  = c(0.99, 0.98, 0.975)))
  }
powerSeq1_null = as.data.frame(powerSeq)
aT1E1 = apply(powerSeq1_null, 2, mean)
powerSeq = NULL
dfalt = dfall[dfall$beta>0,]
for (i in 1:nrow(dfalt)) {
  powerSeq = rbind(powerSeq, 
                   pwr_seq(dfalt$beta[i], 0, exp(dfalt$lambdahat[i]), 
                           ns = c(350, 500, 700), us  = c(0.99, 0.98, 0.975)))
}
powerSeq1_alt = as.data.frame(powerSeq)
aPow1 = apply(powerSeq1_alt, 2, mean)
iRisk1 = irisk(1000, 10, 1, aT1E1[length(aT1E1)], 1-aPow1[length(aPow1)], EN1)
iRisk1
################################################################################
# Design 2: four analysis with decision thresholds 0.99, 0.99, 0.99, and 0.95     
################################################################################
powerSeq = NULL
for (i in 1:nrow(dfall)) {
  powerSeq = rbind(powerSeq, 
                   pwr_seq(dfall$beta[i], 0, exp(dfall$lambdahat[i]), 
                           ns = c(300, 400, 500, 600), 
                           us  = c(0.99, 0.99, 0.99, 0.95)))
}
powerSeq2 = as.data.frame(powerSeq)
assur2 = apply(powerSeq2, 2, mean)
EN2 = round(sum(c(assur2[1], diff(assur2))*c(300, 400, 500, 600)), 0)
names(powerSeq2) = paste(c('p'), 1:ncol(powerSeq), sep = '')
dfdes2 = cbind(dfall, powerSeq2)
bgrid = seq(-0.05, 0.6, 0.01)
dfdes2$beta2 = NULL
for (i in 1:nrow(dfdes2)) {
  dfdes2$beta2[i] = bgrid[which.min(abs(dfdes2$beta[i] - bgrid))]
  }
dfsum2 = dfdes2 %>% group_by(beta2) %>% 
  summarize(ip1 = mean(p1), ip2 = mean(p2), ip3 = mean(p3), ip4 = mean(p4))
dfsum02= reshape(data = dfsum2,
                 idvar = 'beta2',
                 varying = 2:5, 
                 sep= "",
                 timevar= "interim",
                 times = c(1, 2, 3, 4),
                 new.row.names= 1:252,
                 direction = "long")
dfsum02$analysis = as.factor(dfsum02$interim)
ggplot(dfsum02, aes(x = beta2, y = ip, color = analysis)) + 
  geom_line(linewidth = 1) +
  ylab('integrated power') +
  xlab('log(HR)') +
  scale_color_grey(start=0.8, end=0.3) +
  geom_vline(xintercept = log(1.3), linetype = 3, color = 'grey40') +
  geom_vline(xintercept = log(1), linetype = 3, color = 'grey40') +
  geom_hline(yintercept = 0.36, linetype = 3, color = 'grey40') +
  geom_hline(yintercept = 0.52, linetype = 3, color = 'grey40') +
  geom_hline(yintercept = 0.64, linetype = 3, color = 'grey40') +
  geom_hline(yintercept = 0.88, linetype = 3, color = 'grey40') +
  geom_hline(yintercept = 0.055, linetype = 3, color = 'grey40') +
  scale_y_continuous(breaks=c(0.055,0.36,0.52, 0.64, 0.88)) +
  scale_x_continuous(breaks=round(c(0,log(1.3)),2)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14, face = 'bold'))
#ggsave(file = 'plots/design2_ns300_400_500_600_us0.99_0.99_0.99_0.95.png')

powerSeq = NULL
dfnull = dfall[dfall$beta<0,]
for (i in 1:nrow(dfnull)) {
  powerSeq = rbind(powerSeq, 
                   pwr_seq(dfnull$beta[i], 0, exp(dfnull$lambdahat[i]), 
                           ns = c(300, 400, 500, 600), 
                           us  = c(0.99, 0.99, 0.99, 0.95)))
}
powerSeq2_null = as.data.frame(powerSeq)
aT1E2 = apply(powerSeq2_null, 2, mean)
powerSeq = NULL
dfalt = dfall[dfall$beta>0,]
for (i in 1:nrow(dfalt)) {
  powerSeq = rbind(powerSeq, 
                   pwr_seq(dfalt$beta[i], 0, exp(dfalt$lambdahat[i]), 
                           ns = c(300, 400, 500, 600), 
                           us  = c(0.99, 0.99, 0.99, 0.95)))
}
powerSeq2_alt = as.data.frame(powerSeq)
aPow2 = apply(powerSeq2_alt, 2, mean)
iRisk2 = irisk(1000, 10, 1, aT1E2[length(aT1E2)], 1-aPow2[length(aPow2)], EN2)
iRisk2