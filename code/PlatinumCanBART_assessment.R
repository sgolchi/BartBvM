source('code/funs.R')
################################################################################
# get leave-one-out CV for power estimation bias over the sample for the three 
# analysis in D1                                                               
################################################################################
df = read.csv("data/post_sum_train_350.csv")
N.train = 350
dfs = df %>% group_by(h1, h2, h3, h4, cp, beta) %>% 
  summarize(lambda = mean(lambda), p1 = mean(p1>0.99), p2 = mean(p2>0.98), 
            p3 = mean(p3>0.975))


mpred = low = upp = c()
for (i in 1:nrow(dfs)) {
  X2 = cbind(dfs$h1[i], dfs$h2[i], dfs$h3[i], dfs$h4[i], dfs$cp[i], dfs$beta[i])
  X1 = cbind(dfs$h1[-i], dfs$h2[-i], dfs$h3[-i], dfs$h4[-i], dfs$cp[-i], 
             dfs$beta[-i])
  bartfit = SoftBart::softbart(X = X1, 
                               Y = log(dfs$lambda[-i]),
                               X_test = X2)
  mpred[i] = bartfit$y_hat_test_mean
  low[i] = quantile(bartfit$y_hat_test, p = 0.025)
  upp[i] = quantile(bartfit$y_hat_test, p = 0.975)
}
lamdf = data.frame(lambda = log(dfs$lambda), lambdahat = mpred, low = low, 
                   upp = upp)
ggplot(lamdf, aes(x = lambda, y = lambdahat)) + 
  geom_point(colour = 'grey40', size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'grey50') +
  geom_errorbar(aes(ymin=low, ymax=upp), colour = 'grey40') +
  theme_bw() + xlab('simulated') + ylab('LOOCV BART prediction') +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14))
#ggsave(file = 'plots/lambdaBART_PlatCAN.png')

dfs = cbind(dfs, lamdf)
dfs = dfs %>% mutate(pBvM1 = pow(lambdahat, N.train, beta, 0.99),
                     pBvM1_low = pow(low, N.train, beta, 0.99),
                     pBvM1_upp = pow(upp, N.train, beta, 0.99))
dfs = dfs %>% mutate(pBvM2 = pow(lambdahat, N.train*1.5, beta, 0.98),
                     pBvM2_low = pow(low, N.train*1.5, beta, 0.98),
                     pBvM2_upp = pow(upp, N.train*1.5, beta, 0.98))
dfs = dfs %>% mutate(pBvM3 = pow(lambdahat, N.train*2, beta, 0.975),
                     pBvM3_low = pow(low, N.train*2, beta, 0.975),
                     pBvM3_upp = pow(upp, N.train*2, beta, 0.975))


ggplot(dfs, aes(x = p1, y = pBvM1 - p1)) + 
  geom_point(colour = 'grey40', size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'grey50') +
  #geom_errorbar(aes(ymin=pBvM1_low, ymax=pBvM1_upp), colour = 'grey40') +
  theme_bw() + xlab('simulated power') + ylab('LOOCV "bias"') +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14))
#ggsave(file = 'plots/powerbias_PlatCAN1.png')

ggplot(dfs, aes(x = p2, y = pBvM2-p2)) + 
  geom_point(colour = 'grey40', size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'grey50') +
  #geom_errorbar(aes(ymin=pBvM2_low, ymax=pBvM2_upp), colour = 'grey40') +
  theme_bw() + xlab('simulated power') + ylab('LOOCV "bias"') +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14))
#ggsave(file = 'plots/powerbias_PlatCAN2.png')

ggplot(dfs, aes(x = p3, y = pBvM3-p3)) + 
  geom_point(colour = 'grey40', size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, colour = 'grey50') +
  #geom_errorbar(aes(ymin=pBvM2_low, ymax=pBvM2_upp), colour = 'grey40') +
  theme_bw() + xlab('simulated power') + ylab('LOOCV "bias"') +
  theme(axis.title = element_text(size = 14, face = 'bold'),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 14))
#ggsave(file = 'plots/powerbias_PlatCAN3.png')
