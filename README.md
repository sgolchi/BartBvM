
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BartBvM

<!-- badges: start -->

<!-- badges: end -->

The R scripts provided in the `code` folder implements BART-BvM which is
a method for estimating the operating characteristics of a design in
presence of nuisance parameters. The code reproduces the results of the
simulation study and the example of the BART-BvM manuscript.

- All required functions and package dependencies are included in
  `funs.R`.
- The script `sims_logistic.R` reproduces the results for the logistic
  regression model in Section 3.
- The script `PlatinumCAN_sim.R` generates Monte Carlo based estimates
  of the variance parameter $`\lambda`$ and other posterior summaries
  for the PLATINUM-CAN example in Section 4. The simulated “data” are
  saved as `post_sum_train_350.Rdata` in the `data` folder.
- The script `PlatinumCanBARTBvM.R` implements the BART-BvM method to
  generate various operating characteristics such as integrated power
  and the cost associated with the two group sequential designs in
  Section 4.
- The script `PlatinumCanBART_BvM_assessment.R` evaluates leave-one-out
  cross validated deviations for the $`\lambda`$ parameter and power
  estimates for the PLATINUM-CAN example. These results are presented in
  the supplementary file to the manuscript.
