## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 7,
                      fig.height = 4)

## ----setup--------------------------------------------------------------------
library(spaci)

## ----simulate-----------------------------------------------------------------
sim <- simulate_spatial_causal(n = 250, delta_u = 2.0, tau_exp = 0.1, seed = 1)
str(sim, max.level = 1)

## ----estimate-----------------------------------------------------------------
res <- spatial_ate(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, seed = 1)
res

## ----forest, fig.alt = "Forest plot of ATT estimates by method"---------------
plot_ate(res, true_att = sim$true_att)

## ----idaps--------------------------------------------------------------------
fit <- idaps(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1, seed = 1)
fit
fit$weights

## ----recoverU-----------------------------------------------------------------
fp <- recoverUplus(sim$Y, sim$Z, sim$X, sim$coords, tau = 0.1)
fp
head(fp$extras$Uhat)

## ----montecarlo, eval = FALSE-------------------------------------------------
# methods <- c("Naive PS", "DAPS", "iDAPS", "recoverU", "recoverU+")
# nsim <- 200; true <- 2
# store <- matrix(NA, nsim, length(methods), dimnames = list(NULL, methods))
# for (s in seq_len(nsim)) {
#   d <- simulate_spatial_causal(n = 250, delta_u = 2.0, tau_exp = 0.1)
#   store[s, ] <- spatial_ate(d$Y, d$Z, d$X, d$coords, seed = s)$ATT
# }
# data.frame(Method = methods,
#            Bias = round(colMeans(store - true, na.rm = TRUE), 3),
#            MSE  = round(colMeans((store - true)^2, na.rm = TRUE), 3))

