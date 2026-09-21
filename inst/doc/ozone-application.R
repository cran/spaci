## ----setup, include = FALSE---------------------------------------------------
data_path <- system.file("extdata", "analysis_dat.xlsx", package = "spaci")
have_data <- nzchar(data_path) && requireNamespace("readxl", quietly = TRUE)
have_geoR <- requireNamespace("geoR", quietly = TRUE)
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = have_data)

## ----prep---------------------------------------------------------------------
library(spaci)
library(readxl)

dat <- as.data.frame(read_excel(data_path))

outcome     <- "mean4maxOzone"
treatment   <- "SnCR"
NO2         <- "totNOxemissions"        # a mediator, excluded from adjustment
coord_names <- c("Fac.Longitude", "Fac.Latitude")
covariates  <- setdiff(names(dat), c(outcome, treatment, NO2, coord_names))

dat <- dat[complete.cases(
  dat[, c(outcome, treatment, coord_names, covariates)]), ]

Y      <- dat[[outcome]]
Z      <- dat[[treatment]]
coords <- as.matrix(dat[, coord_names])
X      <- as.matrix(dat[, covariates])

c(n = nrow(dat), treated = sum(Z == 1), controls = sum(Z == 0))

## ----units--------------------------------------------------------------------
range(Y)
Y <- Y * 1000          # ppm -> ppb, to match the report's scale

## ----fit----------------------------------------------------------------------
## use geoR (as in the original analysis) when available, else the built-in MLE
engine <- if (have_geoR) "geoR" else "mle"

res <- spatial_ate(Y, Z, X, coords,
                   tau = 0.2, caliper = 0.25,
                   matern_method = engine, seed = 1)
res

## ----forest, fig.alt = "Forest plot of the estimated effect of SCR/SNCR on ozone", fig.width = 7, fig.height = 4----
plot_ate(res, main = "Effect of SCR/SNCR on ozone (ATT)")

## ----seeds, eval = FALSE------------------------------------------------------
# S <- 40
# daps_att  <- vapply(1:S, function(s)
#   daps(Y, Z, X, coords, caliper = 0.25, seed = s)$att, numeric(1))
# idaps_att <- vapply(1:S, function(s)
#   idaps(Y, Z, X, coords, tau = 0.2, caliper = 0.25, seed = s)$att, numeric(1))
# 
# c(DAPS_mean = mean(daps_att),  DAPS_range  = range(daps_att))
# #> DAPS_mean 0.56   range [-0.40, 2.20]   (report: 0.54)
# c(iDAPS_mean = mean(idaps_att), iDAPS_range = range(idaps_att))
# #> iDAPS_mean 0.06  range [-0.58, 0.59]   (report: -0.58)

