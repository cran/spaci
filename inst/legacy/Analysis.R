# =========================================================
# LOAD REQUIRED LIBRARIES
# =========================================================
library(geoR)
library(ggplot2)
library(readxl)
# =========================================================
# CLEAR WORKSPACE
# =========================================================
rm(list = ls())

# =========================================================
# READ DATA
# =========================================================
analysis_dat <- read_excel("analysis_dat.xlsx")

dat <- as.data.frame(analysis_dat)

# =========================================================
# USER SETTINGS
# =========================================================
outcome <- "mean4maxOzone"

treatment <- "SnCR"

NO2 <- "totNOxemissions"

coord_names <- c(
  "Fac.Longitude",
  "Fac.Latitude"
)

# All remaining variables become observed covariates
covariates <- setdiff(
  names(dat),
  c( outcome, treatment, NO2, coord_names
  )
)

# =========================================================
# KEEP COMPLETE CASES
# =========================================================
dat <- dat[
  complete.cases(
    dat[
      ,
      c(
        outcome,
        treatment,
        coord_names,
        covariates
      )
    ]
  ),
]

# =========================================================
# EXTRACT VARIABLES
# =========================================================
Y <- dat[[outcome]]

Z <- dat[[treatment]]

coords <- as.matrix(
  dat[
    ,
    coord_names
  ]
)

Xmat <- as.matrix(
  dat[
    ,
    covariates
  ]
)

n <- nrow(dat)

# =========================================================
# ANALYSIS SETTINGS
# =========================================================
tau_exp <- 0.2

caliper <- 0.25

alpha_grid_step <- 0.1

pi_grid_step <- 0.1

head(dat)
# =========================================================
# FUNCTION: SCALE TO [0,1]
# =========================================================
range01 <- function(x) {
  
  rng <- range(x, na.rm = TRUE)
  
  if (!all(is.finite(rng)) ||
      isTRUE(all.equal(rng[1], rng[2]))) {
    
    if (is.matrix(x)) {
      return(matrix(
        0,
        nrow = nrow(x),
        ncol = ncol(x)
      ))
    }
    
    return(rep(0, length(x)))
  }
  
  (x - rng[1]) / (rng[2] - rng[1])
}

# =========================================================
# FUNCTION: CLIP PROPENSITY SCORES
# =========================================================
clip_ps <- function(p,
                    eps = 1e-6) {
  
  p <- as.numeric(p)
  
  p[!is.finite(p)] <- 0.5
  
  pmin(
    pmax(p, eps),
    1 - eps
  )
}

# =========================================================
# FUNCTION: SAFE SCALE
# =========================================================
safe_scale <- function(x) {
  
  x <- as.numeric(x)
  
  sx <- sd(x, na.rm = TRUE)
  mx <- mean(x, na.rm = TRUE)
  
  if (!is.finite(sx) || sx == 0) {
    return(rep(0, length(x)))
  }
  
  out <- (x - mx) / sx
  
  out[!is.finite(out)] <- 0
  
  out
}

# =========================================================
# FUNCTION: MATCHING ATT
# =========================================================
match_ATT <- function(Dmat,
                      Y,
                      Z,
                      caliper = Inf) {
  
  treated <- which(Z == 1)
  controls <- which(Z == 0)
  
  if (length(treated) == 0 ||
      length(controls) == 0) {
    
    return(list(
      att = NA_real_,
      se = NA_real_,
      n_match = 0,
      n_drop = length(treated),
      pairs = NULL
    ))
  }
  
  available_controls <- controls
  
  diffs <- numeric(0)
  
  pairs <- matrix(
    integer(0),
    ncol = 2
  )
  
  treated_order <- sample(
    treated,
    length(treated),
    replace = FALSE
  )
  
  for (i in treated_order) {
    
    if (length(available_controls) == 0) {
      break
    }
    
    dvec <- Dmat[i, available_controls]
    
    if (all(!is.finite(dvec))) {
      next
    }
    
    jpos <- which.min(dvec)
    dmin <- dvec[jpos]
    
    if (!is.finite(dmin) ||
        dmin > caliper) {
      next
    }
    
    j <- available_controls[jpos]
    
    diffs <- c(
      diffs,
      Y[i] - Y[j]
    )
    
    pairs <- rbind(
      pairs,
      c(i, j)
    )
    
    available_controls <- setdiff(
      available_controls,
      j
    )
  }
  
  n_match <- length(diffs)
  
  n_drop <- length(treated) - n_match
  
  if (n_match < 2) {
    
    return(list(
      att = NA_real_,
      se = NA_real_,
      n_match = n_match,
      n_drop = n_drop,
      pairs = pairs
    ))
  }
  
  list(
    att = mean(diffs),
    se = sd(diffs) / sqrt(n_match),
    n_match = n_match,
    n_drop = n_drop,
    pairs = pairs
  )
}

# =========================================================
# FUNCTION: COVARIATE BALANCE SCORE
# =========================================================
covariate_balance_score <- function(
    pairs,
    Xmat) {
  
  if (is.null(pairs) ||
      nrow(pairs) < 2) {
    return(Inf)
  }
  
  treated_idx <- pairs[, 1]
  control_idx <- pairs[, 2]
  
  X_T <- Xmat[treated_idx, ,
              drop = FALSE]
  
  X_C <- Xmat[control_idx, ,
              drop = FALSE]
  
  sX <- apply(
    Xmat[treated_idx, ,
         drop = FALSE],
    2,
    sd,
    na.rm = TRUE
  )
  
  sX[!is.finite(sX) |
       sX == 0] <- 1
  
  sum(abs(
    (colMeans(X_T) -
       colMeans(X_C)) / sX
  ))
}

# =========================================================
# FUNCTION: IDAPS BALANCE SCORE
# =========================================================
balance_score_idaps <- function(
    pairs,
    Xmat,
    G,
    d_space_raw) {
  
  if (is.null(pairs) ||
      nrow(pairs) < 2) {
    return(Inf)
  }
  
  treated_idx <- pairs[, 1]
  control_idx <- pairs[, 2]
  
  X_T <- Xmat[treated_idx, ,
              drop = FALSE]
  
  X_C <- Xmat[control_idx, ,
              drop = FALSE]
  
  sX <- apply(
    Xmat[treated_idx, ,
         drop = FALSE],
    2,
    sd,
    na.rm = TRUE
  )
  
  sX[!is.finite(sX) |
       sX == 0] <- 1
  
  BX <- sum(abs(
    (colMeans(X_T) -
       colMeans(X_C)) / sX
  ))
  
  dbar <- mean(
    d_space_raw[cbind(
      treated_idx,
      control_idx
    )],
    na.rm = TRUE
  )
  
  if (!is.finite(dbar)) {
    dbar <- Inf
  }
  
  sG <- sd(
    G[treated_idx],
    na.rm = TRUE
  )
  
  if (!is.finite(sG) ||
      sG == 0) {
    sG <- 1
  }
  
  BG <- abs(
    (mean(G[treated_idx]) -
       mean(G[control_idx])) / sG
  )
  
  BX + dbar + BG
}

# =========================================================
# FUNCTION: DAPS TUNING
# =========================================================
select_daps_alpha_covariate_only <- function(
    D_ps,
    d_space,
    Y,
    Z,
    Xmat,
    caliper,
    step = 0.1) {
  
  grid_vals <- seq(
    0,
    1,
    by = step
  )
  
  best_B <- Inf
  best_alpha <- NA_real_
  best_fit <- NULL
  
  for (alpha in grid_vals) {
    
    D_daps_tmp <- alpha * D_ps +
      (1 - alpha) * d_space
    
    diag(D_daps_tmp) <- 0
    
    fit_tmp <- match_ATT(
      D_daps_tmp,
      Y,
      Z,
      caliper
    )
    
    if (is.na(fit_tmp$att) ||
        is.null(fit_tmp$pairs) ||
        nrow(fit_tmp$pairs) < 2) {
      next
    }
    
    B_tmp <- covariate_balance_score(
      fit_tmp$pairs,
      Xmat
    )
    
    if (is.finite(B_tmp) &&
        B_tmp < best_B) {
      
      best_B <- B_tmp
      best_alpha <- alpha
      best_fit <- fit_tmp
    }
  }
  
  if (is.null(best_fit)) {
    
    return(list(
      alpha = NA_real_,
      B = NA_real_,
      fit = list(
        att = NA_real_,
        se = NA_real_,
        n_match = 0,
        n_drop = sum(Z == 1),
        pairs = NULL
      )
    ))
  }
  
  list(
    alpha = best_alpha,
    B = best_B,
    fit = best_fit
  )
}

# =========================================================
# FUNCTION: IDAPS TUNING
# =========================================================
select_idaps_weights <- function(
    D_ps,
    d_space,
    d_exp,
    Y,
    Z,
    Xmat,
    G,
    d_space_raw,
    caliper,
    step = 0.1) {
  
  grid_vals <- seq(
    0,
    1,
    by = step
  )
  
  best_B <- Inf
  best_pi <- c(
    NA_real_,
    NA_real_,
    NA_real_
  )
  
  best_fit <- NULL
  
  for (pi1 in grid_vals) {
    
    for (pi2 in grid_vals) {
      
      pi3 <- 1 - pi1 - pi2
      
      if (pi3 < -1e-12) {
        next
      }
      
      if (pi3 < 0) {
        pi3 <- 0
      }
      
      D_idaps_tmp <- pi1 * D_ps +
        pi2 * d_space +
        pi3 * d_exp
      
      diag(D_idaps_tmp) <- 0
      
      fit_tmp <- match_ATT(
        D_idaps_tmp,
        Y,
        Z,
        caliper
      )
      
      if (is.na(fit_tmp$att) ||
          is.null(fit_tmp$pairs) ||
          nrow(fit_tmp$pairs) < 2) {
        next
      }
      
      B_tmp <- balance_score_idaps(
        fit_tmp$pairs,
        Xmat,
        G,
        d_space_raw
      )
      
      if (is.finite(B_tmp) &&
          B_tmp < best_B) {
        
        best_B <- B_tmp
        
        best_pi <- c(
          pi1,
          pi2,
          pi3
        )
        
        best_fit <- fit_tmp
      }
    }
  }
  
  if (is.null(best_fit)) {
    
    return(list(
      pi = c(
        NA_real_,
        NA_real_,
        NA_real_
      ),
      B = NA_real_,
      fit = list(
        att = NA_real_,
        se = NA_real_,
        n_match = 0,
        n_drop = sum(Z == 1),
        pairs = NULL
      )
    ))
  }
  
  list(
    pi = best_pi,
    B = best_B,
    fit = best_fit
  )
}

# =========================================================
# FUNCTION: MATERN COVARIANCE MATRIX
# =========================================================
matern_cov_matrix <- function(
    dmat,
    sigma2,
    theta,
    nu) {
  
  sigma2 <- max(
    as.numeric(sigma2),
    1e-8
  )
  
  theta <- max(
    as.numeric(theta),
    1e-8
  )
  
  nu <- max(
    as.numeric(nu),
    1e-4
  )
  
  h <- dmat
  
  z <- 2 * sqrt(nu) * h / theta
  
  C <- matrix(
    0,
    nrow = nrow(h),
    ncol = ncol(h)
  )
  
  positive <- z > 0
  
  C[positive] <- sigma2 *
    (z[positive]^nu) *
    besselK(z[positive], nu) /
    (gamma(nu) * 2^(nu - 1))
  
  diag(C) <- sigma2
  
  C[!is.finite(C)] <- 0
  
  C
}

# =========================================================
# FUNCTION: ESTIMATE MATERN PARAMETERS
# =========================================================
estimate_matern_params <- function(
    resid,
    coords) {
  
  resid <- as.numeric(resid)
  
  vres <- var(
    resid,
    na.rm = TRUE
  )
  
  if (!is.finite(vres) ||
      vres <= 0) {
    vres <- 1
  }
  
  fallback <- list(
    sigma2 = 0.7 * vres,
    theta = 0.2,
    nu = 0.5,
    sigma2_eps = 0.3 * vres
  )
  
  fit <- tryCatch({
    
    geodat <- as.geodata(
      data.frame(
        x = coords[, 1],
        y = coords[, 2],
        r = resid
      ),
      coords.col = 1:2,
      data.col = 3
    )
    
    likfit(
      geodat,
      trend = "cte",
      cov.model = "matern",
      ini.cov.pars = c(
        0.7 * vres,
        0.2
      ),
      nugget = 0.3 * vres,
      kappa = 0.5,
      fix.nugget = FALSE,
      fix.kappa = FALSE,
      messages = FALSE
    )
    
  }, error = function(e) NULL,
  warning = function(w) NULL)
  
  if (is.null(fit)) {
    return(fallback)
  }
  
  list(
    sigma2 = fit$cov.pars[1],
    theta = fit$cov.pars[2],
    nu = fit$kappa,
    sigma2_eps = fit$nugget
  )
}

##################
recoverU_core <- function(
    Y,
    Z,
    Xmat,
    G,
    coords,
    d_space_raw,
    include_G_in_PS = FALSE) {
  
  n_obs <- length(Y)
  
  n1 <- sum(Z == 1)
  
  if (n1 < 2 || sum(Z == 0) < 2) {
    
    return(list(
      att = NA_real_,
      se = NA_real_
    ))
    
  }
  
  Xmat <- as.data.frame(Xmat)
  
  dat <- data.frame(
    
    Y = Y,
    
    Z = Z,
    
    Xmat,
    
    G = G
    
  )
  
  ##########################################################
  ## Initial outcome model
  ##########################################################
  
  initial_formula <-
    
    as.formula(
      
      paste(
        
        "Y ~ Z +",
        
        paste(
          
          c(names(Xmat), "G"),
          
          collapse = " + "
          
        )
        
      )
      
    )
  
  fit_initial <- tryCatch(
    
    lm(
      
      initial_formula,
      
      data = dat
      
    ),
    
    error = function(e) NULL
    
  )
  
  if (is.null(fit_initial)) {
    
    return(list(
      
      att = NA_real_,
      
      se = NA_real_
      
    ))
    
  }
  
  resid_initial <- residuals(fit_initial)
  
  ##########################################################
  ## Recover spatial confounder
  ##########################################################
  
  par_hat <- estimate_matern_params(
    
    resid_initial,
    
    coords
    
  )
  
  Sigma_U <- matern_cov_matrix(
    
    d_space_raw,
    
    par_hat$sigma2,
    
    par_hat$theta,
    
    par_hat$nu
    
  )
  
  sigma2_eps_hat <- max(
    
    par_hat$sigma2_eps,
    
    1e-8
    
  )
  
  V <-
    
    Sigma_U +
    
    sigma2_eps_hat *
    
    diag(n_obs)
  
  V <-
    
    V +
    
    1e-6 *
    
    diag(n_obs)
  
  W <- model.matrix(
    
    initial_formula,
    
    data = dat
    
  )
  
  theta_gls <- tryCatch({
    
    Vinv_W <- solve(V, W)
    
    Vinv_Y <- solve(V, Y)
    
    solve(
      
      t(W) %*% Vinv_W,
      
      t(W) %*% Vinv_Y
      
    )
    
  },
  
  error = function(e) NULL)
  
  if (is.null(theta_gls)) {
    
    return(list(
      
      att = NA_real_,
      
      se = NA_real_
      
    ))
    
  }
  
  resid_gls <-
    
    as.numeric(
      
      Y -
        
        W %*% theta_gls
      
    )
  
  Uhat <- tryCatch(
    
    as.vector(
      
      Sigma_U %*%
        
        solve(
          
          V,
          
          resid_gls
          
        )
      
    ),
    
    error = function(e)
      
      rep(
        
        NA_real_,
        
        n_obs
        
      )
    
  )
  
  if (any(!is.finite(Uhat))) {
    
    return(list(
      
      att = NA_real_,
      
      se = NA_real_
      
    ))
    
  }
  
  dat$Uhat <- safe_scale(Uhat)
  
  ##########################################################
  ## RecoverU / RecoverU+
  ##########################################################
  
  rhs <- c(
    
    names(Xmat),
    
    "Uhat"
    
  )
  
  if (include_G_in_PS)
    
    rhs <- c(
      
      rhs,
      
      "G"
      
    )
  
  ps_formula <-
    
    as.formula(
      
      paste(
        
        "Z ~",
        
        paste(
          
          rhs,
          
          collapse = " + "
          
        )
        
      )
      
    )
  
  m0_formula <-
    
    as.formula(
      
      paste(
        
        "Y ~",
        
        paste(
          
          rhs,
          
          collapse = " + "
          
        )
        
      )
      
    )
  
  ps_fit <- tryCatch(
    
    glm(
      
      ps_formula,
      
      family =
        
        binomial(),
      
      data = dat,
      
      control =
        
        glm.control(
          
          maxit = 100
          
        )
      
    ),
    
    error = function(e)
      
      NULL
    
  )
  
  if (is.null(ps_fit)) {
    
    return(list(
      
      att = NA_real_,
      
      se = NA_real_
      
    ))
    
  }
  
  ehat <- clip_ps(
    
    fitted(ps_fit)
    
  )
  
  control_dat <-
    
    dat[
      
      dat$Z == 0,
      
    ]
  
  m0_fit <- tryCatch(
    
    lm(
      
      m0_formula,
      
      data = control_dat
      
    ),
    
    error = function(e)
      
      NULL
    
  )
  
  if (is.null(m0_fit)) {
    
    return(list(
      
      att = NA_real_,
      
      se = NA_real_
      
    ))
    
  }
  
  m0hat <- predict(
    
    m0_fit,
    
    newdata = dat
    
  )
  
  w <-
    
    ehat /
    
    (1 - ehat)
  
  psi <-
    
    (
      
      Z -
        
        (1 - Z) *
        
        w
      
    ) *
    
    (
      
      Y -
        
        m0hat
      
    )
  
  att_hat <-
    
    sum(psi) /
    
    n1
  
  pi_hat <-
    
    n1 /
    
    n_obs
  
  infl <-
    
    psi /
    
    pi_hat -
    
    att_hat
  
  se_hat <-
    
    sd(
      
      infl,
      
      na.rm = TRUE
      
    ) /
    
    sqrt(
      
      n_obs
      
    )
  
  list(
    
    att = att_hat,
    
    se = se_hat
    
  )
  
}

# =========================================================
# RECOVERU
# =========================================================
recoverU_DR_ATT <- function(
    Y,
    Z,
    Xmat,
    G,
    coords,
    d_space_raw) {
  
  recoverU_core(
    Y,
    Z,
    Xmat,
    G,
    coords,
    d_space_raw,
    include_G_in_PS = FALSE
  )
  
}# =========================================================
# RECOVERU+
# =========================================================

recoverUI_DR_ATT <- function(
    Y,
    Z,
    Xmat,
    G,
    coords,
    d_space_raw) {
  
  recoverU_core(
    Y,
    Z,
    Xmat,
    G,
    coords,
    d_space_raw,
    include_G_in_PS = TRUE
  )
  
}

# =========================================================
# COMPUTE SPATIAL DISTANCE
# =========================================================
d_space_raw <- as.matrix(dist(coords))

d_space <- range01(d_space_raw)

diag(d_space) <- 0

# =========================================================
# COMPUTE INTERFERENCE EXPOSURE (G)
# =========================================================
K <- exp(-d_space_raw / tau_exp)

diag(K) <- 0

row_sums <- rowSums(K)

row_sums[row_sums == 0] <- 1

K <- K / row_sums

G <- as.vector(K %*% Z)

# =========================================================
# PROPENSITY SCORE MODEL
# =========================================================
ps_formula <- as.formula(
  
  paste(
    treatment,
    "~",
    paste(covariates,
          collapse = " + ")
  )
  
)

ps_model <- glm(
  
  ps_formula,
  
  family = binomial(link = "logit"),
  
  data = dat
  
)

ps_hat <- clip_ps(
  
  fitted(ps_model)
  
)

D_ps <- abs(
  
  outer(
    ps_hat,
    ps_hat,
    "-"
  )
  
)

D_ps <- range01(D_ps)

diag(D_ps) <- 0

# =========================================================
# EXPOSURE DISTANCE
# =========================================================
d_exp <- as.matrix(
  
  dist(G)
  
)

d_exp <- range01(d_exp)

diag(d_exp) <- 0

# =========================================================
# DAPS
# =========================================================
daps_select <- select_daps_alpha_covariate_only(
  
  D_ps,
  
  d_space,
  
  Y,
  
  Z,
  
  Xmat,
  
  caliper,
  
  alpha_grid_step
  
)

# =========================================================
# iDAPS
# =========================================================
idaps_select <- select_idaps_weights(
  
  D_ps,
  
  d_space,
  
  d_exp,
  
  Y,
  
  Z,
  
  Xmat,
  
  G,
  
  d_space_raw,
  
  caliper,
  
  pi_grid_step
  
)

# =========================================================
# NAIVE PS
# =========================================================
fit_naive <- match_ATT(
  
  D_ps,
  
  Y,
  
  Z,
  
  caliper
  
)

# =========================================================
# DAPS
# =========================================================
fit_daps <- daps_select$fit

# =========================================================
# iDAPS
# =========================================================
fit_idaps <- idaps_select$fit

# =========================================================
# RecoverU
# =========================================================

if(ncol(Xmat) < 2){
  
  stop("RecoverU requires at least two covariates with the current implementation.")
  
}

fit_recoverU <- recoverU_DR_ATT(
  Y,
  Z,
  Xmat,
  G,
  coords,
  d_space_raw
)
# =========================================================
# RecoverU+
# =========================================================
fit_recoverUI <- recoverUI_DR_ATT(
  Y,
  Z,
  Xmat,
  G,
  coords,
  d_space_raw
)

# =========================================================
# RESULTS TABLE
# =========================================================
results <- data.frame(
  
  Method = c(
    
    "Naive PS",
    
    "DAPS",
    
    "iDAPS",
    
    "RecoverU",
    
    "RecoverU+"
    
  ),
  
  ATT = c(
    
    fit_naive$att,
    
    fit_daps$att,
    
    fit_idaps$att,
    
    fit_recoverU$att,
    
    fit_recoverUI$att
    
  ),
  
  SE = c(
    
    fit_naive$se,
    
    fit_daps$se,
    
    fit_idaps$se,
    
    fit_recoverU$se,
    
    fit_recoverUI$se
    
  )
  
)

results$Lower95 <-
  
  results$ATT -
  
  1.96 *
  
  results$SE

results$Upper95 <-
  
  results$ATT +
  
  1.96 *
  
  results$SE


print(results[,-3])

cat("\n")

cat("Optimal DAPS alpha =",
    
    round(daps_select$alpha,4),
    
    "\n")

cat(
  
  "Optimal iDAPS weights =",
  
  paste(
    
    round(idaps_select$pi,4),
    
    collapse = ", "
    
  ),
  
  "\n"
  
)

