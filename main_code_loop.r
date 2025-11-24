setwd('/Users/liyuehan/Desktop/Yale/Change-point model simulation project/Code/version 4')


library(data.table)
library(truncdist)
library(survival)
library(icenReg)
library(ggplot2)
library(pracma)
library(dplyr)
library(lme4)
source("Data_Generation.R")
source("Utility_Functions.R")


##########################
### Specify Parameter ####
##########################
n <- 1000 
a <- 0
b <- 0.3
c <- 2
Delta <- 1.5
sigma <- 7
beta <- c(0, 4)
gamma <- c(0.3, -0.3)
left_censored_rate <- 0.6
right_censored_rate <- 0.3
reltol <- 1e-6

### Pre-specified variable names
xVar <- c("X1", "X2")
yVar <- "yij"
tVar <- "visdyAge"
idVar <- "subjid"
dist_assume <- 'weibull' # used to model the survival part


## ---- A1 wrapper (returns named vector a,b,c,sigma,rho,beta*,Delta) ----
fit_A1 <- function(
    datause0,
    msurv.fit,
    xVar, yVar, tVar, idVar = "subjid",
    Delta_seq = seq(0,5,by=0.1),
    reltol = 1e-6,
    survInf = NULL
){
  setDT(datause0)
  if (is.null(survInf)) {
    maxR <- suppressWarnings(max(datause0$survEndFmt, na.rm = TRUE))
    if (!is.finite(maxR)) maxR <- 90
    survInf <- max(maxR + 10, 100)
  }
  
  init_obj <- cpinitFUN(dt.fl = datause0, Delta_seq = Delta_seq,
                        tVar = tVar, xVar = xVar, yVar = yVar)
  init <- init_obj$coef
  
  obj_model <- function(modelpara0){
    datause0[, cpLikTermFUN.IntlCen(
      .SD,
      modelpara = modelpara0,
      mfit = msurv.fit,
      xVar = xVar,
      yVar = yVar,
      bound.upper = survInf
    ), by = c(idVar)][, -sum(log(V1))]
  }
  
  fit <- optim(par = init, fn = obj_model, method = "BFGS",
               control = list(reltol = reltol))
  
  theta <- fit$par
  beta_names <- if (length(xVar)) paste0("beta", seq_along(xVar)) else character(0)
  names(theta) <- c("a","b","c","sigma","rho", beta_names, "Delta")
  theta
}

## ---- A2 wrapper (global Δ; slope on t, hinge on (t-U)+Δ) ----
fit_A2 <- function(
    dt,
    xVar, yVar, tVar, idVar,
    onsetVar = "OnsetImp",
    eps = 0.01,
    Delta_grid = seq(0, 5, by = 0.1)
){
  d <- as.data.table(dt)
  
  # keep converters only
  # interval bounds or Censoring == "interval" or subjects whose onset U lies within their observed time range
  if ("survStartFmt" %in% names(d) && "survEndFmt" %in% names(d)) {
    conv_ids <- d[, .(is_interval = any(is.finite(survStartFmt) & is.finite(survEndFmt) &
                                          survEndFmt > survStartFmt)), by = idVar][
                                            is_interval == TRUE, get(idVar)]
  } else if ("Censoring" %in% names(d)) {
    conv_ids <- d[, .(C = first(as.character(Censoring))), by = idVar][C == "interval", get(idVar)]
  } else {
    conv_ids <- d[, .(U = first(get(onsetVar)),
                      tmin = min(get(tVar)), tmax = max(get(tVar))), by = idVar][
                        !is.na(U) & U >= tmin & U <= tmax, get(idVar)]
  }
  d <- d[get(idVar) %in% conv_ids]
  if (nrow(d) == 0) stop("No converters found for A2.")
  
  # Anchor time at onset for the hinge
  # Compute centered time t*
  if (!onsetVar %in% names(d)) stop(sprintf("'%s' not found in data.", onsetVar))
  d[, center_t := get(tVar) - get(onsetVar)]
  
  # feasible Δ range
  # Upper bound h_i=max(0, -min center_t) so t*+delta >=0
  hi <- max(0, -min(d$center_t, na.rm = TRUE))
  grid <- Delta_grid[Delta_grid >= 0 & Delta_grid <= hi]
  if (!length(grid)) grid <- seq(0, hi, length.out = max(5, ceiling(hi / 0.25)))
  
  # formula: slope on t, hinge on (t*+Δ)
  # y=a+bt+c(t*+delta)+eta+eps
  # t*=t-U_i
  make_fml <- function() {
    rhs_core <- paste(tVar, " + H")
    rhs <- if (length(xVar)) paste(rhs_core, paste(xVar, collapse=" + "), sep=" + ") else rhs_core
    as.formula(paste0(yVar, " ~ ", rhs, " + (1|", idVar, ")"))
  }
  fml <- make_fml()
  
  # profile ML over Δ
  nll <- function(Delta) {
    d[, H := approx.funcA2(center_t + Delta, eps)]
    fit <- try(lmer(fml, data = d, REML = FALSE), silent = TRUE)
    if (inherits(fit, "try-error")) return(Inf)
    -as.numeric(logLik(fit))
  }
  scan <- vapply(grid, nll, numeric(1))
  opt  <- if (hi > 0) optimize(nll, interval = c(0, hi)) else list(minimum = 0)
  Delta_hat <- as.numeric(opt$minimum)
  
  # final REML fit at Δ̂
  d[, H := approx.funcA2(center_t + Delta_hat, eps)]
  fit_final <- lmer(fml, data = d, REML = TRUE)
  
  fe <- fixef(fit_final)
  a_hat  <- unname(fe["(Intercept)"])
  b_hat  <- unname(fe[tVar])
  c_hat  <- unname(fe["H"])
  betas <- if (length(xVar)) {
    bx <- unname(fe[xVar]); names(bx) <- paste0("beta", seq_along(xVar)); bx
  } else numeric(0)
  sigma_hat <- stats::sigma(fit_final)
  rho_hat   <- {
    vc <- VarCorr(fit_final)[[idVar]]
    as.numeric(attr(vc, "stddev")[1])
  }
  
  c(a = a_hat, b = b_hat, c = c_hat,
    sigma = sigma_hat, rho = rho_hat,
    betas,
    Delta = Delta_hat)
}

## ---- helper to build datause0 + msurv.fit for a generated dataset ----
single_run <- function(dt, xVar, dist_assume){
  setorder(dt, subjid, visdyAge)
  
  datause <- as.data.table(dt[, c("subjid","yij","X1","X2","visdyAge","survStartFmt","survEndFmt")])
  setorder(datause, subjid, visdyAge)
  survInf <- pmax(sort(datause[, survEndFmt], decreasing = TRUE)[1] + 10, 100)
  
  datause_surv <- unique(datause[, .(subjid, X1, X2, survStartFmt, survEndFmt)])
  msurv.fit <- survreg(SurvFormulaGen(xVar = xVar), data = datause_surv, dist = dist_assume)
  datause_surv[, surlinearPred := as.numeric(predict(msurv.fit, newdata = datause_surv, type = "lp"))]
  OnsetImp <- SurvImpute(mfit = msurv.fit, dt = datause_surv)
  datause_surv[, OnsetImp := OnsetImp ]
  
  datause0 <- merge(datause, datause_surv[, .(subjid, surlinearPred, OnsetImp)],
                    by = "subjid", all.x = TRUE)
  list(datause0 = datause0, msurv.fit = msurv.fit, survInf = survInf)
}

## ---- run multiple seeds, collect A1 (res1) and A2 (res2) ----
seeds <- 123 + 0:99
res1_list <- vector("list", length(seeds))  # A1 results
res2_list <- vector("list", length(seeds))  # A2 results

for (i in seq_along(seeds)) {
  s <- seeds[i]
  message(sprintf("Seed %d (%d/%d)", s, i, length(seeds)))
  
  # generate one dataset
  dt <- generate_data(n = n, seed = s,
                      left_censored_rate = left_censored_rate,
                      right_censored_rate = right_censored_rate)
  
  # prep survival + merged long data
  prep <- single_run(dt, xVar = xVar, dist_assume = dist_assume)
  datause0  <- prep$datause0
  msurv.fit <- prep$msurv.fit
  survInf   <- prep$survInf
  
  # A1
  p1 <- fit_A1(datause0 = datause0, msurv.fit = msurv.fit,
               xVar = xVar, yVar = yVar, tVar = tVar,
               idVar = idVar, Delta_seq = seq(0,5,0.1),
               reltol = reltol, survInf = survInf)
  res1_list[[i]] <- as.data.table(as.list(p1))
  
  # A2 (needs onsetVar present -> we pass datause0 which has OnsetImp)
  p2 <- fit_A2(dt = datause0, xVar = xVar, yVar = yVar,
               tVar = tVar, idVar = idVar, onsetVar = "OnsetImp",
               eps = 0.01, Delta_grid = seq(0,5,0.1))
  res2_list[[i]] <- as.data.table(as.list(p2))
}

res1 <- rbindlist(res1_list, use.names = TRUE, fill = TRUE)  # A1 estimates by seed
res2 <- rbindlist(res2_list, use.names = TRUE, fill = TRUE)  # A2 estimates by seed

summarize_params <- function(DT, want = c("a","b","c","delta")) {
  setDT(DT)
  # map case-insensitively
  cols <- names(DT)[match(want, tolower(names(DT)))]
  cols <- cols[!is.na(cols)]
  long <- melt(DT[, ..cols], measure.vars = cols,
               variable.name = "param", value.name = "estimate")
  long[, .(mean = mean(estimate, na.rm = TRUE),
           median = median(estimate, na.rm = TRUE),
           sd = sd(estimate, na.rm = TRUE)),
       by = param][order(match(param, cols))]
}

res1_summary_100 <- summarize_params(res1)  # for A1
res2_summary_100 <- summarize_params(res2)  # for A2
res1_summary
res2_summary

