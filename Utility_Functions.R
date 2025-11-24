##################################################
### function to specify survival formula 
### Input:
### TimeStart: variable for the start time. NA for those left censoring
### TimeEnd  : variable for the end time. NA for those right censoring
### xVar     : a vector of covariate names added to model survival


SurvFormulaGen = function(TimeStart = "survStartFmt",
                          TimeEnd = "survEndFmt",
                          xVar=c("edu","sexm")){     #  Xadd=c("edu","sex") or  Xadd="1"
  m = paste("Surv( time=",TimeStart,", time2=", TimeEnd ,
            ", type =  'interval2') ~", paste(xVar,collapse="+"), sep="" )
  
  return(formula(m))
}


##################################################
### function to impute surv time for those left/right censoring
### Input:
### mfit : complete survival fit result from survreg
### dt: dataset that provided to survreg to obtain mfit

# mfit = msurv.fit
# dt = survdatause

SurvImpute = function(mfit, dt){    
  
  if (mfit$dist == 'weibull'){
    dt[!is.na(survStartFmt) & !is.na(survEndFmt), 
       OnsetImp := extrunc('weibull',a=survStartFmt,b=survEndFmt,                           
                           shape = 1/mfit$scale,scale = exp(surlinearPred)),
       by=subjid]
    
    # right censoring
    dt[!is.na(survStartFmt) & is.na(survEndFmt), 
       OnsetImp := extrunc('weibull',a=survStartFmt,b=Inf,                           
                           shape = 1/mfit$scale,scale = exp(surlinearPred)),
       by=subjid]
    
    # left censoring
    dt[is.na(survStartFmt) & !is.na(survEndFmt), 
       OnsetImp := extrunc('weibull',a=0,b=survEndFmt,                           
                           shape = 1/mfit$scale,scale = exp(surlinearPred)),
       by=subjid]
    
  }
  
  if (mfit$dist == 'lognormal'){
    dt[!is.na(survStartFmt) & !is.na(survEndFmt), 
       OnsetImp := extrunc('lnorm',a=survStartFmt,b=survEndFmt,                           
                           shape = 1/mfit$scale,scale = exp(surlinearPred)),
       by=subjid]
    
    # right censoring
    dt[!is.na(survStartFmt) & is.na(survEndFmt), 
       OnsetImp := extrunc('lnorm',a=survStartFmt,b=Inf, 
                           meanlog = surlinearPred, sdlog = mfit$scale),
       by=subjid]
    
    # left censoring
    dt[is.na(survStartFmt) & !is.na(survEndFmt), 
       OnsetImp := extrunc('lnorm',a=0,b=survEndFmt,                           
                           meanlog = surlinearPred, sdlog = mfit$scale),
       by=subjid]
  }      
  
  return(dt$OnsetImp)
  
}


##################################################
### function to approximate the positive part func

approx.func <- function(x,epsilon =0.01){
  y <- pmax(x, 0)
  y[abs(x)<=epsilon] <- (x[abs(x)<=epsilon]+epsilon)^2/(4*epsilon)
  return(y)
}


##################################################
### function to get derivative of approx.func

approx.deriv <- function(x,epsilon=0.01){
  y <- rep(1, length(x))
  y[x< (-epsilon)] <- 0
  y[abs(x)<=epsilon] <- (x[abs(x)<=epsilon]+epsilon)/(2*epsilon)
  return(y)
}


##################################################
### function to specify lme4 formula 
### Input:
### TimeStart: variable for the start time. NA for those left censoring
### TimeEnd  : variable for the end time. NA for those right censoring
### xVar     : a vector of covariate names added to model survival

# yVar = "motscore"
# idVar="subjid"
# tVar = c("visdyAge", "center.t")
# xVar=c("edu","sexf")

lmerFormulaGen = function(yVar, idVar="subjid", tVar, xVar){     #  Xadd=c("edu","sex") or  Xadd="1"
  
  m = paste(yVar ," ~", paste(tVar ,collapse="+"),"+", 
            paste(xVar,collapse="+"),"+",  paste(" (1|",idVar,")",sep=""), sep="" )
 
  return(formula(m))
}


##################################################
### function to get initial values for the cp model
### Note: Fix Delta, then apply 'lmer'.

### Input:
### dt.fl:  dataset of longitudinal measures
### Delta_seq: The sequence of candidates for delta
### tVar : variable indicates t in dt.fl
### onsetVar: variable that indicates imputed onset time
### xVar : covariates for adjustment
### yVar : longitudinal outcome

### Output: a list with
###  coef: obtained initial value
###  deltaseq: acutal applied grid for deltaseq
###  loglikseq: loglikseq corresponds to deltaseq
  
cpinitFUN <- function(dt.fl,Delta_seq=seq(from=0,to=20,by=0.5),
                      tVar = "visdyAge", onsetVar = "OnsetImp",
                      xVar,yVar= "motscore"){
  
  ###bound the Delta_seq to avoid potential identification problem
  
  t <- dt.fl[[match(tVar,colnames(dt.fl))]]  # return as a vector
  u.val <- dt.fl[[match(onsetVar,colnames(dt.fl))]]
  
  center.t <- t - u.val
  Delta_seq.lower <- -max(center.t)
  Delta_seq.upper <- -min(center.t)
  Delta_seq <- Delta_seq[which(Delta_seq>Delta_seq.lower & Delta_seq<Delta_seq.upper)]
  
  ###lmer:
  require(lme4)
  
  rltList   <- vector("list", length(Delta_seq))
  loglikseq <- numeric(length(Delta_seq)) 
  
  tVar <- c(tVar, "center.t.rev")
  
  lmerFormula.use <- lmerFormulaGen(xVar = xVar,yVar=yVar ,tVar =tVar)
  
  for (j in seq(Delta_seq)){
    dt.fl[,center.t.rev := pmax(center.t + Delta_seq[j], 0)]
    fl.fit <- lmer(formula =  lmerFormula.use, data= dt.fl)
    loglikseq[j] <-   as.numeric(logLik(fl.fit))
    rltList[[j]] <-  fl.fit
  }
  
  j.index <- which.max(loglikseq)    # default to return the first that equals the max value
  fl.fit.index <- rltList[[j.index ]]
  delta.index  <- Delta_seq[j.index]
  
  coef.index <- c(fixef(fl.fit.index)[1:3], sigma(fl.fit.index),  sqrt(as.numeric(VarCorr(fl.fit.index))),
                  fixef(fl.fit.index)[-c(1:3)],delta.index )
  names(coef.index) <- c("a","b", "c", "sigma","rho",xVar,"Delta")
  
  returnList <- list(coef = coef.index, loglikseq = loglikseq,deltaseq = Delta_seq)
  
  return(returnList)
  
}


##################################################
### function to get loglik of the cp model for one subject
### Note: We take left, right, converted subjects 
###    for converted subjects, it uses interval censoring format

# dt.fl = datause[subjid=="R007549795"]
# dt.fl = datause[subjid=="101353"]
# modelpara =init
# xVar = Xadd
# yVar = "motscore"
# tVar = "visdyAge"
# onsetVar = "OnsetImp"
# survStartVar ="survStartFmt.alt"
# survEndVar ="survEndFmt.alt"
# mfit = msurv.fit
# epsilon=0.01
# bound.lower=0
# bound.upper=200

cpLikTermFUN.IntlCen <- function(dt.fl, yVar="motscore",tVar = "visdyAge", 
                                 onsetVar = "OnsetImp", xVar,
                                 survStartVar ="survStartFmt", survEndVar ="survEndFmt",
                                 modelpara,mfit,
                                 epsilon=0.01,bound.lower=0,bound.upper=200){
  
  y <- dt.fl[[match(yVar,colnames(dt.fl))]]  # return as a vector
  t <- dt.fl[[match(tVar,colnames(dt.fl))]]  # return as a vector
  p <- length(y)
  
  rho <- modelpara['rho']
  sigma <- modelpara['sigma']
  taup2 <- rho^2/(sigma^2+p*rho^2)
  a <- modelpara['a']
  b <- modelpara['b']
  c <- modelpara['c']
  Delta <- modelpara['Delta']
  
  if (length(xVar)>0){
    for (j in 1:length(xVar)){
      y <- y - dt.fl[[match(xVar[j],colnames(dt.fl))]]*  modelpara[xVar[j]]
    }
  }
  
  start.t <- dt.fl[[match(survStartVar,colnames(dt.fl))]][1] # all elements are the same
  end.t   <- dt.fl[[match(survEndVar,colnames(dt.fl))]][1] # all elements are the same
  surlinearPred.t <- dt.fl$surlinearPred[1]
  
  # if converted 
  
  if (is.na(start.t)==0 & is.na(end.t)==0){
    
    bound.lower <-  start.t
    bound.upper <- end.t
    
  }else{
    
    if (is.na(end.t)>0){ # right censoring
      bound.lower <-  start.t
    }
    if  (is.na(start.t)>0){ # left censoring
      bound.upper <- end.t
    }
  }
  
  GL <- gaussLegendre(30,bound.lower,bound.upper)   
  
  
  functemp <- function(uval){
    r <- y-a-b*t-approx.func(t-uval+sqrt(Delta^2),epsilon)*c
    
    inexp <- (taup2*(sum(r))^2-sum(r^2))/2/(sigma^2)
    
    temp.val <- dsurvreg(uval,mean=surlinearPred.t,
                         scale=mfit$scale,distribution = mfit$dist)*exp(inexp)
    
    return(temp.val)                                         
  }
  
  
  res <- sum((GL$w) * Vectorize(functemp)(GL$x))
  
  val <- (sqrt(2*pi))^(-p)*sqrt(taup2/(rho^2)*(sigma)^(2-2*p))*res        
  
  return(val)
  
}


##################################################
### function to approximate the positive part func
approx.funcA2 <- function(x, epsilon = 0.01) {
  x <- as.numeric(x)                     
  y <- pmax(x, 0)
  idx <- abs(x) <= epsilon               
  if (any(idx)) {
    y[idx] <- ((x[idx] + epsilon)^2) / (4 * epsilon)
  }
  y
}

# --- helper to compute subject SSE for a proposed Δ_i ---
sse_subject <- function(df_i, a_hat, b_hat, c_hat, beta_hat, eps){
  y  <- df_i[[yVar]]
  t  <- df_i[[tVar]]
  ct <- df_i$center_t
  Xb <- if (length(beta_hat)) as.numeric(as.matrix(df_i[, names(beta_hat), with = FALSE]) %*% beta_hat) else 0
  a_i <- a_hat + df_i$RE_intercept[1]  # subject BLUP intercept
  function(Delta){
    mu <- a_i + b_hat * t + c_hat * approx.funcA2(ct + Delta, eps) + Xb
    sum((y - mu)^2)
  }
}

# --- helper: identify converters (for Delta later; NOT used in the fit) ---
get_converters <- function(dd) {
  if ("Censoring" %in% names(dd)) {
    dd[, .(Cens = first(as.character(Censoring))), by = idVar][
      Cens == "interval", get(idVar)]
  } else if (all(c("survStartFmt","survEndFmt") %in% names(dd))) {
    dd[, .(is_interval = any(is.finite(survStartFmt) & is.finite(survEndFmt) &
                               survEndFmt > survStartFmt)), by = idVar][
                                 is_interval == TRUE, get(idVar)]
  } else {
    # Fallback: U within observed time window
    dd[, .(U = first(get(onsetVar)),
           tmin = min(get(tVar), na.rm = TRUE),
           tmax = max(get(tVar), na.rm = TRUE)), by = idVar][
             !is.na(U) & U >= tmin & U <= tmax, get(idVar)]
  }
}


# --- helper: subject SSE as a function of tau ---
sse_subject_tau <- function(df_i, a_hat, b_hat, c_hat, beta_hat, smooth_eps) {
  y <- df_i[[yVar]]
  t <- df_i[[tVar]]
  Xb <- if (length(beta_hat)) as.numeric(as.matrix(df_i[, names(beta_hat), with = FALSE]) %*% beta_hat) else 0
  a_i <- a_hat + df_i$RE_intercept[1]  # subject BLUP intercept
  function(tau) {
    mu <- a_i + b_hat * t + c_hat * approx.funcA2(t - tau, smooth_eps) + Xb
    sum((y - mu)^2)
  }
}


# Vectorized truncated Weibull draw for (a,b] per row
rtrunc_weibull_vec <- function(a, b, shape, scale) {
  n <- length(a); out <- rep(NA_real_, n)
  
  a <- as.numeric(a); b <- as.numeric(b)
  shape <- as.numeric(shape); scale <- as.numeric(scale)
  
  # masks
  int  <- is.finite(a) & is.finite(b) & (b > a)
  rgt  <- is.finite(a) & !is.finite(b)
  lft  <- !is.finite(a) & is.finite(b)
  tie  <- is.finite(a) & is.finite(b) & (b <= a)   # treat as exact at b
  
  if (any(int)) {
    Fa <- pweibull(a[int], shape = shape[int], scale = scale[int])
    Fb <- pweibull(b[int], shape = shape[int], scale = scale[int])
    u  <- runif(sum(int), Fa, Fb)
    out[int] <- qweibull(u, shape = shape[int], scale = scale[int])
  }
  if (any(rgt)) {
    Fa <- pweibull(a[rgt], shape = shape[rgt], scale = scale[rgt])
    u  <- runif(sum(rgt), Fa, 1)
    out[rgt] <- qweibull(u, shape = shape[rgt], scale = scale[rgt])
  }
  if (any(lft)) {
    Fb <- pweibull(b[lft], shape = shape[lft], scale = scale[lft])
    u  <- runif(sum(lft), 0, Fb)
    out[lft] <- qweibull(u, shape = shape[lft], scale = scale[lft])
  }
  if (any(tie)) {
    out[tie] <- b[tie]  # exact or zero-width interval -> set to right bound
  }
  out
}

# Vectorized truncated Lognormal draw for (a,b] per row
rtrunc_lnorm_vec <- function(a, b, meanlog, sdlog) {
  n <- length(a); out <- rep(NA_real_, n)
  
  a <- as.numeric(a); b <- as.numeric(b)
  meanlog <- as.numeric(meanlog); sdlog <- as.numeric(sdlog)
  
  int  <- is.finite(a) & is.finite(b) & (b > a)
  rgt  <- is.finite(a) & !is.finite(b)
  lft  <- !is.finite(a) & is.finite(b)
  tie  <- is.finite(a) & is.finite(b) & (b <= a)
  
  if (any(int)) {
    Fa <- plnorm(a[int], meanlog[int], sdlog[int])
    Fb <- plnorm(b[int], meanlog[int], sdlog[int])
    u  <- runif(sum(int), Fa, Fb)
    out[int] <- qlnorm(u, meanlog[int], sdlog[int])
  }
  if (any(rgt)) {
    Fa <- plnorm(a[rgt], meanlog[rgt], sdlog[rgt])
    u  <- runif(sum(rgt), Fa, 1)
    out[rgt] <- qlnorm(u, meanlog[rgt], sdlog[rgt])
  }
  if (any(lft)) {
    Fb <- plnorm(b[lft], meanlog[lft], sdlog[lft])
    u  <- runif(sum(lft), 0, Fb)
    out[lft] <- qlnorm(u, meanlog[lft], sdlog[lft])
  }
  if (any(tie)) {
    out[tie] <- b[tie]
  }
  out
}