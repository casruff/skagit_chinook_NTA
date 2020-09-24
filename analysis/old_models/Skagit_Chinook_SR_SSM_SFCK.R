## ----set_options, echo=FALSE, message=FALSE------------------------------
options(width = 100)
knitr::opts_chunk$set(message = FALSE)
set.seed(2017)
if(file.exists("cnt_time.txt")) {
  file.remove("cnt_time.txt")
}

## ----load_pkgs, message=FALSE, warning=FALSE-----------------------------
if(!require("R2jags")) {
  install.packages("R2jags")
  library("R2jags")
}
if(!require("readr")) {
  install.packages("readr")
  library("readr")
}

if(!require("RCurl")) {
  install.packages("RCurl")
  library("RCurl")
}

if(!require("gsl")) {
  install.packages("gsl")
  library("gsl")
}
if(!require("loo")) {
  install.packages("loo")
  library("loo")
}

if(!require("here")) {
  install.packages("here")
  library("here")
}

## set directory locations
datadir <- here("data")
jagsdir <- here("jags")
analdir <- here("analysis")
savedir <- here("analysis/cache")

## ----define_Re2prec------------------------------------------------------
## better round
Re2prec <- function(x,map="round",prec=1) {
  ## 'map' can be round, floor, or ceiling
  ## 'prec' is nearest value (eg, 0.1 means to nearest tenth); default 1 gives normal behavior
  if(prec<=0) {stop("\"prec\" cannot be less than or equal to 0") }
  do.call(map,list(x/prec))*prec
}


## colVars; from Gelman
## returns the column-wise variance of a matrix
colVars <- function(a) {
	n <- dim(a)[[1]]
	c <- dim(a)[[2]]
	mm <- matrix(.colMeans(a, n, c), n, c, byrow = TRUE)
	return(.colMeans(((a - mm) ^ 2), n, c) * n / (n - 1))
}

## waic; from Gelman
## computes WAIC based on pointwise log-like
waic <- function(log_lik) {
  S <- nrow(log_lik)
  n <- ncol(log_lik)
  lpd <- log(colMeans(exp(log_lik)))
  p_waic <- colVars(log_lik)
  elpd_waic <- lpd - p_waic
  waic <- -2*elpd_waic
  loo_weights_raw <- 1/exp(log_lik-max(log_lik))
  loo_weights_normalized <- loo_weights_raw /
    matrix(colMeans(loo_weights_raw),nrow=S,ncol=n,byrow=TRUE)
  loo_weights_regularized <- pmin(loo_weights_normalized, sqrt(S))
  elpd_loo <- log(colMeans(exp(log_lik)*loo_weights_regularized) /
                    colMeans(loo_weights_regularized))
  p_loo <- lpd - elpd_loo
  pointwise <- cbind(waic,lpd,p_waic,elpd_waic,p_loo,elpd_loo)
  total <- colSums(pointwise)
  se <- sqrt(n*colVars(pointwise))
  return(list(waic=total["waic"],
              elpd_waic=total["elpd_waic"],
              p_waic=total["p_waic"],
              elpd_loo=total["elpd_loo"],
              p_loo=total["p_loo"],
              pointwise=pointwise,
              total=total,
              se=se))
}

## ----get_user_inputs-----------------------------------------------------
## 1. file with escapement data
## [n_yrs x 2] matrix of obs counts; 1st col is calendar yr
fn_esc <- "SkagitSFCKEsc.csv"

## 2. file with age comp data
## [n_yrs x (1+A)]; 1st col is calendar yr
fn_age <- "SkagitSFCKAge.csv"
## min & max ages
age_min <- 2
age_max <- 5
## years, if any, of age-comp to skip; see below
age_skip <- 0

## 3. file with harvest data
## [n_yrs x 2] matrix of obs catch; 1st col is calendar yr
fn_harv <- "SkagitSFCKCatch.csv"

## upper threshold for Gelman & Rubin's (1992) potential scale reduction factor (Rhat).
Rhat_thresh <- 1.1

## URL for example data files
## set to NULL if using a local folder/directory
ex_url <- "https://raw.githubusercontent.com/casruff/Skagit-Chinook-run-reconstruction/master/"

## ----get_escapement_data-------------------------------------------------

dat_yrs <- seq(1983,2016,1)

## escapement
dat_esc <- read.csv(file.path(datadir, fn_esc))
dat_esc <- dat_esc[which(dat_esc$year %in% dat_yrs),]

## number of years of data
n_yrs <- length(dat_yrs)
## get first & last years
yr_frst <- min(dat_yrs)
yr_last <- max(dat_yrs)
## log of escapement
ln_dat_esc <- log(dat_esc$escapement)

## ----get_age_data--------------------------------------------------------
dat_age <- read.csv(file.path(datadir,fn_age))
dat_age <- dat_age[which(dat_age$year %in% dat_yrs),]


## drop year col & first age_min+age_skip rows
dat_age <- dat_age[-(1:(age_min+age_skip)),-1]

## num of age classes
A <- age_max-age_min+1
## add row(s) of NA's for forecast years
#dat_age <- rbind(dat_age,matrix(0,n_fore,A,dimnames=list(n_yrs+seq(n_fore),colnames(dat_age))))
## total num of age obs by cal yr
dat_age[,"sum"] <- apply(dat_age,1,sum)
## row indices for any years with no obs age comp
idx_NA_yrs <- which(dat_age$sum<A,TRUE)
## replace 0's in yrs w/o any obs with NA's
dat_age[idx_NA_yrs,(1:A)] <- NA
## change total in yrs w/o any obs from 0 to A to help dmulti()
dat_age[idx_NA_yrs,"sum"] <- A
## convert class
dat_age <- as.matrix(dat_age)

## ----get_harvest---------------------------------------------------------
## harvest
dat_harv <- file.path(datadir,fn_harv))
dat_harv <- dat_harv[which(dat_harv$year %in% dat_yrs),]


## drop year col 
dat_harv <- (dat_harv$ER)

## ----JAGS_RK_AR----------------------------------------------------------
cat("

model {
  
  ##--------
  ## PRIORS
  ##--------
  ## alpha = exp(a) = intrinsic productivity
  alpha ~ dunif(0.1,10);
  mu_Rkr_a <- log(alpha);
  E_Rkr_a <- mu_Rkr_a + sigma_r/(2 - 2*phi^2);
  
  ## strength of dens depend
  beta ~ dunif(0,0.1);
  
  ## AR(1) coef for proc errors
  phi ~ dunif(-0.999,0.999);
  
  ## process variance for recruits model
  sd_r ~ dunif(0.001,20);
  tau_r <- pow(sd_r,-2);
  sigma_r <- pow(sd_r,2);
  
  ## innovation in first year
  innov_1 ~ dnorm(0,tau_r*(1-phi*phi));
  
  ## obs variance for spawners
  sd_s ~ dunif(0.001,20);
  tau_s <- pow(sd_s,-2);
  sigma_s <- pow(sd_s,2);
  
  ## unprojectable early recruits;
  ## hyper mean across all popns
  Rec_mu ~ dnorm(0,0.001);
  ## hyper SD across all popns
  Rec_sig ~ dunif(0,100);
  ## precision across all popns
  Rec_tau <- pow(Rec_sig,-2);
  ## multipliers for unobservable total runs
  ttl_run_mu ~ dunif(1,5);
  ttl_run_tau ~ dunif(1,20);
  
  ## maturity schedule
  ## unif vec for Dirch prior
  for(i in 1:A) { theta[i] <- 1 }
  ## hyper-mean for maturity
  pi_eta ~ ddirch(theta);
  ## hyper-prec for maturity
  pi_tau ~ dunif(0.001,1e3);
  for(t in 1:(n_yrs-age_min)) { pi_vec[t,1:A] ~ ddirch(pi_eta*pi_tau) }
  
  ##------------
  ## LIKELIHOOD
  ##------------
  ## 1st brood yr requires different innovation
  ## predicted recruits in BY t
  E_ln_Rec[1] <- ln_Sp[1] - beta*Sp[1] + mu_Rkr_a + phi*innov_1;
  tot_ln_Rec[1] ~ dnorm(E_ln_Rec[1],tau_r);
  res_ln_Rec[1] <- tot_ln_Rec[1] - E_ln_Rec[1];
  ## median of total recruits
  tot_Rec[1] <- exp(tot_ln_Rec[1]);
  
  ## R/S
  ln_RS[1] <- tot_ln_Rec[1] - ln_Sp[1];
  
  ## brood-yr recruits by age
  for(a in 1:A) {
    Rec[1,a] <- max(1,tot_Rec[1] * pi_vec[1,a]);
  }
  
  ## brood years 2:(n_yrs-age_min)
  for(t in 2:(n_yrs-age_min)) {
    ## predicted recruits in BY t
    E_ln_Rec[t] <- ln_Sp[t]- beta*Sp[t] + mu_Rkr_a + phi*res_ln_Rec[t-1];
    tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t],tau_r);
    res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
    ## median of total recruits
    tot_Rec[t] <- exp(tot_ln_Rec[t]);
    ## R/S
    ln_RS[t] <- tot_ln_Rec[t] - ln_Sp[t];
    ## brood-yr recruits by age
    for(a in 1:A) {
      Rec[t,a] <- max(1,tot_Rec[t] * pi_vec[t,a]);
    }
  } ## end t loop over year
  
  ## get total cal yr returns for first age_min yrs
  for(i in 1:(age_min+age_skip)) {
    ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
    tot_Run[i] <- exp(ln_tot_Run[i]);
  }
  
  ## get predicted calendar year returns by age
  ## matrix Run has dim [(n_yrs-age_min) x A]
  ## step 1: incomplete early broods
  ## first cal yr of this grp is first brood yr + age_min + age_skip
  for(i in 1:(age_max-age_min-age_skip)) {
    ## projected recruits
    for(a in 1:(i+age_skip)) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## imputed recruits
    for(a in (i+1+age_skip):A) {
      lnRec[i,a] ~ dnorm(Rec_mu,Rec_tau);
      Run[i,a] <- exp(lnRec[i,a]);
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    # predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## step 2: info from complete broods
  ## first cal yr of this grp is first brood yr + age_max
  for(i in (A-age_skip):(n_yrs-age_min-age_skip)) {
    for(a in 1:A) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## get predicted calendar year spawners
  ## first cal yr is first brood yr
  for(t in 1:(n_yrs)) {
    ## obs model for spawners
    Sp[t] <- max(1,tot_Run[t] * (1 - dat_harv[t]));
    ln_Sp[t] <- log(Sp[t]);
    ln_dat_esc[t] ~ dnorm(ln_Sp[t], tau_s);
    lp_esc[t] <- logdensity.norm(ln_dat_esc[t],ln_Sp[t], tau_s);
  }
  
} ## end model description

", file="IPM_RK_AR.txt")


## ----JAGS_RK-------------------------------------------------------------
cat("

model {
  
  ##--------
  ## PRIORS
  ##--------
  ## alpha = exp(a) = intrinsic productivity
  alpha ~ dunif(0.1,10);
  mu_Rkr_a <- log(alpha);
  E_Rkr_a <- mu_Rkr_a + sigma_r/2;
  
  ## strength of dens depend
  beta ~ dunif(0,0.1);
  

  ## process variance for recruits model
  sd_r ~ dunif(0.001,20);
  tau_r <- pow(sd_r,-2);
  sigma_r <- pow(sd_r,2);
  
  ## obs variance for spawners
  sd_s ~ dunif(0.001,20);
  tau_s <- pow(sd_s,-2);
  sigma_s <- pow(sd_s,2);
  
  ## unprojectable early recruits;
  ## hyper mean across all popns
  Rec_mu ~ dnorm(0,0.001);
  ## hyper SD across all popns
  Rec_sig ~ dunif(0,100);
  ## precision across all popns
  Rec_tau <- pow(Rec_sig,-2);
  ## multipliers for unobservable total runs
  ttl_run_mu ~ dunif(1,5);
  ttl_run_tau ~ dunif(1,20);
  
  ## maturity schedule
  ## unif vec for Dirch prior
  for(i in 1:A) { theta[i] <- 1 }
  ## hyper-mean for maturity
  pi_eta ~ ddirch(theta);
  ## hyper-prec for maturity
  pi_tau ~ dunif(0.001,1e3);
  for(t in 1:(n_yrs-age_min)) { pi_vec[t,1:A] ~ ddirch(pi_eta*pi_tau) }
  
  ##------------
  ## LIKELIHOOD
  ##------------
  ## 1st brood yr requires different innovation
  ## predicted recruits in BY t
  # covar[1] <- inprod(gamma,scl_cvrs[1,]);
  E_ln_Rec[1] <- ln_Sp[1] - beta*Sp[1] + mu_Rkr_a;
  tot_ln_Rec[1] ~ dnorm(E_ln_Rec[1],tau_r);
  res_ln_Rec[1] <- tot_ln_Rec[1] - E_ln_Rec[1];
  ## median of total recruits
  tot_Rec[1] <- exp(tot_ln_Rec[1]);
  
  ## R/S
  ln_RS[1] <- tot_ln_Rec[1] - ln_Sp[1];
  
  ## brood-yr recruits by age
  for(a in 1:A) {
    Rec[1,a] <- max(1,tot_Rec[1] * pi_vec[1,a]);
  }
  
  ## brood years 2:(n_yrs-age_min)
  for(t in 2:(n_yrs-age_min)) {
    ## predicted recruits in BY t
    E_ln_Rec[t] <- ln_Sp[t] - beta*Sp[t] + mu_Rkr_a;
    tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t],tau_r);
    res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
    ## median of total recruits
    tot_Rec[t] <- exp(tot_ln_Rec[t]);
    ## R/S
    ln_RS[t] <- tot_ln_Rec[t] - ln_Sp[t];
    ## brood-yr recruits by age
    for(a in 1:A) {
      Rec[t,a] <- max(1,tot_Rec[t] * pi_vec[t,a]);
    }
  } ## end t loop over year
  
  ## get total cal yr returns for first age_min yrs
  for(i in 1:(age_min+age_skip)) {
    ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
    tot_Run[i] <- exp(ln_tot_Run[i]);
  }
  
  ## get predicted calendar year returns by age
  ## matrix Run has dim [(n_yrs-age_min) x A]
  ## step 1: incomplete early broods
  ## first cal yr of this grp is first brood yr + age_min + age_skip
  for(i in 1:(age_max-age_min-age_skip)) {
    ## projected recruits
    for(a in 1:(i+age_skip)) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## imputed recruits
    for(a in (i+1+age_skip):A) {
      lnRec[i,a] ~ dnorm(Rec_mu,Rec_tau);
      Run[i,a] <- exp(lnRec[i,a]);
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    # predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## step 2: info from complete broods
  ## first cal yr of this grp is first brood yr + age_max
  for(i in (A-age_skip):(n_yrs-age_min-age_skip)) {
    for(a in 1:A) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## get predicted calendar year spawners
  ## first cal yr is first brood yr
  for(t in 1:(n_yrs)) {
    ## obs model for spawners
    Sp[t] <- max(1,tot_Run[t] * (1 - dat_harv[t]));
    ln_Sp[t] <- log(Sp[t]);
    ln_dat_esc[t] ~ dnorm(ln_Sp[t], tau_s);
    lp_esc[t] <- logdensity.norm(ln_dat_esc[t],ln_Sp[t], tau_s);
  }
  
} ## end model description

", file="IPM_RK.txt")

## ----JAGS_BH_AR----------------------------------------------------------
cat("

model {
  
  ##--------
  ## PRIORS
  ##--------
  ## alpha = exp(a) = intrinsic productivity
  alpha ~ dunif(0.1,10);
  mu_BH_a <- log(alpha);
  E_BH_a <- mu_BH_a + sigma_r/(2 - 2*phi^2);
  
  ## strength of dens depend
  beta ~ dunif(0,0.01);
  
  ## AR(1) coef for proc errors
  phi ~ dunif(-0.999,0.999);
  
  ## process variance for recruits model
  sd_r ~ dunif(0.001,20);
  tau_r <- pow(sd_r,-2);
  sigma_r <- pow(sd_r,2);
  
  ## innovation in first year
  innov_1 ~ dnorm(0,tau_r*(1-phi*phi));
  
  ## obs variance for spawners
  sd_s ~ dunif(0.001,20);
  tau_s <- pow(sd_s,-2);
  sigma_s <- pow(sd_s,2);
  
  ## unprojectable early recruits;
  ## hyper mean across all popns
  Rec_mu ~ dnorm(0,0.001);
  ## hyper SD across all popns
  Rec_sig ~ dunif(0,100);
  ## precision across all popns
  Rec_tau <- pow(Rec_sig,-2);
  ## multipliers for unobservable total runs
  ttl_run_mu ~ dunif(1,5);
  ttl_run_tau ~ dunif(1,20);
  
  ## maturity schedule
  ## unif vec for Dirch prior
  for(i in 1:A) { theta[i] <- 1 }
  ## hyper-mean for maturity
  pi_eta ~ ddirch(theta);
  ## hyper-prec for maturity
  pi_tau ~ dunif(0.001,1e3);
  for(t in 1:(n_yrs-age_min)) { pi_vec[t,1:A] ~ ddirch(pi_eta*pi_tau) }
  
  ##------------
  ## LIKELIHOOD
  ##------------
  ## 1st brood yr requires different innovation
  ## predicted recruits in BY t
  E_ln_Rec[1] <- mu_BH_a + ln_Sp[1] - log(1 + beta*Sp[1]) + phi*innov_1;
  tot_ln_Rec[1] ~ dnorm(E_ln_Rec[1],tau_r);
  res_ln_Rec[1] <- tot_ln_Rec[1] - E_ln_Rec[1];
  ## median of total recruits
  tot_Rec[1] <- exp(tot_ln_Rec[1]);
  
  ## R/S
  ln_RS[1] <- tot_ln_Rec[1] - ln_Sp[1];
  
  ## brood-yr recruits by age
  for(a in 1:A) {
    Rec[1,a] <- max(1,tot_Rec[1] * pi_vec[1,a]);
  }
  
  ## brood years 2:(n_yrs-age_min)
  for(t in 2:(n_yrs-age_min)) {
    ## predicted recruits in BY t
    E_ln_Rec[t] <- mu_BH_a + ln_Sp[t] - log(1 + beta*Sp[t]) + phi*res_ln_Rec[t-1];
    tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t],tau_r);
    res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
    ## median of total recruits
    tot_Rec[t] <- exp(tot_ln_Rec[t]);
    ## R/S
    ln_RS[t] <- tot_ln_Rec[t] - ln_Sp[t];
    ## brood-yr recruits by age
    for(a in 1:A) {
      Rec[t,a] <- max(1,tot_Rec[t] * pi_vec[t,a]);
    }
  } ## end t loop over year
  
  ## get total cal yr returns for first age_min yrs
  for(i in 1:(age_min+age_skip)) {
    ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
    tot_Run[i] <- exp(ln_tot_Run[i]);
  }
  
  ## get predicted calendar year returns by age
  ## matrix Run has dim [(n_yrs-age_min) x A]
  ## step 1: incomplete early broods
  ## first cal yr of this grp is first brood yr + age_min + age_skip
  for(i in 1:(age_max-age_min-age_skip)) {
    ## projected recruits
    for(a in 1:(i+age_skip)) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## imputed recruits
    for(a in (i+1+age_skip):A) {
      lnRec[i,a] ~ dnorm(Rec_mu,Rec_tau);
      Run[i,a] <- exp(lnRec[i,a]);
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    # predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## step 2: info from complete broods
  ## first cal yr of this grp is first brood yr + age_max
  for(i in (A-age_skip):(n_yrs-age_min-age_skip)) {
    for(a in 1:A) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## get predicted calendar year spawners
  ## first cal yr is first brood yr
  for(t in 1:(n_yrs)) {
    ## obs model for spawners
    Sp[t] <- max(1,tot_Run[t] *(1 - dat_harv[t]));
    ln_Sp[t] <- log(Sp[t]);
    ln_dat_esc[t] ~ dnorm(ln_Sp[t], tau_s);
    lp_esc[t] <- logdensity.norm(ln_dat_esc[t],ln_Sp[t], tau_s);
  }
  
} ## end model description

", file="IPM_BH_AR.txt")

## ----JAGS_BH-------------------------------------------------------------
cat("

model {
  
  ##--------
  ## PRIORS
  ##--------
  ## alpha = exp(a) = intrinsic productivity
  alpha ~ dunif(0.1,10);
  mu_BH_a <- log(alpha);
  E_BH_a <- mu_BH_a + sigma_r/2;
  
  ## strength of dens depend
  beta ~ dunif(0,0.01);
  

  ## process variance for recruits model
  sd_r ~ dunif(0.001,20);
  tau_r <- pow(sd_r,-2);
  sigma_r <- pow(sd_r,2);
  
  ## obs variance for spawners
  sd_s ~ dunif(0.001,20);
  tau_s <- pow(sd_s,-2);
  sigma_s <- pow(sd_s,2);
  
  ## unprojectable early recruits;
  ## hyper mean across all popns
  Rec_mu ~ dnorm(0,0.001);
  ## hyper SD across all popns
  Rec_sig ~ dunif(0,100);
  ## precision across all popns
  Rec_tau <- pow(Rec_sig,-2);
  ## multipliers for unobservable total runs
  ttl_run_mu ~ dunif(1,5);
  ttl_run_tau ~ dunif(1,20);
  
  ## maturity schedule
  ## unif vec for Dirch prior
  for(i in 1:A) { theta[i] <- 1 }
  ## hyper-mean for maturity
  pi_eta ~ ddirch(theta);
  ## hyper-prec for maturity
  pi_tau ~ dunif(0.001,1e3);
  for(t in 1:(n_yrs-age_min)) { pi_vec[t,1:A] ~ ddirch(pi_eta*pi_tau) }
  
  ##------------
  ## LIKELIHOOD
  ##------------
  ## predicted recruits in BY t
  E_ln_Rec[1] <- mu_BH_a + ln_Sp[1] - log(1 + beta*Sp[1]);
  tot_ln_Rec[1] ~ dnorm(E_ln_Rec[1],tau_r);
  res_ln_Rec[1] <- tot_ln_Rec[1] - E_ln_Rec[1];
  ## median of total recruits
  tot_Rec[1] <- exp(tot_ln_Rec[1]);
  
  ## R/S
  ln_RS[1] <- tot_ln_Rec[1] - ln_Sp[1];
  
  ## brood-yr recruits by age
  for(a in 1:A) {
    Rec[1,a] <- max(1,tot_Rec[1] * pi_vec[1,a]);
  }
  
  ## brood years 2:(n_yrs-age_min)
  for(t in 2:(n_yrs-age_min)) {
    ## predicted recruits in BY t
    E_ln_Rec[t] <- mu_BH_a + ln_Sp[t] - log(1 + beta*Sp[t]);
    tot_ln_Rec[t] ~ dnorm(E_ln_Rec[t],tau_r);
    res_ln_Rec[t] <- tot_ln_Rec[t] - E_ln_Rec[t];
    ## median of total recruits
    tot_Rec[t] <- exp(tot_ln_Rec[t]);
    ## R/S
    ln_RS[t] <- tot_ln_Rec[t] - ln_Sp[t];
    ## brood-yr recruits by age
    for(a in 1:A) {
      Rec[t,a] <- max(1,tot_Rec[t] * pi_vec[t,a]);
    }
  } ## end t loop over year
  
  ## get total cal yr returns for first age_min yrs
  for(i in 1:(age_min+age_skip)) {
    ln_tot_Run[i] ~ dnorm(ttl_run_mu*Rec_mu,Rec_tau/ttl_run_tau);
    tot_Run[i] <- exp(ln_tot_Run[i]);
  }
  
  ## get predicted calendar year returns by age
  ## matrix Run has dim [(n_yrs-age_min) x A]
  ## step 1: incomplete early broods
  ## first cal yr of this grp is first brood yr + age_min + age_skip
  for(i in 1:(age_max-age_min-age_skip)) {
    ## projected recruits
    for(a in 1:(i+age_skip)) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## imputed recruits
    for(a in (i+1+age_skip):A) {
      lnRec[i,a] ~ dnorm(Rec_mu,Rec_tau);
      Run[i,a] <- exp(lnRec[i,a]);
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    # predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## step 2: info from complete broods
  ## first cal yr of this grp is first brood yr + age_max
  for(i in (A-age_skip):(n_yrs-age_min-age_skip)) {
    for(a in 1:A) {
      Run[i,a] <- Rec[(age_skip+i)-a+1,a];
    }
    ## total run size
    tot_Run[i+age_min+age_skip] <- sum(Run[i,1:A]);
    ## predicted age-prop vec for multinom
    for(a in 1:A) {
      age_v[i,a] <- Run[i,a] / tot_Run[i+age_min];
    }
    ## multinomial for age comp
    dat_age[i,1:A] ~ dmulti(age_v[i,1:A],dat_age[i,A+1]);
    lp_age[i] <- logdensity.multi(dat_age[i,1:A],age_v[i,1:A],dat_age[i,A+1]);
  }
  
  ## get predicted calendar year spawners
  ## first cal yr is first brood yr
  for(t in 1:(n_yrs)) {
    ## obs model for spawners
    Sp[t] <- max(1,tot_Run[t]*(1-dat_harv[t]));
    ln_Sp[t] <- log(Sp[t]);
    ln_dat_esc[t] ~ dnorm(ln_Sp[t], tau_s);
    lp_esc[t] <- logdensity.norm(ln_dat_esc[t],ln_Sp[t], tau_s);
  }
  
} ## end model description

", file="IPM_BH.txt")

## ----start_timer, include=FALSE------------------------------------------
## start timer
timer_start <- proc.time() 

## ----JAGS_IO, message=FALSE, warning=FALSE, cache=TRUE-------------------
## 1. data to pass to JAGS
dat_jags <- c("dat_age","ln_dat_esc","dat_harv",
              "n_yrs","A","age_min","age_max","age_skip") 


## 2. model params/states for JAGS to return
par_jags <- c("lp_age","lp_esc")


## 3. MCMC control params
## MCMC parameters
mcmc_chains <- 4
mcmc_length <- 2e5
mcmc_burn <- 1e5
mcmc_thin <- 100
## total number of MCMC samples
mcmc_samp <- (mcmc_length-mcmc_burn)*mcmc_chains/mcmc_thin

## Eqn 1
init_vals_AR <- function() {
	list(alpha=2,
	     beta=1/exp(mean(ln_dat_esc, na.rm=TRUE)),
	     pi_tau=1, pi_eta=rep(1,A),
	     pi_vec=matrix(c(0.020,0.219,0.581,0.179),n_yrs-age_min,A,byrow=TRUE),
	     Rec_mu=log(1000),
	     Rec_sig=0.1,
	     tot_ln_Rec=rep(log(1000),n_yrs-age_min),
	     innov_1=0,
	     phi=0.5)
}
## Eqn 2
init_vals <- function() {
	list(alpha=2,
	     beta=1/exp(mean(ln_dat_esc, na.rm=TRUE)),
	     pi_tau=1,
	     pi_eta=rep(1,A),
	     pi_vec=matrix(c(0.020,0.219,0.581,0.179),n_yrs-age_min,A,byrow=TRUE),
	     Rec_mu=log(1000),
	     Rec_sig=0.1,
	     tot_ln_Rec=rep(log(1000),n_yrs-age_min))
}


## empty list for fits
mod_fits <- vector("list", 4)

## list of model info for JAGS
mod_jags <- list(data=dat_jags,
                 parameters.to.save=par_jags,
                 n.chains=as.integer(mcmc_chains),
                 n.iter=as.integer(mcmc_length),
                 n.burnin=as.integer(mcmc_burn),
                 n.thin=as.integer(mcmc_thin))

## fit AR(1) models
mod_jags$inits <- init_vals_AR

## Ricker
mod_jags$model.file <- "IPM_RK_AR.txt"
mod_fits[[1]] <- do.call(jags.parallel, mod_jags)
## B-H
mod_jags$model.file <- "IPM_BH_AR.txt"
mod_fits[[2]] <- do.call(jags.parallel, mod_jags)

## fit no AR(1) models
mod_jags$inits <- init_vals

## Ricker
mod_jags$model.file <- "IPM_RK.txt"
mod_fits[[3]] <- do.call(jags.parallel, mod_jags)

## B-H
mod_jags$model.file <- "IPM_BH.txt"
mod_fits[[4]] <- do.call(jags.parallel, mod_jags)


## ----stop_timer, include=FALSE-------------------------------------------
## stop timer
run_time_in_min <- round(((proc.time()-timer_start)/60)["elapsed"], 1)
cat(run_time_in_min, file="run_time_in_min_all.txt")

## ----get_WAIC------------------------------------------------------------
n_mods <- length(mod_fits)
WAIC <- vector("numeric",n_mods)
LOO <- vector("numeric",n_mods)

## extract log densities from JAGS objects
for(i in 1:n_mods) {
  ldens <- cbind(mod_fits[[i]]$BUGSoutput$sims.list$lp_age,
                 mod_fits[[i]]$BUGSoutput$sims.list$lp_esc)
  WAIC[i] <- waic(ldens)$waic
  LOO[i] <- -2*waic(ldens)$elpd_loo
 
}
d_WAIC <- round(WAIC-min(WAIC),1)
wt_WAIC <- exp(-0.5*d_WAIC)/sum(exp(-0.5*d_WAIC))
d_LOO <- round(LOO-min(LOO),1)
wt_LOO <- exp(-0.5*LOO)/sum(exp(-0.5*LOO))
tbl_WAIC <- data.frame(errors=c(2,2,3,3),
                       model=rep(c("Ricker","B-H"),n_mods/2),
                       WAIC=round(WAIC,1),
                       d_WAIC=d_WAIC,
                       wt_WAIC=round(wt_WAIC,3),
                       LOO=round(LOO,1),
                       d_LOO=d_LOO,
                       wt_LOO=round(wt_LOO,3))
 #write.csv(tbl_WAIC,  row.names = FALSE,
#         file = "Skagit_SFChin_SR_model_selection_all.csv")
best_idx <- which(tbl_WAIC$d_WAIC==0)
tbl_WAIC[order(tbl_WAIC$WAIC, tbl_WAIC$LOO),]

## ----refit_best_mod------------------------------------------------------
mod_jags$parameters.to.save <- c("alpha","E_Rkr_a","mu_Rkr_a","beta",
                                 "lp_age","lp_esc",
                                 "Sp","Rec","tot_ln_Rec","ln_RS","pi_vec",
                                 "sigma_r","sigma_s","res_ln_Rec")
mod_jags$inits <- init_vals_AR
mod_jags$model.file <- "IPM_RK_AR.txt"
best_fit <- do.call(jags.parallel, mod_jags)

## load model fits
savedir <- here("Skagit_Chinook_Ricker_SR_SSM_SFCK_cache/latex")
lazyLoad(sub(".rdb","",grep("^JAGS_IO.*.rdb",savedir,value=TRUE)))

library("R2jags")

attach.jags(mod_fit)

## load saved workspace
orig_dir <- getwd()

work_dir <- "/analysis/Skagit_Chinook_SR_SSM_SFCK_cache/latex"

paste0(orig_dir,"/",work_dir)

setwd(paste0(orig_dir,"/",work_dir))

lazyLoad(sub(".rdb","",grep("^JAGS_IO.*.rdb",dir(),value=TRUE)))

library("R2jags")

attach.jags(mod_fit)


## ----model_diagnostics---------------------------------------------------
mod_fit <- best_fit
## Rhat values for all parameters
rh <- mod_fit$BUGSoutput$summary[,"Rhat"]
## histogram of Rhat values for all parameters
par(mai=c(0.9,0.9,0.3,0.1))
hist(rh, breaks=seq(1,ceiling(max(rh)/0.01)*0.01,by=0.0001),main="",
     col=rgb(0,0,255,alpha=50,maxColorValue=255),border="blue3",xlab=expression(italic(R[hat])))
## Rhat values > threshold
bad_Rhat <- rh[rh>Rhat_thresh]
## prop of params with Rhat > threshold
round(length(bad_Rhat)/length(rh),3)
## param names
par_names <- sub("\\[.*","",names(bad_Rhat))
## number of Rhat > threshold by param name
table(par_names)
## index values for offenders
idx <- as.integer(sub("(^.*\\[)([0-9]{1,3})(.*)","\\2",names(bad_Rhat)))
## data frame of offenders
(df <- data.frame(par=par_names, index=idx))


## ----tbl_summary_stats---------------------------------------------------
tbl_smry <- mod_fit$BUGSoutput$summary[c("alpha","E_Rkr_a","beta"),
                                       c("mean","sd","2.5%","50%","97.5%")]
print(tbl_smry,digits=3,quote=FALSE,justify="right")

## ----plot_S_R, fig.width=8, fig.height=5, fig.pos="placeHere"------------
layout(matrix(c(1,1,2,3),2,2),c(3,2),c(1,1))
CI_vec <- c(0.025,0.5,0.975)
offSet <- 0.06

MC <- 100
set.seed(123)
idx <- sample(seq(mcmc_samp),MC)

## posterior of spawners
sDat <- apply(mod_fit$BUGSoutput$sims.list$Sp,2,quantile,CI_vec)
sDat <- sDat[,1:(n_yrs-age_min)]
## posterior of recruits
rDat <- exp(apply(mod_fit$BUGSoutput$sims.list$tot_ln_Rec,2,quantile,CI_vec))
## median values for a & b
aa <- apply(mod_fit$BUGSoutput$sims.list$mu_Rkr_a,1,median)
bb <- mod_fit$BUGSoutput$sims.list$beta
# aa <- median(mod_fit$BUGSoutput$sims.list$alpha)
## empty plot space for spawner-recruit relationships
dd <- 3000
yM <- Re2prec(max(rDat),"ceiling",dd)
#yM <- 30000
xM <- Re2prec(max(sDat),"ceiling",dd)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0,0,0))
plot(sDat[2,],rDat[2,], xlim=c(0,xM), ylim=c(0,yM), pch=16, col="blue3", type="n",
	 xaxs="i", yaxs="i", ylab="Recruits (1000s)", xlab="Spawners (1000s)", cex.lab=1.2,
	 xaxt="n", yaxt="n")
axis(1, at=seq(0,xM,dd*2), labels=seq(0,xM,dd*2)/1000)
axis(2, at=seq(0,yM,dd*2), labels=seq(0,yM,dd*2)/1000)
for(i in 1:MC) { lines((seq(xM)*exp(aa[idx[i]]-bb[idx[i]]*seq(xM))), col="darkgray") }
# lines(aa*seq(0,xM)/(1+bb*seq(0,xM)), col="darkgray")
## add S-R estimates and medians
abline(a=0,b=1,lty="dashed")
nCB <- n_yrs-age_max
points(sDat[2,1:nCB],rDat[2,1:nCB], xlim=c(0,xM), ylim=c(0,yM), pch=16, col="blue3")
segments(sDat[2,1:nCB],rDat[1,1:nCB],sDat[2,1:nCB],rDat[3,1:nCB], col="blue3")
segments(sDat[1,1:nCB],rDat[2,1:nCB],sDat[3,1:nCB],rDat[2,1:nCB], col="blue3")
nTB <- dim(sDat)[2]
clr <- rgb(100, 0, 200, alpha=seq(200,100,length.out=age_max-age_min), maxColorValue=255)
segments(sDat[2,(nCB+1):nTB],rDat[1,(nCB+1):nTB],sDat[2,(nCB+1):nTB],rDat[3,(nCB+1):nTB], col=clr)
segments(sDat[1,(nCB+1):nTB],rDat[2,(nCB+1):nTB],sDat[3,(nCB+1):nTB],rDat[2,(nCB+1):nTB], col=clr)
points(sDat[2,(nCB+1):nTB],rDat[2,(nCB+1):nTB],
       xlim=c(0,xM), ylim=c(0,yM), pch=16, col=clr)
text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
	 y=par()$usr[4]-offSet*diff(par()$usr[3:4]),"(a)")

## posterior for alpha
clr <- rgb(0, 0, 255, alpha = 50, maxColorValue = 255)
a_thresh <- 99
par(mai=c(0.8,0.4,0.3,0.1))
## Ricker alpha
R_alpha_est <- mod_fit$BUGSoutput$sims.list$alpha
alphaCI <- quantile(R_alpha_est,c(0.025,0.5,0.975))
R_alpha_est[R_alpha_est>a_thresh] <- a_thresh
hist(R_alpha_est,freq=FALSE,xlab="",main="",breaks=seq(0,10,0.2),
     col=clr, border="blue3", ylab="", cex.lab=1.2, yaxt="n")
aHt <- (par()$usr[4]-par()$usr[3])/12
arrows(alphaCI,par()$usr[3],alphaCI,par()$usr[3]-aHt,
       code=1,length=0.05,xpd=NA,col="blue3",lwd=1.5)
mtext(expression(Instrinsic~productivity~(alpha)), 1, line=3, cex=1)
text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
	 y=par()$usr[4]-offSet*diff(par()$usr[3:4]),"(b)")

## posterior for K
par(mai=c(0.8,0.4,0.3,0.1))
aa <- matrix(mod_fit$BUGSoutput$sims.array[,,"E_Rkr_a"],ncol=1)
bb <- matrix(mod_fit$BUGSoutput$sims.array[,,"beta"],ncol=1)
R_b_est <- (aa)/bb
R_b_est <- R_b_est[R_b_est > 0]
R_b_CI <- quantile(R_b_est,c(0.025,0.5,0.975))
R_b_est[R_b_est>1e5] <- 1e5
brks <- seq(Re2prec(min(R_b_est),"floor",2000),1e5,length.out=length(seq(0,9,0.2)))
hist(R_b_est, freq=FALSE, breaks=brks, col=clr, border="blue3",
	 xlab="", xaxt="n", yaxt="n",
	 main="", ylab="", cex.lab=1.2)
axis(1, at=seq(Re2prec(min(R_b_est),"floor",2000),1e5,20000))
aHt <- (par()$usr[4]-par()$usr[3])/12
arrows(R_b_CI,par()$usr[3],R_b_CI,par()$usr[3]-aHt,
       code=1,length=0.05,xpd=NA,col="blue3",lwd=1.5)
mtext(expression(Carrying~capacity~(italic(K))), 1, line=3, cex=1)
text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
	 y=par()$usr[4]-offSet*diff(par()$usr[3:4]),"(c)")

## ----ref_pts-------------------------------------------------------------
# abbreviations for ref points
ref_names <- c("MSY","Smsy","Umsy","Umax","Seq","Scrit")
# proportions of MSY to consider
yld_prop <- c(0.75,0.85,0.95)
aa <- matrix(mod_fit$BUGSoutput$sims.array[,,"E_Rkr_a"],ncol=1)
alpha <- matrix(mod_fit$BUGSoutput$sims.array[,,"alpha"],ncol=1)
mcmc <- length(aa)
# empty matrix for ref pts
ref.pts <- matrix(NA,mcmc,length(ref_names))
colnames(ref.pts) <- ref_names
# spawner series for optimal yield profile
SS <- seq(100,3e4,100)
# empty matrix for optimal yield profiles
OYP <- matrix(0,length(SS),length(yld_prop))
for(i in 1:mcmc) {
  # spawners at MSY
  ref.pts[i,"Smsy"] <- (1 - lambert_W0(exp(1-aa[i]))) / bb[i]
  # MSY
  ref.pts[i,"MSY"] <- ref.pts[i,"Smsy"]*((exp(aa[i]-bb[i]*ref.pts[i,"Smsy"])) - 1)
  # harvest rate at MSY
  ref.pts[i,"Umsy"] <- (1 - lambert_W0(exp(1-aa[i])))
  # max harvest rate
  ref.pts[i,"Umax"] <- 1 - 1/alpha[i]
  # equilibrium escapement
  ref.pts[i,"Seq"] <- aa[i]/bb[i]
  # critical escapement
  ref.pts[i,"Scrit"] <- .05*ref.pts[i,"Seq"]
  
  # yield over varying S
  yield <- SS*(exp(aa[i]-bb[i]*SS) - 1)
  for(j in 1:length(yld_prop)) {
    OYP[,j] <- OYP[,j] + 1*(yield > yld_prop[j]*ref.pts[i,"MSY"])
  }
}
OYP <- OYP/mcmc


## Prob of overfishing
hh <- seq(100)
Pr_over <- cbind(hh,hh,hh)
colnames(Pr_over) <- c("Umsy75","Umsy","Umax")
for(i in hh) {
  Pr_over[i,"Umsy75"] <- sum(ref.pts[,"Umsy"]*0.75 < i/100)/mcmc_samp
  Pr_over[i,"Umsy"] <- sum(ref.pts[,"Umsy"] < i/100)/mcmc_samp
  Pr_over[i,"Umax"] <- sum(ref.pts[,"Umax"] < i/100)/mcmc_samp
}

## observed exploitation rate & posterior spawner abundance
Sp_ts <- mod_fit$BUGSoutput$sims.list$Sp[,1:n_yrs]


## ----plot_ref_pts, fig.width=8, fig.height=4, warning=FALSE--------------

layout(matrix(c(2,1,4,3),2,2),heights=c(1,5))

## OYP
par(mai=c(0.9,0.9,0,0), omi=c(0,0,0.1,0.1))
x_lp <- yld_prop
for(i in 1:length(x_lp)) {
  x_lp[i] <- SS[max(which(OYP[,i] == max(OYP[,i]) | abs(OYP[,i] - (yld_prop[i]-0.3)) <= 0.05))]
}
matplot(SS, OYP, type="l", lty="solid", las=1, col=c("slateblue","blue","darkblue"), lwd=2,
        xlab="Spawners", ylab="Probability of X% of MSY", cex.lab=1.2,
        main="", ylim=c(0,1))
points(x=x_lp, y=yld_prop-0.3, pch=21, cex=3.5, col="white", bg="white")
text(x=x_lp, y=yld_prop-0.3, paste0(yld_prop*100,"%"),
     col=c("slateblue","blue","darkblue"), cex=0.7)
text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
     y=par()$usr[4]-offSet*diff(par()$usr[3:4]),"(a)")
## posterior spawner abundance over all years
par(mai=c(0,0.9,0.05,0))
hist(Sp_ts[Sp_ts<3e4], col=clr, border="blue3", breaks=40,
     main="", yaxs="i", xaxt="n",yaxt="n",ylab="")

## prob of overfishing
par(mai=c(0.9,0.9,0,0))
matplot(Pr_over, type="l", las=1, lwd=2, lty="solid", col=c("slateblue","blue","darkblue"),
        ylab="Probability of overfishing", cex.lab=1.2,
        xlab="Harvest rate", xaxt="n")
axis(1,seq(0,100,20),seq(0,100,20)/100)
x_lp <- c(0,0,0)
for(i in 1:length(x_lp)) {
  x_lp[i] <- max(which(abs(Pr_over[,i] - 0.5) <= 0.05))
}
points(x=x_lp, y=rep(0.5,3), pch=21, cex=4, col="white", bg="white")
text(x=x_lp, y=0.5, expression(U[M75], U[MSY], U[Max]),
     col=c("slateblue","blue","darkblue"), cex=0.8)
text(x=par()$usr[1]+par()$pin[2]/par()$pin[1]*offSet*diff(par()$usr[1:2]),
     y=par()$usr[4]-offSet*diff(par()$usr[3:4]),"(b)")



## ----table_ref_pts-------------------------------------------------------

tbl_refpt_smry <- apply(ref.pts,2,quantile,CI_vec)
print(tbl_refpt_smry,digits=3,quote=FALSE,justify="right")

## ----plot_escapement, fig.width=6, fig.height=4.5, fig.pos="placeHere"----
pDat <- apply(mod_fit$BUGSoutput$sims.list$Sp,2,quantile,CI_vec)
ypMin <- min(pDat[,1:n_yrs])
ypMax <- max(pDat[,1:n_yrs])
t_idx_T <- seq(yr_frst,length.out=n_yrs)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(t_idx_T,pDat[3,1:n_yrs], ylim=c(ypMin,ypMax), type="n", log="y", xaxt="n", yaxt="n",
     xlab="Year", ylab="Spawners", main="", cex.lab=1.2)
polygon(c(t_idx_T,rev(t_idx_T)),c(pDat[3,1:n_yrs],rev(pDat[1,1:n_yrs])), col=clr, border=NA)
lines(t_idx_T, pDat[2,1:n_yrs], col="blue3", lwd=2)
points(seq(yr_frst,length.out=n_yrs), exp(ln_dat_esc), pch=16, cex=1)
axis(1,at=seq(1980,2015,5))
axis(2,at=c(5000,9000,20000))

## ----plot_recruits, fig.width=6, fig.height=4.5, fig.pos="placeHere"-----

CI_vec <- c(0.025,0.5,0.975)
pDat <- apply(mod_fit$BUGSoutput$sims.list$Rec,c(1,2),sum)
pDat <- apply(apply(pDat,2,sort),2,function(x) { x[mcmc_samp*CI_vec] })
ypMin <- min(pDat)
ypMax <- max(pDat)
t_idx_a <- seq(yr_frst,length.out=n_yrs-age_min)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(t_idx_a,pDat[3,], ylim=c(ypMin,ypMax), type="n", log="y", yaxt="n",
     xlab="Brood year", ylab="Recruits", main="", cex.lab=1.2)
axis(2,at=c(9000,15000,27000))
polygon(c(t_idx_a,rev(t_idx_a)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(t_idx_a, pDat[2,], col="blue3", lwd=2)



## ----plot_R_per_S, fig.width=6, fig.height=4.5, fig.pos="placeHere"------

pDat <- apply(mod_fit$BUGSoutput$sims.list$ln_RS,2,sort)
pDat <- apply(pDat,2,function(x) { x[mcmc_samp*CI_vec] })
pDat[2,] <- apply(mod_fit$BUGSoutput$sims.list$ln_RS,2,median)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(t_idx_a,pDat[3,], ylim=c(ypMin,ypMax), type="n", #log="y",
     xlab="Brood year", ylab="ln(R/S)", main="", cex.lab=1.2)
abline(h=0, lty="dashed")
polygon(c(t_idx_a,rev(t_idx_a)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(t_idx_a, pDat[2,], col="blue3", lwd=2)


## ----plot_age_comp, fig.width=6, fig.height=4, fig.pos="placeHere"-------

par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.1,0.2,0.2))
clr <- rgb(0, 0, 255, alpha = 40, maxColorValue = 255)
age_est <- t(apply(apply(mod_fit$BUGSoutput$sims.list$pi_vec,c(3,2),mean),2,cumsum))
nRec <- n_yrs-age_min
plot(t_idx_a, rep(1,nRec), ylab="Proportion", xlab="Brood year", ylim=c(0,1), las=1,
     xaxs="i", yaxs="i", type="n", lty="solid", col="blue3", cex.lab=1.2)
for(i in c(1,2,3,4)) {
  polygon(c(t_idx_a,rev(t_idx_a)),c(age_est[,i],rep(0,nRec)), col=clr, border=NA)
}
lbl <- apply(cbind(c(0,age_est[nRec,-A]),age_est[nRec,]),1,mean)
text(par()$usr[2],par()$usr[4]*1.05,"Age", xpd=NA, pos=4, offset=0.05, col="black", cex=0.8)
text(par()$usr[2],lbl[1:4],seq(2,5), xpd=NA, pos=4, col="black", cex=0.7)


## ----plot_innovations, fig.width=6, fig.height=4, fig.pos="placeHere"----
t_idx_a <- seq(yr_frst,length.out=n_yrs-age_min)
pDat <- apply(mod_fit$BUGSoutput$sims.list$res_ln_Rec,2,quantile,CI_vec)
ypMin <- min(pDat)
ypMax <- max(pDat)
par(mai=c(0.8,0.8,0.1,0.1), omi=c(0,0.2,0.1,0.2))
plot(t_idx_a,pDat[3,], ylim=c(ypMin,ypMax), type="n", #log="y",
	 xlab="Brood year", ylab="Innovations", main="", cex.lab=1.2)
abline(h=0, lty="dashed")
polygon(c(t_idx_a,rev(t_idx_a)),c(pDat[3,],rev(pDat[1,])), col=clr, border=NA)
lines(t_idx_a, pDat[2,], col="blue3", lwd=2)

## ----VRAP_sim------------------------------------------------------------

## ----VRAP_sim------------------------------------------------------------

## posterior S/R parameters
aa <- mod_fit$BUGSoutput$sims.list$E_Rkr_a
bb <- mod_fit$BUGSoutput$sims.list$beta

## management error mean and std dev. (FRAM post season ER/pre- season ER)
me_mean <- 1.07
me_sd <- 0.19

## bound observed ER by a minimum and maximum to avoid unrealistic harvest scenarios
ER_max <- 0.95
ER_min <- 0.15 

##select whether to apply management error or not (1 = include, 0 = exclude)
me <- 0

## number of years for each simulation: starting with 25 since this is what NOAA requires
numYears <- 25

## number of simulations
numSim <- 1000

## range of target harvest rates to explore
ER_range <- seq(0,0.8,.01)

## start with most recent five spawner years in time series
S_start <- apply(mod_fit$BUGSoutput$sims.list$Sp,2,quantile,probs=CI_vec)[2,(n_yrs - 4):n_yrs]  

## SMSY under average conditions
Smsy <- quantile(ref.pts[,"Smsy"],CI_vec) 
Scrit <- quantile(ref.pts[,"Scrit"],CI_vec) 

## maturation schedule
matSched_2 <- apply(mod_fit$BUGSoutput$sims.list$pi_vec[,,1],2,quantile,CI_vec)
matSched_3 <- apply(mod_fit$BUGSoutput$sims.list$pi_vec[,,2],2,quantile,CI_vec)
matSched_4 <- apply(mod_fit$BUGSoutput$sims.list$pi_vec[,,3],2,quantile,CI_vec)
matSched_5 <- apply(mod_fit$BUGSoutput$sims.list$pi_vec[,,4],2,quantile,CI_vec)

matSched <- rbind(matSched_2[2,],matSched_3[2,],matSched_4[2,],matSched_5[2,])


## take 100 random samples of S/R parameters from posterior distribution
sample <- 100
samplePosterior <- sample(1:mcmc,sample,replace = TRUE)

## these are the posterior samples that will be used. Note that for each paired sample,
## 1000 25 year simulations will be conducted for each target exploitation rate
a <- aa[samplePosterior]
b <- bb[samplePosterior]


## median posterior process standard deviation
median_sd_r <- median(mod_fit$BUGSoutput$sims.list$sigma_r)^0.5

## output arrays for RER criteria
p_UET <- array(data = NA,dim = c(length(ER_range),(sample + 1)))
p_UET_10pct <- array(data = NA,dim = c(length(ER_range),(sample + 1)))
p_Crit_5pct <- array(data = NA,dim = c(length(ER_range),(sample + 1)))

for(i in 1:(sample + 1)){

  ## set up arrays for keeping track of abundance
  sp<-array(data = NA, dim = c(numSim,(numYears+5),length(ER_range)))
  catch<-array(data = NA, dim = c(numSim,(numYears+5),length(ER_range)))
  mature_Run<-array(data = NA, dim = c(numSim,(numYears+5),length(ER_range)))
  
  ## SR parameters used for this set of simulations.Median is evaluated first.
  if(i == 1){a_Sim <- median(aa);b_Sim <- median(bb)}else{
    a_Sim <- a[i-1];b_Sim <- b[i-1]
  }
  
  ## reset counter for each round of simulations
  c <- 1
  for (ER_target in ER_range){
    for (sim in 1:numSim){
      #sim <- 2
      ## output vector for total age specific recruitment
      age2Rec <- NULL
      age3Rec <- NULL
      age4Rec <- NULL
      age5Rec <- NULL
      
      ## sample from yearly estimates of maturation schedule
      matSchedule_sim <- matSched[,sample(1:dim(matSched)[2],numYears + 5, replace = TRUE)]
      
      if(me == 1){
        
        ## apply management error to target exploitation rates within bounds of
        ## 15% - 95% total ER
        mgmtError <- rnorm(numYears + 5,me_mean,me_sd)
        ER_Obs <- ER_target*mgmtError
        ER_Obs[which(ER_Obs > ER_max)] <- ER_max
        ER_Obs[which(ER_Obs < ER_min)] <- ER_min
          
      }else{ER_Obs <- rep(ER_target,numYears + 5)}
      
      
      for(year in 1:(numYears + 5)){
        
        if(year <= 5){
          
          lnRec <- (log(S_start[year]) + a_Sim- b_Sim*S_start[year])
          
          totRec <- exp(rnorm(1,lnRec,median_sd_r))
          
          ## apply maturation schedule to project recruits forward
          totRec <- totRec * matSchedule_sim[,year]
          age2Rec[year+2]<-totRec[1]
          age3Rec[year+3]<-totRec[2]
          age4Rec[year+4]<-totRec[3]
          age5Rec[year+5]<-totRec[4]
          
        }
        
        if(year > 5){
          #year <- 6
          mature_Run[sim,year,c] <- age2Rec[year] + age3Rec[year] + age4Rec[year] + age5Rec[year]
          catch[sim,year,c]<-ER_Obs[year]*(mature_Run[sim,year,c])
          sp[sim,year,c] <- mature_Run[sim,year,c] - catch[sim,year,c]
          
          lnRec <- (log(sp[sim,year,c]) + a_Sim - b_Sim*sp[sim,year,c])
          
          totRec <- exp(rnorm(1,lnRec,median_sd_r))
          
          ## apply maturation schedule to project recruits forward
          totRec <- totRec * matSchedule_sim[,year]
          age2Rec[year+2]<-totRec[1]
          age3Rec[year+3]<-totRec[2]
          age4Rec[year+4]<-totRec[3]
          age5Rec[year+5]<-totRec[4]
          
        }
        
      }#next time step
    }# next sim
    c <- c+1
  }#next Er
  
  p_UET[,i] <- apply(sp[,27:30,],3,function(x){length(which(x > Smsy[2]))/length(x)})
  p_UET_10pct[,i] <- p_UET[1,i] - p_UET[,i]

  p_Crit_temp <- apply(sp[,5:30,],3,function(x){length(which(x < Scrit[2]))/length(x)})
  p_Crit_5pct[,i] <- p_Crit_temp - p_Crit_temp[1]
  
}## next sample

## compute median RER based on criterion
RER_UET_median <- ER_range[length(which(p_UET[,1] >= 0.80))]
RER_UET_10pct_median <- ER_range[length(which(p_UET_10pct[,1] <= 0.10))]
RER_Crit_5pct_median <- ER_range[length(which(p_Crit_5pct[,1] <= 0.05))]



## ----plot_RER_profile, fig.width=8, fig.height=4, fig.pos="placeHere"----

## plot RER profiles for each criterion
layout(matrix(seq(1,3,1),1,3),heights = 4,widths = c(4,4,4,4)) 
par(mai=c(0.15,0.15,0.15,0.15), omi=c(0.5,0.5,0.5,0.5))

plot(p_Crit_5pct[,1]~ER_range,xlab = "",ylab = "",type = "l", bty = "n",lwd = 2,main = "% > LET.base - % < LET")
points(RER_Crit_5pct_median,0.05,pch = 16,col = "red",cex = 1.5)
abline(h = 0.05,v =RER_Crit_5pct_median,lty = 2, col = "grey" )

plot(p_UET[,1]~ER_range,xlab = "",ylab = "",type = "l",bty = "n",lwd = 2, main = "% > UET")
points(RER_UET_median,0.80,pch = 16,col = "red",cex = 1.5)
abline(h = 0.80,v = RER_UET_median,lty = 2, col = "grey")

plot(p_UET_10pct[,1]~ER_range,xlab = "",ylab = "",type = "l",bty = "n", lwd = 2, main = "% > UET.base - % > UET")
points(RER_UET_10pct_median,0.10,pch = 16,col = "red",cex = 1.5)
abline(h = 0.10,v = RER_UET_10pct_median,lty = 2, col = "grey")

mtext("Target exploitation rate",1,outer = TRUE,cex = 1.2,line = 2)


## ----RER_posterior_dist--------------------------------------------------

##calculate posterior RER distribution
RER_UET <- NULL
RER_UET_10pct <- NULL
RER_Crit_5pct <- NULL

for(i in 2:dim(p_UET)[2]){
  RER_UET[i] <- ER_range[length(which(p_UET[,i] >= 0.80))]
  RER_UET_10pct[i] <- ER_range[length(which(p_UET_10pct[,i] <= 0.10))]
  RER_Crit_5pct[i] <- ER_range[length(which(p_Crit_5pct[,i] <= 0.05))]
  
}


## ----plot_RER_posterior_dist, fig.width=8, fig.height=4, fig.pos="placeHere"----

## plot posterior RER's
layout(matrix(seq(1,3,1),1,3),heights = 4,widths = c(4,4,4,4)) 
par(mai=c(0.15,0.15,0.15,0.15), omi=c(0.5,0.5,0.5,0.5))

hist(RER_Crit_5pct,col = "grey",main = "RER(<5%)", xlab = "", ylab = "", ylim = c(0,40))
abline(v = RER_Crit_5pct_median,lwd = 2)

hist(RER_UET,col = "grey", main = "RER(>80%)",xlab = "", ylab = "")
abline(v = RER_UET_median,lwd = 2)

hist(RER_UET_10pct,col = "grey", main = "RER(<10%",xlab = "",ylab = "")
abline(v = RER_UET_10pct_median,lwd = 2)

mtext("Target exploitation rate",1,outer = TRUE,cex = 1.2,line = 2)

## ----tbl_RER_CI----------------------------------------------------------
RER_range <- quantile(RER_UET,c(0.025,0.975),na.rm = TRUE)
RER_range <- c(RER_range[1],RER_UET_median,RER_range[2])
#RER_range <- round(RER_range,2)

## ----summary_tbl
print(rbind(Scrit = round(Scrit),Smsy = round(Smsy),RER_range = RER_range))

