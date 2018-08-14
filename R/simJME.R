#' Function to simulate a data structure suitable for joint modeling!
#'
#' @param n number of individuals
#' @param n_i initial number of longitudinal measurements per individual
#' @param betal coefficient vector for the longitudinal model
#' @param betas coefficient vector for the survival model
#' @param betals coefficient vector for the shared model
#' @param betat value for the time effect
#' @param int value for the intercept
#' @param alpha value for the association parameter
#' @param lambda value for the baseline hazard
#' @param sigma2 value for the model error
#' @param noninfi number of additional, non-informative covariates per predictor

simJM <- function(n, n_i, betal=0, betas=0, betat=0, betals=0,
                  int, alpha, lambda, sigma2, noninfl = 0, noninfs = 0, noninfls = 0, high.dim = FALSE){

  ### generate id vector
  id = rep(1:n, each = n_i)
  first = rep(c(1, rep(0,n_i-1)), n)

  ### generate covariate matrices, standardize columns to mean = 0, sd = 1
  ### note, that Xls contains duplicated values of the same measurement, Xls_un not
  if(is.null(betal)){Xl = 0}else{
    pl = length(betal)
    Xl = matrix(0, n*n_i, pl)
    for(i in 1:pl){
      Xl[,i] = runif(n*n_i, -1, 1)
      Xl[,i] = Xl[,i]/sd(Xl[,i])
    }
  }
  if(is.null(betas)){Xs = 0}else{
    ps = length(betas)
    Xs = matrix(0, n, ps)
    for(i in 1:ps){
      Xs[,i] = runif(n, -1, 1)
      Xs[,i] = Xs[,i]/sd(Xs[,i])
    }
  }
  if(is.null(betals)){Xls = 0; Xls_un = 0}else{
    pls = length(betals)
    Xls = matrix(0, n*n_i, pls)
    for(i in 1:pls){
      Xls[,i] = rep(runif(n, -1, 1), each = n_i)
      Xls[,i] = Xls[,i]/sd(Xls[,i])
    }
    Xls_un = Xls[first==1,]
  }

  ### generate time points
  day_in_year = sample(1:365, n*n_i, replace = TRUE)
  time = rep(seq(0,(n_i-1)*365, 365), n)
  time = time + day_in_year
  for(i in 1:n){time[id==i] = time[id==i] - min(time[id==i])}
  T_long = time/(n_i*365)

  ### generate random effects and compute predictor vectors
  gamma0 = rnorm(n, 0, .01)
  gamma1 = rnorm(n, 0, .01)
  etal = int + Xl%*%betal
  etas = Xs%*%betas
  etals = Xls%*%betals + rep(gamma0, each=n_i) + T_long*(rep(gamma1, each=n_i) + betat)

  ### simulate longitudinal outcome
  y = rnorm(n*n_i, etal + etals, sqrt(sigma2))

  ###### simulate event times with censoring via inversion sampling
  ### define shortcuts for time-dependent and time-independent part of etals times alpha
  a.etals.ti = alpha*(gamma0 + Xls_un%*%betals)
  a.etals.td = alpha*(gamma1 + betat)

  ### draw U([0,1]) random numbers and plug them in the inverse distribution function
  u = runif(n)
  T_surv = log(-log(1-u)*a.etals.td/(lambda*exp(etas + a.etals.ti)) + 1)/a.etals.td
  T_surv[is.nan(T_surv)] = 2

  ### create censoring index, censoring occurs if event-time exceeds last measurement time
  time_mat <- matrix(nrow = n_i, data = T_long)
  delta = rep(1, n)
  for(i in 1:n){
    if(T_surv[i]>max(time_mat[,i])){T_surv[i] = max(time_mat[,i]); delta[i] = 0}
    else if(which.max(time_mat[,i]>T_surv[i])<=n_i){time_mat[which.max(time_mat[,i]>T_surv[i]):n_i,i] = 42}
  }

  ### remove all values corresponding to measurement times after an event has occured
  time_zero = which(as.vector(time_mat)==42)

  if (length(time_zero == 0) > 0) {
    id = id[-time_zero]
    y = y[-time_zero]
    Xl = Xl[-time_zero,]
    Xls = Xls[-time_zero,]
    T_long = T_long[-time_zero]
  }

  if(high.dim){
    n_hd = length(id)
    p_hd = length(betal) + length(betas) + length(betals) + 1 + as.numeric(betat!=0)
    r_hd = ceiling((n_hd - p_hd)/3)
    noninfl = r_hd
    noninfs = r_hd
    noninfls = r_hd
  }

  ### add non-informative covariates, if selected
  if(noninfl > 0 | noninfls>0 | noninfs>0){
    nni = length(id)
    for(i in 1:max(c(noninfl, noninfls, noninfs))){
      if(i <= noninfl){
        x = runif(nni, -1, 1); x = x/sd(x)
        Xl = cbind(Xl, x)}
      if(i <= noninfs){
        x = runif(n, -1, 1); x = x/sd(x)
        Xs = cbind(Xs, x)}
      if(i <= noninfls){
        x = c()
        for(j in 1:n){
          x = c(x, rep(runif(1, -1, 1), sum(id==j)))
        }
        x = x/sd(x)
        Xls = cbind(Xls, x)}
    }
  }

  return(list(
    #### longitudinal outcome
    "y" = y,
    #### longitudinal predictor fixed effect covariates
    "Xl" = Xl,
    #### survival predictor fixed effect covariates
    "Xs" = Xs,
    #### shared predictor fixed effect covariates
    "Xls" = Xls,
    ### random intercept
    "gamma0" = gamma0,
    ### random slope
    "gamma1" = gamma1,
    #### id indicator
    "id" = id,
    #### measurement times
    "T_long" = T_long,
    #### Event times
    "T_surv" = T_surv,
    #### Censoring indicator
    "delta" = delta))
}
