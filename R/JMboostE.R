#' Main function to carry out boosting for a joint model!
#'
#' @param y vector containing the longitudinal outcome
#' @param Xl design matrix for the longitudinal part containing time-varying covariates
#' @param Xs design matrix for the survival part containing one measurement per individual and covariate
#' @param Xls design matrix for the shared part containing time-independent covariates (duplicated values per individual)
#' @param T_long longitudinal time points
#' @param T_surv observed survival times
#' @param delta censoring indicator
#' @param id id vector labeling the longitudinal outcome/design matrices with corresponding individuals
#' @param alpha starting value for the association parameter
#' @param lambda starting value for the baseline hazard
#' @param int starting value for the intercept in the longitudinal submodel
#' @param time.effect logical, whether a shared time.effect shall be included
#' @param mstop_i number of boosting iterations per boosting step
#' @param nyi step length per boosting step

JMboost = function(y, Xl = NULL, Xs = NULL, Xls = NULL, T_long, T_surv, delta, id,
                   alpha = 1, lambda = 1, time.effect = TRUE,
                   mstop_l, mstop_s, mstop_ls, nyl = .1, nys = .3, nyls = .1, nyr = .1){

  mstop = max(mstop_l, mstop_s, mstop_ls)
  n = length(id)
  N = length(unique(id))

  ############### construct random effects ######################
  Xr = matrix(ncol=2*N, nrow = n, data=0)
  unid = order(unique(id))
  id = rep(unid, as.vector(table(id)))
  for(i in 1:N){
    Xr[which(id==as.character(i)),i] = 1
    Xr[which(id==as.character(i)),N+i] = T_long[which(id==as.character(i))]
  }
  XrA = Xr[,1:N]
  lambdaran = mboost:::df2lambda(XrA, 4, weights=1)[2]
  XrAt = t(XrA)
  SA = solve(XrAt%*%XrA + lambdaran*diag(N))%*%XrAt
  XrB = Xr[,-(1:N)]
  lambdaran = mboost:::df2lambda(XrB, 4, weights=1)[2]
  XrBt = t(XrB)
  SB = solve(XrBt%*%XrB + lambdaran*diag(N))%*%XrBt
  ###############################################################

  ### set offset fixed and random intercept according to mle
  offset = lme(y ~ 1, random = ~ 1 | id)
  int = 0 #offset$coefficients$fixed
  gamma0 = rep(0, N)#offset$coefficients$random$id
  gamma1 = rep(0, N)
  sigma2 = offset$sigma^2
  print(sigma2)

  ### set starting values based on chosen sets of covariates
  betal = 0
  betas = 0
  betat = 0
  betals = 0
  if(is.null(Xl)){pl = 0; Xl = 0}else{pl = ncol(Xl); betal = rep(0, pl)}
  if(is.null(Xs)){ps = 0; Xs = 0}else{ps = ncol(Xs); betas = rep(0, ps)}
  if(is.null(Xls)){
    pls = 0; Xls = 0; Xls_un = 0
    }else{
    pls = ncol(Xls)
    first = rep(FALSE, n)
    for(i in 1:n){
      first[which.max(id==i)] = TRUE
    }
    Xls_un = as.matrix(Xls[first==1,])
    betals = rep(0, pls)
  }

  ### define storing matrices/vectors
  GAMMA0 = matrix(0, ncol=mstop, nrow=N)
  GAMMA1 = matrix(0, ncol=mstop, nrow=N)
  BETAL = matrix(0, ncol=mstop, nrow=pl)
  BETAS = matrix(0, ncol=mstop, nrow=ps)
  BETALS = matrix(0, ncol=mstop, nrow=pls)
  BETAT = rep(0, mstop)
  INT = rep(0, mstop)
  ALPHA = rep(0, mstop)
  LAMBDA = rep(0, mstop)
  SIGMA2 = rep(0, mstop)

  for(m in 1:mstop){
    ###############################################################
    #### S1 #######################################################
    ###############################################################
    if(m <= mstop_l){
      etal = as.vector(int + Xl%*%betal)
      etas = as.vector(Xs%*%betas)
      etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
      etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
      ###################COMPUTING THE GRADIENT######################
      u = (y - etal - etals)/sigma2
      ##################/COMPUTING THE GRADIENT######################
      fits = matrix(0, 3, pl + 1 + as.numeric(time.effect))
      if(pl>0){
        for(i in 1:pl){
          fit = mylm(u, Xl[,i])
          fits[1,i] = fit$int
          fits[2,i] = fit$slp
          fits[3,i] = fit$RSS
        }
      }else if(!time.effect){
        int = int + nyl*mean(u)
      }
      rfit = rbind(SA, SB)%*%u
      fits[3, pl+1] = sum((u-(Xr%*%rfit))^2)
      if(time.effect){
        fit = mylm(u, T_long)
        fits[1, pl+2] = fit$int
        fits[2, pl+2] = fit$slp
        fits[3, pl+2] = fit$RSS
      }
      best = which.min(fits[3,])
      if(best==pl+1){
        gamma0 = gamma0 + nyr*rfit[1:N]
        gamma1 = gamma1 + nyr*rfit[-(1:N)]
      }else if(best==pl+2){
        betat = betat + nyl*fits[2,best]
        int = int + nyl*fits[1,best]
      }else{
        betal[best] = betal[best] + nyl*fits[2,best]
        int = int + nyl*fits[1,best]
      }
    }
    INT[m] = int
    BETAT[m] = betat
    BETAL[,m] = betal
    GAMMA0[,m] = gamma0
    GAMMA1[,m] = gamma1

    ###############################################################
    #### S2 #######################################################
    ###############################################################
    if(m<=mstop_s && ps>0){
      etal = as.vector(int + Xl%*%betal)
      etas = as.vector(Xs%*%betas)
      etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
      etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
      ###################COMPUTING THE GRADIENT######################
      if(time.effect){
        if(sum(gamma1)!=0 || betat!=0){
          u = delta - lambda*exp(etas) * (exp(alpha*etals_un) - exp(alpha*(gamma0 + Xls_un%*%betals)))/(alpha*(betat+gamma1))
        }else{
          u = delta - lambda*exp(etas + alpha*etals_un)*T_surv
        }
      }else{
        if(sum(gamma1)!=0){
          u = delta - lambda*exp(etas) * (exp(alpha*etals_un) - exp(alpha*(etals_un - gamma1*T_surv)))/(alpha*gamma1)
        }else{
          u = delta - lambda*exp(etas + alpha*etals_un)*T_surv
        }
      }
      ##################/COMPUTING THE GRADIENT######################
      fits = matrix(0, 3, ps)
      for(i in 1:ps){
        fit = mylm(u, Xs[,i])
        fits[1,i] = fit$int
        fits[2,i] = fit$slp
        fits[3,i] = fit$RSS
      }
      best = which.min(fits[3,])
      betas[best] = betas[best] + nys*fits[2,best]
      lambda = lambda * exp(nys*fits[1,best])
    }
    LAMBDA[m] = lambda
    BETAS[,m] = betas

    ###############################################################
    #### S3 #######################################################
    ###############################################################
    if(m<=mstop_ls && pls>0){
      etal = as.vector(int + Xl%*%betal)
      etas = as.vector(Xs%*%betas)
      etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
      etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)
      ###################COMPUTING THE GRADIENT######################
      u_l = (y - etal - etals)/sigma2
      if(time.effect){
        if(sum(gamma1)!=0 || betat!=0){
          u_s = delta*alpha - lambda*exp(etas)*(exp(alpha*etals_un) - exp(alpha*(gamma0 + Xls_un%*%betals)))/(gamma1 + betat)
        }else{
          u_s = delta*alpha - alpha*lambda*exp(etas)*exp(alpha*etals_un)*T_surv
        }
      }else{
        if(sum(gamma1)!=0){
          u_s = delta*alpha - alpha*lambda*exp(etas + alpha*(Xls_un%*%betals + Xr[,1:N]%*%gamma0))*((exp(gamma1*T_surv) - 1)/gamma1)
        }else{
          u_s = delta*alpha - alpha*exp(etas)*lambda*exp(alpha*etals_un)*T_surv
        }
      }
      u = c(u_l,u_s)
      ##################/COMPUTING THE GRADIENT######################
      fits = matrix(0, 3, pls)
      for(i in 1:pls){
        fit = mylm(u, c(Xls[,i], Xls_un[,i]))
        fits[1,i] = fit$int
        fits[2,i] = fit$slp
        fits[3,i] = fit$RSS
      }
      best = which.min(fits[3,])
      betals[best] = betals[best] + nyls*fits[2,best]
      int = int + nyls*fits[1,best]
      lambda = lambda*exp(nyls*alpha*fits[1,best])
    }
    INT[m] = int
    BETALS[,m] = betals

    ###############################################################
    #### S4 #######################################################
    ###############################################################
    etal = as.vector(int + Xl%*%betal)
    etas = as.vector(Xs%*%betas)
    etals = as.vector(Xr%*%c(gamma0, gamma1) + as.vector(Xls%*%betals) + T_long*betat)
    etals_un = as.vector(gamma0 + gamma1*T_surv + as.vector(Xls_un%*%betals) + betat*T_surv)

    oi.min = min(alpha - 0.1*abs(alpha), alpha - 0.1)
    oi.max = max(alpha + 0.1*abs(alpha), alpha + 0.1)
    optim.int = c(oi.min, oi.max)
    alpha = optimize(surv.risk, lambda=lambda, etals=etals_un, etas=etas, delta=delta, gamma1=gamma1,
                     betals=betals, betat=betat, time.effect=time.effect, T_surv=T_surv, interval=optim.int)$minimum

    ALPHA[m] = alpha
    if(m%%100==0){print(m)}
  }

  p = pl + pls + as.numeric(time.effect) + 3
  sigma2 = sum((y - etal - etals)^2)/(n - p)

  structure(list(GAMMA0 = GAMMA0, GAMMA1 = GAMMA1, BETAL = BETAL, BETAS = BETAS,
                 BETALS = BETALS, ALPHA = ALPHA, LAMBDA = LAMBDA, SIGMA2 = SIGMA2,
                 gamma0 = gamma0, gamma1 = gamma1, betal = betal, betas = betas, betals = betals,
                 betat=betat, int=int, BETAT=BETAT, INT=INT, alpha = alpha, lambda = lambda,  sigma2 = sigma2))
}

surv.risk = function(alpha, lambda, etas, etals_un, T_surv, delta, gamma1, betals, time.effect, betat){
  if(time.effect){
    if(sum(gamma1)!=0||betat!=0){
      integral = lambda*exp(etas)*(exp(alpha*etals_un) - exp(alpha*(etals_un - (gamma1 + betat)*T_surv)))/(alpha*(betat+gamma1))
    }else{
      integral = lambda*exp(etas)*exp(alpha*etals_un)*T_surv
    }
  }else{
    if(sum(gamma1)!=0){
      integral = lambda*exp(etas)*(exp(alpha*etals_un) - exp(alpha*(etals_un - gamma1*T_surv)))/(alpha*gamma1)
    }else{
      integral = lambda*exp(etas)*exp(alpha*etals_un)*T_surv
    }
  }

  risk = -sum(delta*(log(lambda) + etas + alpha*etals_un) - integral)
  return(risk)
}

mylm = function(y,x){
  X = cbind(1,x)
  beta = solve(t(X)%*%X)%*%t(X)%*%y
  RSS = sum((y-X%*%beta)^2)
  return(list("int" = beta[1], "slp" = beta[2], "RSS" = RSS))
}
