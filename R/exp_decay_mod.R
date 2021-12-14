library(Matrix)
exp_decay_mod <- function(params_ini=c(0, 0.001), sigma, t, u){
  n <- length(t)
  sigma <- sigma
  #t-u vector #Xi(ti) vector
  compare_ti_ui <- (t-u)
  Xi_ti <- compare_ti_ui
  Xi_ti[Xi_ti <= 0] <- 0
  Xi_ti[Xi_ti > 0 ] <- 1

  #compare ti_uj matrix
  #Xj_ti matrix
  compare_ti_uj <- outer(t,u,FUN=function(x,y)(x-y)) #n*n matrix
  Xj_ti <- compare_ti_uj
  Xj_ti[Xj_ti <= 0] <- 0
  Xj_ti[Xj_ti > 0 ] <- 1

  #Yj_ti matrix
  Yj_ti <- outer(t, t, FUN=function(x,y)(x-y))
  Yj_ti[Yj_ti < 0] <- 0
  Yj_ti[Yj_ti > 0] <- 1
  diag(Yj_ti) <- ifelse(sigma==0, 1, 0)

  ###############
  #making Xj_ti and Yj_ti into sparse matrices
  ###############
  Xj_ti = as(Xj_ti,"RsparseMatrix")
  Yj_ti = as(Yj_ti,"RsparseMatrix")

  ####################################
  #logLikelihood
  logL <- function(params){
    #n <- length(t) defined before the three functions
    b <- params[1]
    r <- params[2]

    Zj_ti <- Xj_ti*exp(-r * (compare_ti_uj))
    Zj_ti = as(Zj_ti,"RsparseMatrix")

    #first part of log Likelihood
    #first_part <- sum( sigma * b * Xi_ti * exp(-r*compare_ti_ui) )
    first_part <- sum( sigma * b * diag(Zj_ti) )

    #second part of log Likelihood
    inner_part <- colSums(Yj_ti * exp(b*Zj_ti))
    log_part <- log(inner_part)
    #print(inner_part)

    second_part <- sum(sigma*log_part)
    #print(second_part)

    return (-first_part + second_part)
    #return (first_part - second_part)
  }

  ######################################
  #gradient with respect to gamma
  g_gamma = function(params){
    #n <- length(t)
    b <- params[1]
    r <- params[2]

    Zj_ti <- Xj_ti*exp(-r * (compare_ti_uj))
    Zj_ti = as(Zj_ti,"RsparseMatrix")

    first_part <- (-1) * sum(sigma * b * diag(Zj_ti) * compare_ti_ui )

    num <- colSums(Yj_ti * Xj_ti * exp(b * Zj_ti - r * compare_ti_uj))
    denom <- colSums(Yj_ti * exp(b * Zj_ti))
    sec_part <- b * sum(sigma* (num / denom))

    return(first_part + sec_part)
  }

  ########################################
  #gradient with respect to beta
  b.derivative = function(params){
    b = params[1]
    r = params[2]

    Zj_ti = exp(-r*compare_ti_uj)*Xj_ti
    Zj_ti = as(Zj_ti,"RsparseMatrix")
    #first part of beta derivative
    A = b*Zj_ti
    first_part = sum(sigma*diag(Zj_ti))

    exp.A = exp(A)
    yexp.A = Yj_ti*exp.A
    #denominator of the second part
    #denom = apply(yexp.A,2,sum)
    denom = colSums(yexp.A)

    #numerator of the second part
    #num = apply(yexp.A*Zj_ti,2,sum)
    num = colSums(yexp.A*Zj_ti)

    #second part of beta derivative
    second_part = sum(sigma*num/denom)

    return (first_part - second_part)

  }

  gradients <- function(params){
    #b <- params[1]
    #r <- params[2]
    return ( c(b.derivative(params), g_gamma(params)) )
  }

  #########apply nlminb optimization
  opt <- nlminb(start = params_ini, objective = logL, gradient = gradients, lower=c(-Inf,0),upper=c(+Inf, 1), control = list(step.min=0.01, step.max=0.05))
  #opt <- nlminb(start = params_ini, objective = logL, gradient = gradients)
  #lower=c(-Inf,0),upper=c(+Inf, 1), control = list(step.min=0.01, step.max=0.05))
  #opt <- nlminb(start = params_ini, objective = logL, lower=0, control=list(trace=TRUE))

  print(opt$par)
  return (opt)


}
