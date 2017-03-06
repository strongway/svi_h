library(dplyr)
library(GA) # genetic algo
library(parallel) # parallel processing

#' generic SVI function with parameters: 
#' a,b,m, rho, sigma, and H
svi_fun <- function(k,param){
  a = param[1];
  b = param[2];
  m = param[3];
  rho = param[4];
  sigma = param[5];
  H = param[6];
  result = (a + b*(rho*(k-m) + sqrt((k-m)^2 + sigma^2)))^H
}

#' Estimate svi parameers with GA algorithm
#' 
#' This require ga package, GA gives robust result, though 
#' it is relative slow. Alternative, use constOptim
#' It takes two-step optimization:
#' 1. obtain inner parameters: a,b, rho
#' 2. estimate m and sigma, H using GA
#' @param k log(K/F) -a vector
#' @param iv implied volatility - a vector
#' @param tau maturity (scalar)
#' 
fit_ga_svi <- function(k, iv, tau){
  # total iv
  w = iv*sqrt(tau[1])
  
  # optimization function: x (m, sigma, H)
  ofitfun <- function(x) { # m, sigma, H
    y = (k-x[1])/x[2] # m, sigma
    # find inner parameters: a,b, rho, H
    par <- fit_inner(y, w, x[2],x[3]) # a, b, rho, H
    param = c(par[1],par[2], x[1], par[3], x[2], x[3]) #a, b, m, rho, sigma, H
    raw <- svi_fun(k, param)
    value = 1/mean((raw - w)^2)
    return(value)
  }
  
  # genetic Algorith, set boxes for m, sigma, and H
  op <- GA::ga("real-valued", fitness = ofitfun,
               min = c(min(k),0.00001, 0.4),
               max = c(max(k),10, 1.1),
               maxiter = 80, run =80, optim = TRUE,
               monitor = FALSE)
  # possible multiple solutions, only select the first one
  opar = op@solution[1,]
  
  y = (k-opar[1])/opar[2]
  par <- fit_inner(y, w, opar[2], opar[3])
  # parameters: a,b, m, rho, sigma, H
  parameters = data.frame(tau =tau[1], a=par[1],b=par[2],
                          m =opar[1], rho =par[3], sigma=opar[2],H= opar[3])
  
}

#' Inner optimization for parameter: c,d,a
#'
#' optimization of the inner parameters
#' let y(k) = (k - m)/sig
#' svi w = a + b*sig*(rho*y(k) + sqrt(y(x)^2+1))
#'        = a + dy(x) + cz(x)
#' d = roh*b*sig, c = b*sig, z = sqrt(y(x)^2+1)
#' Constrains:
#'     0 <= c <= 4sig
#'     |d| <= c
#'     |d| <= 4sig -c
#'     0 <= a <= max_i(w_i)
#' @param y (k-m)/sigma
#' @param w  total iv, iv*sqrt(tau)
#' @param sig sigma
#' @param H fractional exponent
#' @return a,b,rho parameters
#' @export
fit_inner <- function(y, w, sig, H){
  
  # constrain optimization for (a, d, c, H)
  fitfun <- function(x) {
    yhat = (x[1] + x[2]*y + x[3]*sqrt(y^2+1))^H
    sum(( yhat - w )^2)
  }
  
  # constrains (see above)
  ui = rbind(c(0,0,1),c(0,0,-1),
             c(0,-1,1), c(0, 1, 1),
             c(1, 0, 0), c(-1, 0, 0))
  ci = c(0, -4*sig, 0, 0, 0, -max(w))
  
  # initial guess
  x0 = c(max(w)/2, -sig, 2*sig)
  
  # constrained optimization
  op <- constrOptim(x0,fitfun, NULL, ui, ci)
  parameters <- op$par
  
  a = parameters[1]
  b = parameters[3]/sig
  rho = parameters[2]/parameters[3]
  return(c(a,b,rho))
}


#' Estimated svi parameters with a given H
#'
#' Estimate modified svi parameters with a given H.
#' @param k log moneyness
#' @param iv iv
#' @param tau tau
#' @param H exponent H for the modified SVI function
#' @export
fit_svih <- function(k, iv, tau, H) {
  # two-step optimization:
  # 1. obtain inner parametplers: a,b, rho
  # 2. estimate m and sigma
  w = iv*sqrt(tau[1])
  ofitfun <- function(x) { # m, sigma, exponent
    y = (k-x[1])/x[2] # m, sigma
    par <- fit_inner(y, w, x[2],H) # a, b, rho, H
    param = c(par[1],par[2], x[1], par[3], x[2], H) #a, b, m, rho, sigma, H
    raw <- svi_fun(k, param)
    sum((raw - w)^2)
  }
  
  # intial guess of m, sigma
  x0 = c( 0.0, 0.02)
  
  ui = rbind(c(1,0),c(-1,0), c(0,1), c(0,-1))
  ci = c( min(k), -max(k), 0.0001, -5)
  
  op <- constrOptim(x0, ofitfun, NULL, ui, ci)
  opar = op$par
  y = (k-opar[1])/opar[2]
  par <- fit_inner(y, w, opar[2],H)
  # parameters: a,b, m, rho, sigma, H
  parameters = data.frame(tau =tau[1], a=par[1],b=par[2],
                          m =opar[1], rho =par[3], sigma=opar[2],H= H)
}

#' parallel computing a surface with a given H
#' @param chain surface table contains date, tau, k, iv
#' @param H H exponent
#' 
par_fit_svi <- function(chain, H){

  # use parallel computing
  # detect number of cores
  no_cores <- detectCores() -1
  # Initiate cluster
  cl <- makeCluster(no_cores)
  clusterEvalQ(cl, {library(dplyr)})
  clusterExport(cl, c('fit_svih','fit_inner','svi_fun'))
  # fit function
  fit_curve <- function(c, h){
    param = fit_svih(c$k,c$iv, c$tau[1],h)
    return(param)
  }
  
  # split chains to each dte
  lchain <- split(chain, chain[,'tau'])

  # use paralell computing
  paras <- parLapply(cl, lchain,fit_curve, H)
  # convert list to data.frame
  surface = do.call(rbind, paras)
  # combine back date and dte
  ch <- chain %>% distinct(date, tau)
  surf <- left_join(surface, ch, by = c('tau'))
  # stop parallel
  stopCluster(cl)
  return(surf)
}
