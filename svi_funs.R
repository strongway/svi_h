# functions of stochastic volatility inspired paramater estimations
# based on Phillip Rindler's matlab code
# Strongway 23.09.2016

#require('Hmisc', quiet = TRUE) # approxExtrap
#require('RQuantLib', quiet = TRUE)

#require('GA', quiet = TRUE)

#' SVI Raw formulation
#'
#' @import Hmisc
#' @import RQuantLib
#' the raw formulation from Gatheral's book equation 3.20
#' @param k log-moneyness at which to evaluate the total implied variance
#' @param (5x1) parameters:
#'    a = scalar = level of variance
#'    b = scalar = slope of wings
#'    m = scalar = translates smile to right
#'    rho = scalar = counter-clockwise rotation of smile
#'  sigma = scalar = reduces ATM curvature of the smile
#' @return A list contain
#'    totalvariance = (Nx1) = estimated total variance at k
#'    impliedvolatility = (Nx1) = estimated implied volatility at (k,tau)
#' @export
svi_raw <- function(k,param,tau=1){
  if (length(param) != 5)
    stop('There have to be five parameters: a, b, m, rho, sigma')

  a = param[1];
  b = param[2];
  m = param[3];
  rho = param[4];
  sigma = param[5];

  # valid parameters
  if (b <0) stop('b has to be non-negative')
  if (abs(rho)>1) stop('|rho| has to be smaller than 1')
  if (sigma < 0)  stop( 'sigma has to be positive')
  if(a + b*sigma*sqrt(1-rho^2)<0) stop( 'a + b sigma (1-rho^2)^0.5 has to be non-negative');

  # calculate total variance
  totalvariance = a + b*(rho*(k-m) + sqrt((k-m)^2 + sigma^2))
  impliedvolatility = sqrt(totalvariance/tau)
  result <- list(totalvariance=totalvariance, impliedvolatility=impliedvolatility)
  result
}

#' SVI convert parameters
#'
#' svi_convertparameters converts the parameter set of one type of SVI
#' formulation to another. The parameterizations are assumed to be:
#' * raw =(a,b,m,rho, sigma)
#' * natural = (delta, mu, rho, omega, zeta)
#' * jumpwing = (v, psi, p, c, vt)
#'
#'
#' @param param_old (5x1) original parameters
#' @param from string = formulation of original parameters (raw, natural, jumpwing)
#' @param to  string = formulation of new parameters (raw, natural, jumpwings)
#'
#' @return
#' param_new = (5x1) = new parameters
#' @export
svi_convertparameters <- function(param_old, from, to, tau=1) {
  if (length(param_old) != 5)
    stop('There have to be five parameters.')

  if (!(from %in% c('raw','natural','jumpwing')))
      stop('from has to be one of: raw, natural, jumpwing')

  if (!(to %in% c('raw','natural','jumpwing')))
    stop('from has to be one of: raw, natural, jumpwing')

  switch(from,
       'raw' = {
         a = param_old[1]
         b = param_old[2]
         m = param_old[3]
         rho = param_old[4]
         sigma = param_old[5]
         switch(to,
           'raw' = {
             param_new = param_old },
           'natural' = {
               omega = 2*b*sigma/sqrt(1-rho^2)
               delta = a - omega/2*(1-rho^2)
               mu = m + rho*sigma/sqrt(1-rho^2)
               zeta = sqrt(1-rho^2)/sigma
               param_new = c(delta, mu, rho, omega, zeta)  },
           'jumpwing' ={
               w = a + b*(-rho*m + sqrt(m^2+sigma^2))
               v = w/tau
               psi = 1/sqrt(w) * b/2 * (-m/sqrt(m^2+sigma^2) + rho)
               p = 1/sqrt(w) * b * (1-rho)
               c = 1/sqrt(w) * b * (1+rho)
               vt = 1/tau * (a + b*sigma*sqrt(1-rho^2))
               param_new = c(v, psi, p, c, vt)
           }
         )
       }, # end of 'raw'
      'natural'={
        switch(to,
          'raw' = {
              delta = param_old[1]
              mu = param_old[2]
              rho = param_old[3]
              omega = param_old[4]
              zeta = param_old[5]

              a = delta + omega/2*(1-rho^2)
              b = omega*zeta/2
              m = mu - rho/zeta
              sigma = sqrt(1-rho^2)/zeta
              param_new = c(a, b, m, rho, sigma) },
          'natural' = {
              param_new = param_old },
          'jumpwing' = {
              param_temp = svi_convertparameters( param_old, 'natural', 'raw', tau)
              param_new = svi_convertparameters( param_temp, 'raw', 'jumpwing', tau)
          }
        )
      }, # end of 'natural'
      'jumpwing'={
        switch(to,
          'raw' = {
              v = param_old[1]
              psi = param_old[2]
              p = param_old[3]
              c = param_old[4]
              vt = param_old[5]
              w = v*tau

              b = sqrt(w)/2 * (c+p)
              rho = 1-p*sqrt(w)/b
              beta = rho - 2*psi*sqrt(w)/b
              alpha = sign(beta)*sqrt(1/beta^2-1)
              m = (v-vt)*tau/(b*(-rho+sign(alpha)*sqrt(1+alpha^2)-alpha*sqrt(1-rho^2)))
              if (is.nan(m) | abs(m)< 1e-5)
                sigma = (vt*tau - w)/b / (sqrt(1-rho^2)-1)
              else
                sigma = alpha*m

              a = vt*tau - b*sigma*sqrt(1-rho^2);

              if (sigma < 0)
                sigma = 0
              param_new = c(a, b, m, rho, sigma)
          },
        'natural' = {
            param_temp = svi_convertparameters( param_old, 'jumpwing', 'raw', tau)
            param_new = svi_convertparameters( param_temp, 'raw', 'natural', tau)
        },
        'jumpwing' = {
          param_new = param_old
          }
        )
      }
  ) # end of switch 'from'

  param_new
}

#' SVI Natural formulation
#'
#'SVI - Stochastic Volatility Inspired parameterization of the implied
#'volatility smile. This function implements the natural formulation.
#'
#' @param k (Nx1) log-moneyness at which to evaluate the total implied
#' variance
#' @param param = (5x1) = parameters:
#'   * delta = scalar = level of variance
#'   * mu = scalar = slope of wings
#'   * rho = scalar = translates smile to right,
#'   * omega = scalar = counter-clockwise rotation of smile
#'   * zeta = scalar = reduces ATM curvature of the smile
#' @return
#' * totalvariance = (Nx1) = estimated total variance at k
#' * impliedvolatility = (Nx1) = estimated implied volatility at (k,tau)
#' @export
svi_natural <- function(k, param, tau=1) {

    # check that input is consistent
    if (length(param) != 5)
      stop('There have to be five parameters.')

      delta = param[1]
      mu = param[2]
      rho = param[3]
      omega = param[4]
      zeta = param[5]

      # make sure that parameter restrictions are satisfied
      if(omega < 0) stop( 'omega has to be non-negative')
      if(abs(rho) >= 1) stop( '|rho| has to be smaller than 1')
      if(zeta <= 0) stop( 'zeta has to be positive')
      if(delta + omega*(1-rho^2) < 0) stop( 'delta + omega (1-rho^2) has to be non-negative')

      # calculate total variance
      totalvariance = delta + omega/2*(1+ zeta*rho*(k-mu) + sqrt((zeta*(k-mu)+rho)^2+(1-rho^2)))
      impliedvolatility = sqrt(totalvariance/tau);
      result <- list(totalvariance=totalvariance, impliedvolatility=impliedvolatility)
      result
}

#' SVI Jump wing formulation
#'
#' SVI - Stochastic Volatility Inspired parameterization of the implied
#' volatility smile. This function implements the jump-wings formulation.
#'
#' @param k = (Nx1) = log-moneyness at which to evaluate the total implied
#' variance
#' @param param = (5x1) = parameters:
#'   * v = scalar = ATM variance
#'   * psi = scalar = ATM skew
#'   * p = scalar = slope of left/put wing
#'   * c = scalar = slope of right/call wing
#'   * vt = scalar = minimum implied variance
#' @return
#' * totalvariance = (Nx1) = estimated total variance at k
#' * impliedvolatility = (Nx1) = estimated implied volatility at (k,tau)
#' @export
svi_jumpwing <- function(k, param, tau=1){

    # check that input is consistent
      if (length(param) != 5)
        stop('There have to be five parameters.')
    # convert parameters to raw formulation
      param_raw = svi_convertparameters( param, 'jumpwing', 'raw', tau)

      # calculate total variance
      result = svi_raw(k, param_raw, tau)
      result

}

#'SVI Surface
#'
#'svi_surface calcualtes the surface SVI free of statis arbitrage.
#'
#'
#' @param k = (Kx1) = moneyness at which to evaluate the surface
#' @param theta = (Tx1) = ATM variance time at which to evaluate the surface.
#' If theta and k have the same dimensions, then the surface is evaluated
#' for each given (k,theta) pair. If the dimensions are different, the
#' function evaluates the surface on the grid given by k and theta
#' @param param = (Xx1) = parameters of the surface
#'
#' @return
#' * totalvariance = (KxT) or (Kx1) = estimated total variance
#' * impliedvolatility = (KxT) or (Kx1) = estimated implied volatility
#' @export
svi_surface <- function(k, theta, rho, phifun, phi_param, tau=1) {

    if (length(rho)>1) stop('rho has to be a scalar')

    #evaluate phi
    phi = switch (phifun,
        heston_like =  heston_like(theta, phi_param),
        power_law = power_law(theta, phi_param) )

    if (length(k) ==length(theta))
      totalvariance = theta/2 * (1 + rho * phi*k + sqrt((phi*k + rho)^2 + (1-rho^2)))
    else {
      T = length(theta)
      K = length(k)
      totalvariance = matrix(0,K,T)
    for (t in c(1:T)) {
        totalvariance[,t] = theta[t]/2* (1 + rho * phi[t]*k +
                                        sqrt((phi[t]*k + rho)^2 + (1-rho^2)))
      }
    }
    if (sum(abs(dim(totalvariance)-dim(tau)))==0)
      impliedvolatility = sqrt(totalvariance/tau)
    else {
      # expand tau
      tau_expanded = matrix(0, dim(totalvariance)[1],dim(totalvariance)[2])
      maturities = unique(tau)
      T = length(maturities)
      for (t in c( 1:T))
        tau_expanded[,t] = maturities[t]

      impliedvolatility = sqrt(totalvariance/tau_expanded)
    }
    result <- list(totalvariance=totalvariance, impliedvolatility=impliedvolatility)
    result

}


#' Heston like parameterization
#'
#' input: theta, param
#' Return: heston parameter
heston_like <- function(theta, param){
    # Heston-like parameterization
    lambda = param[1]
    value = 1./(lambda*theta)* (1 - (1 - exp(-lambda*theta))/(lambda*theta))
    value
}


#' Power-law parameterization
#'
#' Input: theta, param
#' Return: powerlaw parameter
power_law <- function(theta, param){
    # Power-law parameterization
    eta = param[1]
    gamma = param[2];

    value = eta/(theta^(gamma) * (1+theta)^(1-gamma))
    value
}

#' Stochastic Volatility Inspired parameterization
#'
#' SVI - Stochastic Volatility Inspired parameterization of the implied
#' volatility smile. This function implements combination of three types formulation.
#'
#' @param k = (Kx1) = moneyness at which to evaluate the surface
#' @param param = parameters of the surface
#' @param type = 'raw', 'natural','jumpwing'
#' @param tau expiration
#' @export
svi <- function(k, param, type, tau = 1) {
  # combination of three types of svi methods
  result <- switch(type,
          raw = svi_raw(k,param, tau),
          natural = svi_natural(k,param,tau),
          jumpwing = svi_jumpwing(k,param,tau))
  result
}

#' SVI parameter bounds
#'
#' svi_parameter_bounds returns the parameter bounds of an SVI parameterization. The parameters are
#' assumed to be in the following order:
#' raw = [a b m rho sigma]
#' natural = [delta mu rho omega zeta]
#' jumpwing = [v psi p c vt]
svi_parameter_bounds <- function(parameterization){

  large = 1e5

  bounds <- switch(parameterization,
            raw = list(lb = c(-large, 0,-large,-1, 0),
                       ub =c(large,large,large,1, large)),
            natural = list( lb = c(-large, -large, -1, 0, 0),
                            ub = c(large,large,  1, large, large)),
            jumpwing = list(lb = c(0,-large, 0,0,0),
                            ub = c(large, large, large, large, large)))
  bounds
}

#' A simplified version of svi_interpolation
#'
#' Output:
#' * implied_volatility (Kx1) = implied volatilities corresponding to total_implied_variance
#' para is the output of the fitted surface function
#' @export
svi_iv_interp <- function(log_moneyness, tau_interp=1, para) {
  parameters = t(as.matrix(para[,c('X1','X2','X3','X4','X5')]))
  theta = as.numeric(para$theta)
  maturities = as.numeric(para$maturity)
  #ensure scalar input
  if (!is.null(dim(tau_interp))) stop('tau_interp has to be scalar')

  if (tau_interp < min(maturities) | tau_interp > max(maturities)) stop('maturity exceeds the range.')

  # estimate theta for interpolated maturity
  #  theta_interp = interp1(maturities, theta, tau_interp, 'linear')
  theta_interp = approxExtrap(maturities,theta, xout = tau_interp)$y
  if (tau_interp %in% maturities){
    iv = svi_jumpwing(log_moneyness, parameters[,maturities==tau_interp], tau_interp)
    total_implied_variance = iv$totalvariance
    implied_volatility = sqrt(total_implied_variance/tau_interp)
  } else {
      # interpolation
      idx = which.min(abs(tau_interp-maturities))
      # if closest maturity is smaller than tau_interp, make idx one unit larger --> idx is index of
      # smallest maturity larger than tau_interp
      if (maturities[idx] < tau_interp)
        idx = idx+1

      alpha_t = (sqrt(theta[idx]) - sqrt(theta_interp))/(sqrt(theta[idx]) - sqrt(theta[idx-1]))
      param_interp = alpha_t*parameters[,idx-1] + (1-alpha_t)*parameters[,idx]
      iv = svi_jumpwing(log_moneyness, param_interp, tau_interp)
      total_implied_variance = iv$totalvariance
      implied_volatility = sqrt(total_implied_variance/tau_interp)
  }
  implied_volatility
}


#' Butterfly arbitrage free constrain
#'
constrain1 <- function(x,phifun = 'power_law') {
# parameters = [rho, lambda] subject to
#  constraints: in heston_like: (1 + |rho|) <= 4 lambda, in power-law: eta(1+|rho|) <= 2
  switch(phifun,
         heston_like = 1 + abs(x[1]) - 4*x[2],
         power_law =  x[2] * (1 + abs(x[1])) - 2)
}

#' fit function for GA optimization
fitfunction1 <- function(x, phifun='power_law', log_moneyness, theta, total_implied_variance) {
  # extract parameters from x
  rho = x[1]
  phi_param = switch(phifun,
                     heston_like = x[2],
                     power_law = c(x[2],0.5))

  iv = svi_surface(log_moneyness, theta, rho,phifun, phi_param)
  model_total_implied_variance = iv$totalvariance
  # inverse the rmse for maximization
  value = 1/norm(as.matrix(total_implied_variance - model_total_implied_variance),
                         type = '2')
  # add penalty for the parameter constrains
  pen <- sqrt(.Machine$double.xmax)/4  # penalty term
  penalty1 <- max(constrain1(x),0)*pen
  penalty2 <- (x[2]<0)*pen # x2 must be positive
  penalty3 <- (x[1]<=-1)*pen # rho must be >-1
  value - penalty1 -penalty2 - penalty3        # fitness function value

}

#' fit function for initial optimization
fitfunction3 <- function(x, phifun='power_law', log_moneyness, theta, total_implied_variance) {
  # extract parameters from x
  rho = x[1]
  phi_param = switch(phifun,
                     heston_like = x[2],
                     power_law = c(x[2],0.5))

  iv = svi_surface(log_moneyness, theta, rho,phifun, phi_param)
  model_total_implied_variance = iv$totalvariance
  #
  value = norm(as.matrix(sqrt(total_implied_variance) - sqrt(model_total_implied_variance)),
                 type = '2')
  # add penalty for the parameter constrains
  pen <- sqrt(.Machine$double.xmax)/4  # penalty term
  penalty1 <- max(constrain1(x),0)*pen
#  penalty2 <- (x[2]<0)*pen # x2 must be positive
#  penalty3 <- (x[1]<=-1)*pen # rho must be >-1
  value +penalty1        # fitness function value

}

#' Fit SSVI
#'
#' Inputs:
#' phifun = default 'power_law'
#' log_moneyness
#' theta
#' total implied variance
fit_ssvi <- function(phifun='power_law', log_moneyness, theta_expanded, total_implied_variance){

# method 1: genetic algorithm (see volModel)
# method 2: constrain optimization
  fitfun <- function(x) fitfunction3(x, phifun, log_moneyness , theta_expanded, total_implied_variance)

  ui = rbind(c(1,0),c(-1,0), c(0,1))
  ci = c(-1, -1, 0)

  x0 = c(-0.3, 1.5)

  # constrained optimization
  op <- constrOptim(x0,fitfun, NULL, ui, ci)
  parameters <- op$par

  if (phifun == 'power_law')
    parameters = c(parameters,0.5)

  parameters

}

#' fit surface with all vol
#'
#' Inputs:
#'    implied volatility
#'    maturity
#'    log moneyness
#'    phifun
fit_raw_surface <- function(implied_volatility, maturity, log_moneyness,
                            phifun='power_law'){
  # step one: estimate total implied variance
  total_implied_variance = implied_volatility^2* maturity;

  # step two: use linear interpolation for ATM total implied variance
  maturities = sort(unique(maturity))
  T = length(maturities)
  theta = rep(0,T)
  theta_expanded = rep(0, length(maturity))
  # find atm total-implied-variance (i.e. at log_moneyness=0)
  for (t in c(1:T)){
    pos = which(maturity==maturities[t])
    tiv_t = total_implied_variance[pos]
    if (0 %in% log_moneyness[pos]){
      theta[t] = tiv_t[log_moneyness[pos]==0]
    } else{
#      theta[t] = interp1(log_moneyness[pos], tiv_t, 0, 'linear')
      theta[t] = approxExtrap(log_moneyness[pos], tiv_t, xout =0)$y
    }
    theta_expanded[pos] = theta[t]
  }


  # step three: fit SVI surface by estimating parameters = [rho, lambda] subject to parameter bounds:
  # -1 < rho < 1, 0 < lambda
  # and constraints: in heston_like: (1 + |rho|) <= 4 lambda, in power-law: eta(1+|rho|) <= 2
  parameters = fit_ssvi(phifun, log_moneyness, theta_expanded, total_implied_variance)
  list(parameters = parameters, theta = theta, theta_expanded = theta_expanded)
}

#' Fit SVI surface
#'
#' fit_svi_surface calibrates the SVI surface to market data. First, the entire Surface SVI is fitted
#' to all log-moneyness-theta observations. Second, each slice is fitted again using the SSVI fit as
#' initial guess.
#' Input:
#' * impliedvolatility, maturity, moneyness, phi
#' Output:
#' * parameters = (2x1) = parameters of SSVI = [rho, lambda]
#' * theta      = (Tx1) = ATM total implied variance
#' * maturities = (Tx1) = corresponding time to maturity
#' @export
fit_svi_surface <- function(implied_volatility, maturity, log_moneyness,
                            phifun='power_law') {

  # step one: estimate total implied variance
  total_implied_variance = implied_volatility^2* maturity;
  maturities = sort(unique(maturity))
  T = length(maturities)
  # step two: fit entire surface using GA algo
  fit1 <- fit_raw_surface(implied_volatility, maturity, log_moneyness, phifun)

  parameters = fit1$parameters
  theta = fit1$theta

  rho = parameters[1]
  phi_param = parameters[-1]

#  print(parameters) # debug

  # step four: transform SSVI parameters to SVI-JW parameters
  phi = switch(phifun,
               heston_like = heston_like(theta, phi_param),
              power_law = power_law(theta, phi_param))

  v = theta/maturities
  psi = 0.5*rho*sqrt(theta)*phi
  p = 0.5*sqrt(theta)*phi*(1-rho)
  c = p + 2*psi
  vt = v*(4*p*c)/((p+c)^2)

  # step five: iterate through each maturity and fit c and vt for best fit
  parameters = matrix(0,5,T)

  for (t in seq(T,1, by=-1)) {
      pos <-  maturity==maturities[t]
      log_moneyness_t = log_moneyness[pos]
      total_implied_variance_t = total_implied_variance[pos]

      if (t==T){
        param_before = c(v[t-1], psi[t-1], p[t-1], c[t-1], vt[t-1])
        iv = svi_jumpwing(log_moneyness_t, param_before, maturities[t-1])
        slice_before = iv$totalvariance
        slice_after = NULL
      }
      else if (t==1) {
        slice_before = NULL
        param_after = c(v[t+1], psi[t+1], p[t+1], c[t+1], vt[t+1])
        iv = svi_jumpwing(log_moneyness_t, param_after, maturities[t+1])
        slice_after = iv$totalvariance
      } else {
        param_before = c(v[t-1], psi[t-1], p[t-1], c[t-1], vt[t-1])
        iv = svi_jumpwing(log_moneyness_t, param_before, maturities[t-1])
        slice_before = iv$totalvariance
        param_after = c(v[t+1], psi[t+1], p[t+1], c[t+1], vt[t+1])
        iv2 = svi_jumpwing(log_moneyness_t, param_after, maturities[t+1])
        slice_after = iv2$totalvariance
      }
      param0 = c(v[t], psi[t], p[t], c[t], vt[t])

 #     print(param0) # debug

      parameters[1:3,t] = fit_svi(param0, log_moneyness_t, total_implied_variance_t,
                                  slice_before, slice_after,maturities[t])
      parameters[4,t] = parameters[3,t] + 2*parameters[2,t]
      parameters[5,t] = parameters[1,t]*4*parameters[3,t]*parameters[4,t]/(parameters[3,t]+parameters[4,t])^2
  }
  list(parameters = parameters, theta = theta, maturities=maturities)
}

#' FIT SVI
#'
#' fit_svi fits a volatility slice to observed implied volatilities.
#' Input:
#'   x0, k, total implied variance, slice before, slice after, tau
#' Output:
#'   estimated surface parameters
#'
#' @export
fit_svi <- function(x0, k, total_implied_variance, slice_before, slice_after, tau){
# using constrOptim instead of GA, given x0 is available
  # get variable bounds
  large = 1e5
  small = 1e-6
  # linear inequality: -p <= 2*psi, i.e., 2*x[2]

  ui = rbind(c(1,0,0),c(-1,0,0), c(0,1,0),c(0,-1,0),c(0,0,1),c(0,0,-1),c(0,2,1))
  ci = c(small, -large, -large, -large,small,-large, small)
  #ui = rbind(c(1,0,0),c(-1,0,0), c(0,1,0),c(0,-1,0),c(0,0,1),c(0,0,-1))
  #ci = c(small, -large, -large, -large,small,-large)

  # only optimize first three variables, final two are set by no-arbitrage condition
  x0 = x0[1:3]

  targetfun <- function(x) fitfunction2(x, k, total_implied_variance, slice_before, slice_after, tau)

  # constrained optimization
  op <- constrOptim(x0,targetfun, NULL, ui, ci)
  parameters <- op$par
}

#' The objective function of the minimization of fit_svi.
#'
#' The objective is the 2-norm
#' of the error in total implied variance. If specified, the function tests whether the current total
#' variance slice lies between the prior and later slice. If the arbitrage bound is violated, the
#' function is penalized with a very large value.
#' Input:
#' * x = (5x1) = parameters of SVI model
#' * k = (Nx1) = log moneyness at which total implied variance is calculated
#' * total_implied_variance = (Nx1) = market observed total implied variance
#' * slice_before = (Nx1) = total implied variance with smaller maturity than current slice
#' * slice_after = (Nx1) = total implied variance with larter maturity than current slice
fitfunction2 <- function(x, k, total_implied_variance, slice_before, slice_after, tau){

  # check that input is consistent
  # assert(isequal(size(k), size(total_implied_variance)), 'moneyness and total implied variance need to have the same size');
  # if ~isempty(slice_before)
  # assert(isequal(size(k), size(slice_before)), 'moneyness and slice_before need to have the same size');
  # end
  # if ~isempty(slice_after)
  # assert(isequal(size(k), size(slice_after)), 'moneyness and slice_after need to have the same size');
  # end

  # expand x
  x[4] = x[3] + 2*x[2]
  x[5] = x[1]*4*x[3]*x[4]/(x[3]+x[4])^2

  # calculate model total implied variance given current parameter estimate
  iv = svi_jumpwing(k, x, tau)
  model_total_implied_variance = iv$totalvariance
  # calculate value of objective function
  # add weighted normalization (strongway) - no improvement
#  weights = (max(total_implied_variance)-total_implied_variance)/sd(tot_var)
  value = norm(as.matrix((sqrt(total_implied_variance)-sqrt(model_total_implied_variance)) ),
              type = '2')

  # if the current model total implied variance crosses the earlier or later slice, set value to 1e6
  if (!is.null(slice_before) && any(model_total_implied_variance<slice_before))
   value = value+1e6

  if (!is.null(slice_after) && any(model_total_implied_variance>slice_after))
   value = value+1e6

  if (is.nan(value))    value = +1e6

  value

}



#' interpolation of svi parameters
#'
#' @param para surface parameters
#' @param tau_interp to-be-interpolated tau
#' @return
#' interpoplated parameter at tau_interp
#' Note: Interpolation works only within the maturity range.
#' In this version, interpolation also for the linear part
#' @return values contains 5 of svi parameters, 3 of linear model
#' @export
svi.interparam <- function(para, tau_interp) {
  if ("slope" %in% colnames(para))
    # hybrid model
    parameters = t(as.matrix(para[,c('X1','X2','X3','X4','X5','intercept', 'slope', 'square','theta')]))
  else # only svi parameters
    parameters = t(as.matrix(para[,c('X1','X2','X3','X4','X5')]))

  theta = as.numeric(para$theta)
  maturities = as.numeric(para$maturity)
  #ensure scalar input
  if (!is.null(dim(tau_interp))) stop('tau_interp has to be scalar')
  if (tau_interp < min(maturities) | tau_interp > max(maturities)) stop('maturity exceeds the range.')

  if (tau_interp %in% maturities){
    param_interp = parameters[,maturities==tau_interp]
  } else {
    # estimate theta for interpolated maturity
    theta_interp = approxExtrap(maturities,theta, xout = tau_interp)$y
    # interpolation
    idx = which.min(abs(tau_interp-maturities))
    # if closest maturity is smaller than tau_interp, make idx one unit larger --> idx is index of
    # smallest maturity larger than tau_interp
    if (maturities[idx] < tau_interp)
      idx = idx+1

    alpha_t = (sqrt(theta[idx]) - sqrt(theta_interp))/(sqrt(theta[idx]) - sqrt(theta[idx-1]))
    param_interp = alpha_t*parameters[,idx-1] + (1-alpha_t)*parameters[,idx]
  }
  param_interp
}

#' Fit Chain SVI parameters
#'
#' Fit iv surface using svi approach with a given day chain
#'
#' @param chain an iv option chain containing iv, maturity, logmoneyness
#' Note: chain should be from one single day
#' @return surface parameter
#'
#' @export
svi.fitparam <- function(chain){
  # fit surface, and return data.frame parameters
  chain.Para = fit_svi_surface(chain$iv, chain$maturity, chain$logmoneyness)

  # convert parameters matrix to data.frame
  para = data.frame(list(theta = chain.Para$theta,
                         maturity = chain.Para$maturities, t(chain.Para$parameters)))
  #  para$symbol = chain$symbol[1]
  #  para$date = chain$date[1]
  ch <- chain %>% distinct(dte, maturity, mF, R, U, spot)

  left_join(para, ch, by = c('maturity'))
}

#' Fit Chain surface with hybrid svi model
#'
#' Fit iv surface using svi approach with a given day chain,
#' left side with a quadra-linear model, right side: svi
#'
#' @param chain an iv option chain containing iv, maturity, logmoneyness
#' Note: chain should be from one single day
#' @return surface parameter
#'
#' @export
svi.hybridparam <- function(chain, boundary = -0.02){
  # fit surface, and return data.frame parameters
  # part 1: surface with logmoneyness > - 0.05
  chain1 <- chain %>% dplyr::filter(logmoneyness >= boundary)
  chain.Para = fit_svi_surface(chain1$iv, chain1$maturity, chain1$logmoneyness)

  # part 2: surface with logmoneyness < -0.03 (linear part)
  chain2 <- chain %>% dplyr::filter(logmoneyness <= boundary/2)

  lm.Para <- chain2 %>% group_by(maturity) %>%
    do(tidy(lm(iv*sqrt(maturity) ~ poly(logmoneyness,2, raw = TRUE) , data = .))) %>%
    select(maturity, term, estimate) %>% spread(term, estimate)
  names(lm.Para) <- c('maturity', 'intercept', 'slope', 'square')

  # convert parameters matrix to data.frame
  para = data.frame(list(theta = chain.Para$theta,
                         maturity = chain.Para$maturities, t(chain.Para$parameters)))
  #  para$symbol = chain$symbol[1]
  #  para$date = chain$date[1]
  para <- left_join(para, lm.Para, by = c('maturity'))
  para$boundary <-  boundary

  ch <- chain %>% distinct(dte, maturity, mF, R, U, spot)

  left_join(para, ch, by = c('maturity'))
}

