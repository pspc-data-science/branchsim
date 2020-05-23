#' Renewal equation solution to the expectation of the Crump-Mode-Jagers (CMJ) over a random characteristic.
#'
#' This function numerically solves the renewal equation for a CMJ process using the properties of the triple
#' (lambda_chi, xi_chi, chi). See appendix C of Branching Processes in Biology by Kimmel and Axelrod for more details.
#' Z = G(t) + int_0^t Z(t - u) mu(du)
#'  
#' @param dMu The function dMu(t) that defines mu(dt) = dMu(t) dt.
#' @param G A function that that defines the expectation of the random characteristic of interest. Set to G=1 
#' by default (total number born).
#' @param T The maximum of the time interval for integration 0,T.
#' @param nstep The number of steps for the integration. 
#' 
#' @return A tibble containing the time steps with the solution to the renewal equation.
#'
#' @export
renewal_function<- function(dMu, G=1, T=100, nstep = 10000){
  
  times<- seq(0,T, length.out = nstep)
  
  g_dt<- dMu(times)*(times[2]- times[1])
  
  if(is.function(G)){
    G_t<- G(times) 
  }else{
    G_t<- rep(G, length(times))
  }
  
  Z<- rep(0, length(times))
  Z[1]<-G_t[1]
  for(k in 2:length(Z)){
    Z[k]<- G_t[k] + sum(Z[seq(k-1,1)]*g_dt[1:k-1])
  }
  
  return(tibble(time = times, solution = Z))
  
}

#' A function for the infinitesimal expectation in births: mu(dt) = dMu(t) dt
#'
#' @param t A vector of time evaluations
#' @param a The shape parameter of the gamma life time distribution. Default a =10
#' @param b The rate parameter of the gamma life time distribution. Default b = 1
#' @param lambda The arrival rate of infectious interactions. Default lambda = .11
#' @param p The parameter of the logarithmic distribution for the number of infected during an event.
#' Default p=0.5
#' 
#' @return Function evaluation at time vector t.
#'
#' @export
dMu<- function(t, a=10, b=1, lambda = .11, p = 0.5){
  
  r = -lambda/log(1-p)
  K<- r*p/(1-p)
  
  gamma_fun<- function(x){
    K*(1 - pracma::gammainc(b*x,a)[3])
  }
  
  
  result<- map_dbl(t, gamma_fun)
  result[t==0]<-K
  return(result)
  
}

#' A function for the expecation for the number of idividuals existing in the process at moment t. 
#' Used with the renewal eqaution, it will give the expectation for the instantaneous number infected 
#' as a function of time.
#'
#' @param t A vector of time evaluations
#' @param a The shape parameter of the gamma life time distribution. Default a =10
#' @param b The rate parameter of the gamma life time distribution. Default b = 1
#' 
#' @return Function evaluation at time vector t.
#'
#' @export
G<- function(t, a=10, b=1){
  
  gamma_fun<- function(x){
    (1 - pracma::gammainc(b*x,a)[3] )
  }
  
  result<- map_dbl(t, gamma_fun)
  result[t==0]<-1
  return(result)
}

#' A function for obtaining the Malthusian parameter of the CMJ process. It also yel
#' Used with the renewal eqaution, it will give the expectation for the instantaneous number infected 
#' as a function of time.
#'
#' @param a The shape parameter of the gamma life time distribution. Default a =10
#' @param b The rate parameter of the gamma life time distribution. Default b = 1
#' @param lambda The arrival rate of infectious interactions. Default lambda = .11
#' @param p The parameter of the logarithmic distribution for the number of infected during an event.
#' Default p=0.5
#' 
#' @return The malthusian parameter alpha, and the beta parameter.
#'
#' @export
get_malthusian<- function(a=10,b=1, lambda=.11, p=.5){
  
  r = -lambda/log(1-p)
  
  find_malthusian<-function(alpha){
    x<- r*p/(1-p)*(1 - (b/(alpha + b))^a) - alpha
  return(x)
  }
  
  alpha<- uniroot(find_malthusian, c(.001,5))$root
  
  beta<- 1/alpha*(1- a*r*p/(b*(1-p))*(b/(alpha +b))^(a+1) )
  
  return(c(alpha,beta))
  
}

#' A function for obtaining the extinction probability over all time for the process.
#'
#' @param a The shape parameter of the gamma life time distribution. Default a =10
#' @param b The rate parameter of the gamma life time distribution. Default b = 1
#' @param lambda The arrival rate of infectious interactions. Default lambda = .11
#' @param p The parameter of the logarithmic distribution for the number of infected during an event.
#' Default p=0.5
#' 
#' @return The extinction probability.
#'
#' @export
get_extinct_prob<- function(a=10, b=1, lambda = .11, p=.5){
  
  r = -lambda/log(1-p)
  
  extinct_prob<- function(s){
    (1-r/b*log((1-p)/(1-p*s)))^(-a) - s
  }
  
  result<- tryCatch({uniroot(extinct_prob, c(0,.999))$root},
                    error = function(e){
                      1
                    })
  
  return(result)
  
}

#' The generationg function of the branching process f(s) = sum_k p_k s^k.
#'
#' @param s The argument of the geberating function. 0<s<1.
#' @param a The shape parameter of the gamma life time distribution. Default a =10
#' @param b The rate parameter of the gamma life time distribution. Default b = 1
#' @param lambda The arrival rate of infectious interactions. Default lambda = .11
#' @param p The parameter of the logarithmic distribution for the number of infected during an event.
#' Default p=0.5
#' 
#' @return The evaluation the generating function at s.
#'
#' @export
generating_function<- function(s, a=10, b=1, lambda = .11, p=.5){
  r = -lambda/log(1-p)

  result<- (1-r/b*log((1-p)/(1-p*s)))^(-a)
  
  return(result)
}

# derivative_generating_function<- function(s, a=10, b=1, lambda = .11, p=.5){
#   r = -lambda/log(1-p)
# 
#   result<- a*(1-r/b*log((1-p)/(1-p*s)))^(-a-1)*r/b*p/(1-p*s)
# 
#   return(result)
# }
# 
# ave_comp<- function(u, A=10, B=1, Lambda = .11, P=.5){
# 
#   R = -Lambda/log(1-P)
#   z<- R*P/(1-P)*A
# 
#   if(z>=1){
#   s<- seq(u,1, length.out = 1000)
#   S<- pracma::trapz(s, generating_function(s, a = A, b = B, lambda = Lambda, p = P))*z
# 
#   y<- 1 - derivative_generating_function(u, a = A, b = B, lambda = Lambda, p = P)
# 
#   result<- 1 + (z*u^2)/((1-S)*y)}else{
#   result<- 1/(1-z)
#   }
# 
#   return(result)
# 
# 
# }
# 


