#' Plot the solution to the renewal equation solution for expectations of the Crump-Mode-Jagers (CMJ) over a random characteristic.
#'
#' This function plots the numerical soution to the renewal equation for a CMJ process using the properties of the triple
#' (lambda_chi, xi_chi, chi). See appendix C of Branching Processes in Biology by Kimmel and Axelrod for more details.
#' Z = G(t) + int_0^t Z(t - u) mu(du)
#'  
#' @param Lambda The arrival rate of infectious interactions.
#' @param mu The average number infected per infectious interaction
#' @param A The shape parameter of the communicable preiod gamma function.
#' @param B The rate parameter of the communicable preiod gamma function.
#' @param Time The end of the time interval of interest \code{[0,T]}.
#' @param type The type or random characteristic Total for total infected Infectious for currently infected
#' 
#' @return A plot showing the solution of the renewal equation with the malthusian parameter.
#'
#' @export
renewal_plot<- function(Lambda = .11, mu = 1.5, A=5.5, B=0.85, Time=100, type = "Total"){
  
  get_p<- Vectorize(find_p)
  p<- get_p(mu)
  
  malthusian_parameter<-  get_malthusian(a=A, b=B, lambda = Lambda, p=p)[1]
  
  if(malthusian_parameter == 0){
    malthusian_parameter<- as.character(malthusian_parameter)
  }else{
    malthusian_parameter<- sprintf('%.4f', malthusian_parameter)
  }
  
  g_title<- str_c("Malthusian parameter: ", malthusian_parameter)
  
  if(type == "Total"){
    dmu_fun<- dMu(A=A, B=B, Lambda = Lambda, P = p)
    eval_tibble<-renewal_function(dmu_fun, G=1, Time_limit=Time, nstep = 10000)
    G<- ggplot(eval_tibble, aes(time, solution)) + geom_line(col = "blue") + 
      theme_bw()+ theme(text = element_text(size=15)) + ylab("Expected total infected")+ xlab("Time (days)")+
      scale_x_continuous(breaks = pretty_breaks(10)) +
      scale_y_continuous(breaks = pretty_breaks(10)) +
      ggtitle(g_title)
  }else{
    dmu_fun<- dMu(A=A, B=B, Lambda = Lambda, P = p)
    eval_tibble<-renewal_function(dmu_fun, G=G(A=A,B=B), Time_limit=Time, nstep = 10000)
    G<- ggplot(eval_tibble, aes(time, solution)) + geom_line(col = "blue") + 
      theme_bw()+ theme(text = element_text(size=15)) + ylab("Expected infectious")+ xlab("Time (days)")+
      scale_x_continuous(breaks = pretty_breaks(10)) +
      scale_y_continuous(breaks = pretty_breaks(10))+
      ggtitle(g_title)
    
  }
  
  print(G)
  
}

#' @export
asymptotic_plot<- function(lambda_limits, mu_limits, a=5.5, b=.85){
  
  mu<- seq(mu_limits[1],mu_limits[2], length.out = 40)
  lambda<- seq(lambda_limits[1],lambda_limits[2], length.out = 40)
  
  get_p<- Vectorize(find_p)
  p<- get_p(mu)
  
  data_lp<- expand.grid(lambda = lambda, p=p) %>% as_tibble()
  
  Z<- map2(data_lp$lambda, data_lp$p, get_malthusian, a=a, b=b) %>% reduce(rbind) %>% .[,1]
  data_lp$malthusian<- Z
  data_lp<- data_lp %>% mutate(mu = -p/((1-p)*log(1-p)))
  
  X<- tibble(x = seq(lambda_limits[1],lambda_limits[2], length.out = 400), y = 1/(x*a/b))
  X<- X %>% filter(y<= max(data_lp$mu)) %>% filter(y>= min(data_lp$mu))
  
  #my_breaks<- min
  ggplot(data_lp, aes(lambda,mu, fill=malthusian)) + geom_raster(interpolate = T) + 
    scale_fill_viridis_c(breaks = pretty_breaks(10)) + 
    theme_bw()+ theme(text = element_text(size=15))+
    xlab("Average number of infectious events per day") +
    ylab("Average number infected per infectious event") + labs(fill="malthusian parameter")+
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) + 
    geom_line(inherit.aes = FALSE,data = X, aes(x,y))+
    scale_x_continuous(expand=c(0,0),breaks = pretty_breaks(10)) + 
    scale_y_continuous(expand=c(0,0), breaks = pretty_breaks(10))
}

#' @export
extinct_plot<- function(lambda_limits, mu_limits, a=5.5, b=.85){
  
  mu<- seq(mu_limits[1],mu_limits[2], length.out = 40)
  lambda<- seq(lambda_limits[1],lambda_limits[2], length.out = 40)
  
  
  get_p<- Vectorize(find_p)
  p<- get_p(mu)
  
  data_lp<- expand.grid(lambda = lambda, p=p) %>% as_tibble()
  
  Z<- map2_dbl(data_lp$lambda, data_lp$p, get_extinct_prob, a=a, b=b)
  data_lp$extinct<- Z
  data_lp<- data_lp %>% mutate(mu = -p/((1-p)*log(1-p)))
  
  X<- tibble(x = seq(lambda_limits[1],lambda_limits[2], length.out = 400), y = 1/(x*a/b))
  X<- X %>% filter(y<= max(data_lp$mu)) %>% filter(y>= min(data_lp$mu))
  
  ggplot(data_lp, aes(lambda,mu, fill=extinct)) + geom_raster(interpolate = T) + 
    scale_fill_viridis_c(breaks = pretty_breaks(10)) + 
    theme_bw()+ theme(text = element_text(size=15))+
    #scale_x_continuous(expand=c(0,0), breaks = pretty_breaks(10)) + 
    #scale_y_continuous(expand=c(0,0), breaks = pretty_breaks(10)) + 
    geom_line(inherit.aes = FALSE,data = X, aes(x,y))+
    scale_x_continuous(expand=c(0,0),breaks = pretty_breaks(10)) + 
    scale_y_continuous(expand=c(0,0), breaks = pretty_breaks(10))+
    xlab("Average number of infectious events per day") +
    ylab("Average number infected per infectious event") + labs(fill="extinction probability")+
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20))
}

#' @export
component_size_plot<- function(lambda_limits, mu_limits, a=5.5, b=.85){
  
  mu<- seq(mu_limits[1],mu_limits[2], length.out = 50)
  lambda<- seq(lambda_limits[1],lambda_limits[2], length.out = 50)
  
  
  get_p<- Vectorize(find_p)
  p<- get_p(mu)
  
  data_lp<- expand.grid(lambda = lambda, p=p) %>% as_tibble()
  
  Z<- map2_dbl(data_lp$lambda, data_lp$p, get_extinct_prob, a=a, b=b)
  data_lp$extinct<- Z
  data_lp<- data_lp %>% mutate(mu = -p/((1-p)*log(1-p)))
  
  args_list<- list(u = as.list(Z), A = as.list(rep(a, nrow(data_lp))), 
                   B = as.list(rep(b, nrow(data_lp))), Lambda = as.list(data_lp$lambda), P=as.list(data_lp$p))
  
  V<- pmap_dbl(args_list, average_component_size)
  
  data_lp$comp_size<- V
  #data_lp<- data_lp %>% filter(extinct ==1)
  
  X<- tibble(x = seq(lambda_limits[1],lambda_limits[2], length.out = 400), y = 1/(x*a/b))
  X<- X %>% filter(y<= max(data_lp$mu)) %>% filter(y>= min(data_lp$mu))
  my_breaks<- c(3, 10, 30, 100)
  ggplot(data_lp, aes(lambda,mu, fill=comp_size)) + geom_raster(interpolate = TRUE) + 
    scale_fill_viridis_c(breaks = my_breaks, trans = "log10", na.value = "white") + 
    theme_bw()+ theme(text = element_text(size=15))+
    xlab("Average number of infectious events per day") +
    ylab("Average number infected per infectious event") + labs(fill="extinct mean\ncomponent size")+
    guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) + 
    geom_line(inherit.aes = FALSE,data = X, aes(x,y, col = ""), size=3)+
    scale_x_continuous(expand=c(0,0),breaks = pretty_breaks(10)) + 
    scale_y_continuous(expand=c(0,0), breaks = pretty_breaks(10))+labs(colour="phase transition")
  
  
}

#' @export
plot_single_mother_dist<- function(a=5.5,b=.85,lambda=.11, p=.5){
  
  X<- single_mother_distribution(a=a,b=b,lambda=lambda, p=p) %>% filter(prob>3*10^-4)
  
  
  ggplot(X, aes(count, prob)) + geom_segment(aes(x= count, xend = count, y = 0, yend=prob), size =1.5) + 
    geom_point(size = 5, color = "red", fill=alpha("orange", 0.3), alpha=0.7, shape=21, stroke=2) +
    theme_bw() + theme(text = element_text(size=15))+
    xlab("Number infected during communicable period") + ylab("Probability")+
    scale_x_continuous(breaks = pretty_breaks(10)) + scale_y_continuous(breaks = pretty_breaks(10))
  
}

#' @export
plot_gamma_dist<- function(a=5.5, b= .85){

  x<- seq(0, 3*a, length.out = 1000)
  y<- dgamma(x, shape = a, rate = b)
  X<- tibble(x=x, y=y)
  
  
  ggplot(X, aes(x,y)) + geom_line(col = "blue") + theme_bw() + 
    geom_ribbon(aes(ymax=y),ymin=0,fill="blue",colour=NA,alpha=0.3)+
    theme(text = element_text(size=15)) +
    xlab("Communicable period") + ylab("Density") + scale_x_continuous(expand = c(0,0), breaks = pretty_breaks(10))+
    scale_y_continuous(expand = c(0,0),breaks = pretty_breaks(10)) + geom_vline(xintercept = a/b, col = "red")
  
}



