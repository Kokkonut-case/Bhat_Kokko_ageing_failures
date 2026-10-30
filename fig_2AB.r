############################
# Make plots comparing Gompertz-Makeham and logistic failure/mortality curves
# For the ageing networks manuscript
# Written by Shikhara Bhat
# JGU Mainz, Germany
# Date created: Fri Apr 11 15:04:02 2025
# Date last modified: Mon Oct 06 17:25:33 2025
###########################
library(ggplot2) #for plotting
library(svglite) #for saving figures as vectorized SVGs (resolution free)


##########################
#Define the curves we want to plot

#Gompertz-Makeham curve
# mu(t) = M + a e^{bt}
gompertz_makeham <- function(t,phi,k,mu_e,chi,f0){
  
  #composite variables for convenience
  alpha <- phi + k*f0
  beta <- phi + k
  
  #from our analytics
  M <- mu_e - (chi*phi)/(k+phi)
  a <- chi*alpha/((k+phi)*(1-f0))
    
  return (M + a*exp(beta*t))
}

#define the solution to our logistic failures model
# the model is df/dt = (1-f)(r(f0) - kr*f0 + kr*f), approximating the more general model df/dt = (1-f)r(f)
#f0 is initial number of failures.
# r0 is r(f0)
#kr is dr/df at f = f0
#M is extrinsic mortality (analogous to Makeham parameter in GM)
logistic_fail_linear_mort <- function(t,phi,k,mu_e,chi,f0){
  
  #composite variables for convenience
  alpha <- phi + k*f0
  beta <- phi + k
  
  numerator <- (alpha*exp(beta*t))-phi*(1-f0)
  denominator <- beta + alpha*(exp(beta*t)-1)
  
  #analytical solution
  f <- numerator/denominator
  
  return (mu_e + chi*f)
}

logistic_fail_geometric_mort <- function(t,phi,k,mu_e,chi,f0){
  
  #composite variables for convenience
  alpha <- phi + k*f0
  beta <- phi + k
  
  numerator <- (alpha*exp(beta*t))-phi*(1-f0)
  denominator <- beta + alpha*(exp(beta*t)-1)
  
  #analytical solution
  f <- numerator/denominator
  
  intr_mort <- f/(1-f) 
  return (mu_e + chi*intr_mort)
}

##Find the plateau time
plateau_point <- function(phi,kr){
  return (0.5*(1-(phi/kr)))
}

plateau_time <- function(phi,f0,kr){
  beta <- phi + kr
  fstar <- plateau_point(phi,kr)
  log_num <- (1-f0)*(phi+kr*fstar)
  log_den <- (1-fstar)*(phi+kr*f0)
  log_term <- log(log_num/log_den)
  return( (1/beta)*log_term)
}

###############################
#Plotting

#define parameters
#X and Y axis scales are unparameterized/arbitrary anyway.
# my_params <- list(phi = 0.001, #intrinsic failure (constant component)
#                   k = 0.1, #interdependence rate (coefficient of linear component)
#                   mu_e = 1e-5, #extrinsic mortality
#                   chi = 0.2,   #how strongly failures impact mortality
#                   f0 = 0.01 #initial condition
# )


my_params <- list(phi = 0.001, #intrinsic failure (constant component)
                  k = 0.1, #interdependence rate (coefficient of linear component)
                  mu_e = 1e-5, #extrinsic mortality
                  chi = 1,   #how strongly failures impact mortality
                  f0 = 0.001 #initial condition
)


# my_params_linear <- list(phi = 0.01, #intrinsic failure (constant component)
#                          k = 0.1, #interdependence rate (coefficient of linear component)
#                          mu_e = 1e-5, #extrinsic mortality
#                          chi = 0.5,   #how strongly failures impact mortality
#                          f0 = 0.01 #initial condition
# )


  
#base layer. Control the axis limits and padding here
lineplot <- ggplot() + scale_x_continuous(limits = c(10,70), expand = c(0,0))
#lineplot <- ggplot() + scale_x_continuous(limits = c(0.1,0.85), expand = c(0,0))
lineplot <- lineplot + scale_y_continuous(limits = c(-0.5,2.3), expand = c(0,0))


###########################
#Plot on a linear scale

#add Gompertz curve with specified parameters
lineplot <- lineplot + geom_function(fun=gompertz_makeham, args = my_params,aes(color='Gompertz-Makeham'),linetype='dashed',linewidth=2)
#add logistic curve with specified parameters
lineplot <- lineplot + geom_function(fun=logistic_fail_linear_mort, args = my_params,aes(color='Logistic failure, linear mortality'),linewidth=1.5)
lineplot <- lineplot + geom_function(fun=logistic_fail_geometric_mort, args = my_params,aes(color='Logistic failure, geometric mortality'),linewidth=1.5, alpha = 0.6)
#add a vertical line for plateau time
#lineplot <- lineplot + geom_vline(xintercept = plateau_time(phi=0.08,f0=0,kr=7.92),linetype='dotted',linewidth=2,alpha=0.6, color='#666666')

#add a legend
lineplot <- lineplot + scale_color_manual(values =  c('Gompertz-Makeham'='#000000',
                                                      'Logistic failure, linear mortality'='#40B0A6',
                                                     # 'Linear failure, linear mortality'='#cc9922',
                                                      'Logistic failure, geometric mortality'='#e17f1e'))


# #aesthetics
lineplot <- lineplot + theme_light() + xlab('Age') + ylab("Mortality rate")
lineplot <- lineplot + theme(panel.border = element_blank(), axis.line = element_line(colour = "black",linewidth=1), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
#lineplot <- lineplot + theme(axis.text = element_text(face = 'bold', color = 'black',size = 30)) #in case we want numbers on the axes
lineplot <- lineplot + theme(axis.title = element_text(face = 'bold',color = 'black',size = 35))
lineplot <- lineplot + theme(legend.text = element_text(face = 'bold', color = 'black',size = 10))
lineplot <- lineplot + theme(aspect.ratio = 1) #set aspect ratio
lineplot <- lineplot + theme(legend.title=element_blank()) #remove legend title
lineplot <- lineplot + theme(axis.ticks.length = unit(.15, "cm")) + theme(axis.ticks = element_line(color='black',linewidth=0.5))
#plot margins are to be provided in the order top, right, bottom, left
lineplot <- lineplot + theme(plot.margin = margin(0.1, 1, ,0.8 , "cm"))

#Remove the numeric values on the x and y axes (because we don't care about the actual numbers)
lineplot <- lineplot + theme(axis.text.x=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank())

lineplot


#uncomment below to save the figure as an SVG
ggsave(filename = './output/fig_2A.svg',plot=lineplot,width=12,height=12)



#################################
#also plot the same curves on a log scale

#transform to log scale
logplot <- lineplot + scale_y_continuous(trans='log') + ylab("log mortality rate")
logplot

#uncomment below to save the figure as an SVG
ggsave(filename = './output/fig_2B.svg',plot=logplot,width=12,height=12)