############################
# Compare how the deterministic predictions compare to
# expected trajectory of finite N (stochastic) systems that we got from Ito
# for the Tithonus model
# Written by Shikhara Bhat
# JGU Mainz, Germany
# Date: Mon Oct 06 20:33:41 2025
###########################
library(ggplot2) #for plotting


##########################
#Define the curves we want to plot
# we are using failures: f ---> f + 1/N at rate r_+ = (1-f)*(2phi + kr f)
#              repair:   f ---> f - 1/N at rate r_- = (1-f)*phi
#
# thus the failure rule is r_+ - r_- = (1-f)*(phi + kr f)
# and the turnover is (r_+ + r_-) = (1-f)*(3phi + kr).

#The Gompertz-Makeham curve M + a*exp(bt)
gompertz_makeham <- function(t,phi,k,mu_e,chi,f0){
  
  #composite variables for convenience
  alpha <- phi + k*f0
  beta <- phi + k
  
  #from our analytics
  M <- mu_e - (chi*phi)/(k+phi)
  a <- chi*alpha/((k+phi)*(1-f0))
  
  return (M + a*exp(beta*t))
}


# the N to infinity prediction
inf_N_solution <- function(t,phi,k,mu_e,chi,f0,N){
  
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


dmu <- function(t,phi,k,chi,f0,N){
  
  #composite variables for convenience
  alpha <- phi + k*f0
  beta <- phi + k
  
  numerator <- (alpha*exp(beta*t))-phi*(1-f0)
  denominator <- beta + alpha*(exp(beta*t)-1)
  
  #analytical solution
  f <- numerator/denominator
  
  r <- (1-f)*(phi+k*f)
  tau <- (1-f)*(3*phi + k)
  
  muprime <- chi/((1-f)^2)
  
  mudoubleprime <- 2*chi/((1-f)^3)
  
  
  return (r*muprime + (tau*mudoubleprime/(2*N)))
}

#for N < \infty we'll just solve it numerically
expected_finite_N_solution <- function(t,phi,k,mu_e,chi,f0,N){
  
  intrinsic_mort <- integrate(dmu, lower=0, upper = t, phi=phi,k=k,chi=chi,f0=f0,N=N)$value
  
  return (mu_e +  intrinsic_mort)
}

expected_finite_N_solution <- Vectorize(expected_finite_N_solution)


###############################
#Plotting

#define parameters
#X and Y axis scales are unparameterized/arbitrary anyway.
my_params <- list(phi = 0.001, #intrinsic failure (constant component)
                  k = 0.1, #interdependence rate (coefficient of linear component)
                  mu_e = 1e-5, #extrinsic mortality
                  chi = 0.2,   #how strongly failures impact mortality
                  f0 = 0.01 #initial condition
)

#base layer. Control the axis limits and padding here. Limits chosen here are arbitrary.
lineplot <- ggplot() + scale_x_continuous(limits = c(50,95), expand = c(0,0))
lineplot <- lineplot + scale_y_continuous(limits = c(-0.5,50), expand = c(0,0))


###########################
#Plot on a linear scale


#lineplot <- lineplot + geom_function(fun=gompertz_makeham, args = my_params,aes(color='Gompertz-Makeham'),linetype='dashed',linewidth=2,linetype = 'dashed')
lineplot <- lineplot + geom_function(fun=inf_N_solution, args = c(my_params,N=Inf),aes(color='N = infinity'),linewidth=1.5, alpha = 1)
lineplot <- lineplot + geom_function(fun=expected_finite_N_solution, args = c(my_params,N=200),aes(color='N = 200'),linewidth=1.5, alpha = 0.6)
lineplot <- lineplot + geom_function(fun=expected_finite_N_solution, args = c(my_params,N=100),aes(color='N = 100'),linewidth=1.5, alpha = 0.6)
lineplot <- lineplot + geom_function(fun=expected_finite_N_solution, args = c(my_params,N=500),aes(color='N = 500'),linewidth=1.5, alpha = 0.6)
lineplot <- lineplot + geom_function(fun=expected_finite_N_solution, args = c(my_params,N=50),aes(color='N = 50'),linewidth=1.5, alpha = 0.6)

#add a vertical line for plateau time
#lineplot <- lineplot + geom_vline(xintercept = plateau_time(phi=0.08,f0=0,kr=7.92),linetype='dotted',linewidth=2,alpha=0.6, color='#666666')

#add a legend
lineplot <- lineplot + scale_color_manual(values =  c('N = infinity'='#000000',
                                                      'N = 50'='#ff0000',
                                                      'N = 100'='#ff5500',
                                                      'N = 200'='#ff8800',
                                                      'N = 500'='#ffbb00'
                                                      ))


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
# 
# suppressWarnings(print(lineplot)) #view the plot
#
# #uncomment below to save the figure as an SVG
#ggsave(file='./output/fig_2C.svg', plot=lineplot, width=10, height=10)
 
# #################################
# #also plot the same curves on a log scale
# 
#transform to log scale
logplot <- lineplot + scale_y_continuous(trans='log') + ylab("log mortality rate")
logplot
# 
# #suppressWarnings(print(logplot)) #view the plot
# 
# #uncomment below to save the figure as an SVG
#ggsave(file='./output/fig_2D.svg', plot=logplot, width=10, height=10)
