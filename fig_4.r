############################
# Compute expected survivorship at the point where the deviation between
# the model with disappearance and the model without disappearance first
# exceeds a fixed threshold
# Written by Shikhara Bhat
# JGU Mainz, Germany
# Date: Sat Mar 28 23:27:01 2026
###########################

library(dplyr) #for data wrangling
library(parallel) #for parallelisation
library(ggplot2) #for plotting

set.seed(314) #for reproducibility

save_figs <- F #should you save the plots to file or just view them?
#########################
#define parameters

#failure parameters
# we are using failures: f ---> f + 1/N at rate r_+ = (1-f)*(2phi + kr f)
#              repair:   f ---> f - 1/N at rate r_- = (1-f)*phi
#
# thus the drift is r_+ - r_- = (1-f)*(phi + kr f)
# and the diffusion is (r_+ + r_-)/N = (1-f)*(3phi + kr)/N.
phi <- 0.001 #intrinsic failure (constant component)
k <- 0.1 #interdependence rate (coefficient of linear component)

#mortality parameters. We will use the mortality rule mu_e + chi*f/(1-f)
mu_e <- 1e-5 #extrinsic mortality
chi<- 0.2   #how strongly failures impact mortality

f0 <- 0.01 #initial condition

# simulation parameters
n_paths <- 20000 #number of realizations to sum over (larger is better)
n_steps <- 1000 #discretisation of the time interval (larger is better, but very time and space/memory intensive)
T_total <- 500 #total time for which to calculate our quantities of interest

dt <- T_total / n_steps #discretisation of the time (will be used for the Brownian motion)


# parameter grid to loop over
N_vals <- 10^seq(2, 5, length.out = 25) #system size parameter values. Making this finer will make the code slower
theta_vals <- seq(0.01, 1, length.out = 100) #threshold difference in deviation, \delta_mu/E[\mu(\mathfrak{f})]


########################
#Define functions that define the model (failure rule, mortality rule, turnover)
# We will simulate the SDE dxt = r(x_t)dt +sqrt{tau(x_t)}dW_t

# failure, repair, and mortality
# we are using failures: f ---> f + 1/N at rate r_+ = (1-f)*(2phi + kr f)
#              repair:   f ---> f - 1/N at rate r_- = (1-f)*phi
#
# so that we have   r_+ - r_-   =   (1-f)*(phi + kr f)
#     and          (r_+ + r_-)/N   = (1-f)*(3phi + kr)/N.


r <- function(x){#expected failure rate r_+ - r_-
  return ((1-x)*(phi + k*x))
}

tau <- function(x,N){#variance in failure rate (r_+ + r_-)/N. Note that in the main text, we define tau as r_+ + r_- and the SDE has sqrt{tau}, so they are equivalent.
  return ((1-x)*(3*phi + k*x)/N)
}

mu <- function(x){#mortality
  
  #finite order truncation in case x/(1-x) gets too large for R
  # return (mu_e + chi*(x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 + x^8 + x^9 + x^10)) 
  
  return (mu_e + chi*(x/(1-x)))
} 

#keep the process in [0,1]
reflect_at_boundaries <- function(x) {
  if (x < 0) {#0 boundary
    return (-x)
  } else if (x > 1) {#1 boundary
    return (2 - x)
  } else {
    return (x)
  }
}

#we vectorise this function when parallelising
reflect_at_boundaries <- Vectorize(reflect_at_boundaries)

# put it all together. 
#This will loop over many thetas for fixed N and return expected survivorship
#at the point where \delta first crosses each threshold \theta
#unlike feynmankac_expectations.R, we do NOT store the entire path here,
#we just keep the expected survivorship at the threshold point. This makes
#simulation much more memory efficient.
Sbar_for_one_N <- function(N, theta_vals) {
  
  ft <- rep(f0, n_paths) #start with initial condition
  integral_mu <- numeric(n_paths)
  
  Sbar_at_crossing <- rep(NA, length(theta_vals)) #survivorship at the crossing point. Initiate at NA so we know if it never crossed.
  theta_crossed <- rep(FALSE, length(theta_vals)) #boolean. check if we've crossed the detection threshold
  
  #simulate
  for (step in 1:n_steps) {
    
    # the Wiener process
    dWt <- rnorm(n_paths, sd = sqrt(dt))
    
    # Simulate Tithonus SDE using the Euler-Maruyama method.
    ft <- reflect_at_boundaries(ft + r(ft) * dt + sqrt(tau(ft,N)) * dWt)
    
    # compute mortality and survivorship
    mu_vals <- mu(ft) 
    integral_mu <- integral_mu + mu_vals * dt
    surv <- exp(-integral_mu)
    
    # expectations
    Sbar <- mean(surv)
    Expmu <- mean(mu_vals)
    ExpmuSbar <- mean(mu_vals * surv)
    
    #the covariance process Cov[X,Y] \coloneqq E[XY] - E[X]E[Y]
    cov <- ExpmuSbar - Expmu*Sbar
   
    #\delta_\mu, from Feynman-Kac
    delta <- -cov/Sbar
    proportional_delta <- delta/(Expmu) #scale by Tithonus mortality
    
    # check if proportional_delta has crossed our threshold
    crossed_thetas <- which((!theta_crossed) & (proportional_delta >= theta_vals))
    if (length(crossed_thetas) > 0) {
      Sbar_at_crossing[crossed_thetas] <- Sbar
      theta_crossed[crossed_thetas] <- TRUE
    }
    
    # stop the loop if we've crossed all the relevant thresholds
    if ((all(theta_crossed))) {
      break
    }
  }
  
  return(Sbar_at_crossing)
}

####################
#parallelisation

my_clusters  <- makeCluster(detectCores()-1,
                            type='FORK') #use type 'FORK' if on a Mac or Unix system and 'PSOCK' if on windows

#import all the variables into the clusters
clusterExport(my_clusters, varlist = c("r", "tau", "mu", "reflect_at_boundaries",
                                       "phi","k","mu_e","chi","f0",
                                       "n_paths", "n_steps","T_total", "dt", "theta_vals"))

#set seed
clusterSetRNGStream(my_clusters, 314)

#parallelise and loop over both N and theta
results_list <- parLapply(my_clusters, N_vals, Sbar_for_one_N, theta_vals)

stopCluster(my_clusters)
####################

#collate the results
results <- do.call(cbind, results_list)

#put the data in a format that ggplot likes
expand.grid(theta = theta_vals,N = N_vals) %>%
  mutate(S = as.vector(results)) -> plotdata

############################
#Plotting 
#for formatting the axis labels on log plots. From https://stackoverflow.com/a/10036038
fmt_dcimals <- function(decimals=0,scientific=FALSE){
  function(x) format(x,nsmall = decimals,scientific = scientific)
}

#make the exponent labels look prettier. Slight modification of https://stackoverflow.com/a/18530540
pretty_powers_of_10 <- function(x) {
  parse(text=gsub("1e", "10^", scales::scientific_format()(x)))
}


#make the plot
heatmap_plot <- ggplot(plotdata, aes(x = N, y = theta, fill = S)) + geom_tile()

heatmap_plot <- heatmap_plot + scale_fill_viridis_c(name = "Expected\nsurvivorship",
                     option = 'viridis',
                    breaks = c(1e-5,1e-4,1e-3,1e-2,1e-1,0.5,1),
                     transform =scales::pseudo_log_trans(sigma =1e-12),
                     na.value = "grey", #NA means either we never crossed the threshold or survivorship is under the limit
                     limits = c(1e-5,1),
                     #oob = scales::squish, #this uses the scales package to squish out of bound points
                     labels=c(pretty_powers_of_10(c(1e-5,1e-4,1e-3,1e-2)),0.1,0.5,1)
)

#convert X axis to log scale 
heatmap_plot <- heatmap_plot +  scale_x_continuous(labels = fmt_dcimals(0),transform = 'log10')

#aesthetics
heatmap_plot <- heatmap_plot + theme_light() + xlab("Number of \nsub-systems (N)") + ylab("Detection Threshold \n(theta)")
heatmap_plot <- heatmap_plot + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
heatmap_plot <- heatmap_plot + theme(axis.text = element_text(face = 'bold', color = 'black',size = 30))
heatmap_plot <- heatmap_plot + theme(axis.title = element_text(face = 'bold',color = 'black',size = 35))
#failplot <- failplot + theme(legend.position = 'none') #remove legend
#legend formatting
heatmap_plot <- heatmap_plot + theme(
  legend.key.height = unit(1, "null"),
  legend.margin = margin(0, 0, 0, 0),
  legend.text = element_text(face = 'bold',color = 'black', size = 30),
  legend.title = element_text(face = 'bold',color = 'black', size = 30)
)
heatmap_plot <- heatmap_plot + theme(aspect.ratio = 1) #set aspect ratio
heatmap_plot <- heatmap_plot + theme(axis.ticks.length = unit(.15, "cm")) + theme(axis.ticks = element_line(color='black',linewidth=0.5))
#plot margins are to be provided in the order top, right, bottom, left
heatmap_plot <- heatmap_plot + theme(plot.margin = margin(0.8, 0.5, 0.5 ,0.5 , "cm"))

#add param values in title
heatmap_plot <-  heatmap_plot + ggtitle(paste('chi = ',chi, ', mu_e = ',mu_e,'\nphi = ',phi,',kr = ',k,sep='')) 
heatmap_plot<- heatmap_plot + theme(title = element_text(face='bold',size=15))

#zoom in and make x and y limits right
heatmap_plot <- heatmap_plot + coord_cartesian(xlim = c(min(N_vals),max(N_vals)),ylim = c(0,1),expand=FALSE)  

#view the plot
heatmap_plot


if(save_figs){
  
  #save as SVG
  ggsave('./output/fig_4.svg',heatmap_plot,width=12,height=12)
  
}