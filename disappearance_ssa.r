############################
# Run Gillespie simulations illustrating stochastic failure accumulation with disappearance
# Written by Shikhara Bhat
# JGU Mainz, Germany
# Date: Tue Feb 24 17:03:47 2026
###########################

library(GillespieSSA2) #for stochastic simulations
library(dplyr) #for manipulating dataframes
library(parallel) #for parallelisation

setwd('/Users/shikharabhat/Desktop/work/PhD/projects/AmNat_aging_note/code/')

#change parameters as appropriate here 
model_params <- c(
  
  #failure parameters. We will use the effective failure rule (1-f)*(phi + k*f)  
  
  phi = 0.001, #intrinsic failure (constant component)
  k = 0.1, #interdependence rate (coefficient of linear component)
  
  
  #mortality parameters. We will use the mortality rule mu_e + chi*f/(1-f)
 mu_e = 1e-5, #extrinsic mortality
 chi = 0.2, #how strongly failures impact mortality
  # mu_e = 0, #extrinsic mortality
  # chi = 0, #how strongly failures impact mortality
  
  #system size parameter
  N = 200 #how many processes exist that can fail. Our state variable will go from 0 to N.
)


#Define total birth and death rates
reactions <- list(
  
  #                     rate                                    effect
  
  
  #failure
  
  reaction("(N-F1)*(2*phi + k*(F1/N))"     ,                    c(F1 = +1)),
  
  #repair
  reaction(       "(N-F1)*phi"             ,                   c(F1 = -1)),
  
  #mortality (jerry rigged: send the state to a large negative number, gillespieSSA2 will then immediately stop the simulation.)
  reaction("mu_e + (chi*(F1)/(N-F1))"      ,                   c(F1 = -1e9))
)

#Initial condition
initial_state <- c(
  F1 = round(as.integer(model_params['N']/100)) # 1% of systems failed to start with (arbitrary)
)


final_time <- 85     #time at which the Gillespie simulation stops
runs <- 9500         #no. of independent realizations to simulate


#Run the SSA

start_time = Sys.time() #for timing how long running the code takes

# Create cluster (leave 1 core free)
my_parallels <- makeCluster(detectCores() - 1,
                            type = 'FORK') #'FORK' only works on UNIX-like systems. Use PSOCK on Windows.

# Load package on each core
#clusterEvalQ(my_parallels, library(GillespieSSA2)) #uncomment if using PSOCK 

# Export the objects we defined above to each core
#uncomment if using PSOCK 
# clusterExport(my_parallels, c("initial_state", 
#                     "reactions",
#                     "model_params",
#                     "final_time"))


clusterSetRNGStream(my_parallels, iseed = 314) #set seeds on each cluster (for reproducibility)

# Run simulations in parallel
sim_list <- parLapply(my_parallels, 1:runs, function(i) {

  sim_run <- ssa(
    initial_state = initial_state,
    reactions = reactions,
    params = model_params,
    final_time = final_time,
    method = ssa_exact(),
    stop_on_neg_state = TRUE,
    census_interval = 0,
    sim_name = "disappearance"
  )

  data.frame(
    run  = i,
    time = sim_run$time,
    F1   = sim_run$state[, "F1"]
  )
})

stopCluster(my_parallels)

# Combine results
sim_data <- do.call(rbind, sim_list)


#Report how much time it took to run the simulations
running_time = Sys.time() - start_time

running_time

#Collate to dataframe

colnames(sim_data) <- c('run','time','F1')



#move from number of failed subsystems to proportion of failed subsystems
sim_data %>%
  filter(F1 >= 0) %>%  #remove dead/disappeared individuals at the point at which they died/disappearaed
  mutate(f = F1/(model_params['N'])) %>% #proportion of failed systems
  select(run,time,f) -> ssa_data

#save data to file
write.csv(ssa_data,paste('./output/ssa_data/N_',as.integer(model_params['N']),'_chi_',as.numeric(model_params['chi']),'_phi_',as.numeric(model_params['phi']),'_kr_',as.numeric(model_params['k']),'_num_runs_',max(sim_data$run),'.csv',sep=''),row.names=FALSE)

######################################
#plotting
# library(ggplot2) #for plotting
# 
# lineplot <- ggplot(ssa_data,aes(x=time,y=f,group=run)) + geom_line(color='darkgreen',alpha=0.6)
# 
# #x and y axis limits
# lineplot <- lineplot +  scale_x_continuous(limits=c(0,max(ssa_data$time)),expand=c(0,0)) + scale_y_continuous(limits=c(0,1),expand=c(0,0))
# 
# #aesthetics
# lineplot <- lineplot + theme_light() + xlab('Age (t)') + ylab("Failures (f)")
# lineplot <- lineplot + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# lineplot <- lineplot + theme(axis.text = element_text(face = 'bold', color = 'black',size = 30))
# lineplot <- lineplot + theme(axis.title = element_text(face = 'bold',color = 'black',size = 35))
# lineplot <- lineplot + theme(legend.position = 'none') #remove legend
# lineplot <- lineplot + theme(aspect.ratio = 1) #set aspect ratio
# lineplot <- lineplot + theme(axis.ticks.length = unit(.15, "cm")) + theme(axis.ticks = element_line(color='black',linewidth=0.5))
# #plot margins are to be provided in the order top, right, bottom, left
# lineplot <- lineplot + theme(plot.margin = margin(0.8, 0.5, 0.5 ,0.5 , "cm"))
# lineplot #view the plot
# 
# running_time #report running time
# 
# # #for saving animations
# library(gganimate)
# 
# animation_dir <- paste(paste('./output/ssa_anims/',runs,'_paths',sep=''),'/',sep='')
# 
# #create the directory if it does not exist yet
# if (!file.exists(animation_dir)){
#   dir.create(animation_dir)
# }
# 
# 
# setwd(animation_dir)
# my_params <- paste('N_',as.integer(model_params['N']),'_chi_',as.numeric(model_params['chi']),'_phi_',as.numeric(model_params['phi']),'_kr_',as.numeric(model_params['k']),'_num_runs_',max(ssa_data$run),sep='')
# #save the final state to file as a PNG
# ggsave(paste(my_params,'_final_state.png',sep=''),lineplot,width=10,height=10)
# 
# #animate the plot and save as a gif
# anim <- lineplot + transition_reveal(time)
# anim_save(paste(my_params,'_ssa.gif',sep=''),animation=animate(anim,nframes=150,height=800,width=800))

