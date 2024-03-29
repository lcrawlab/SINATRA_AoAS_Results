#!/usr/bin/env Rscript

### Clear Console ###
cat("\014")

### Clear Environment ###
#rm(list = ls(all = TRUE))
### Load in the R Libraries ###
library(truncnorm)
library(doParallel)
library(svd)
library(numbers)

### Load SINATRA functions ###
source("SINATRA_Code/ec_computation.R")
source("SINATRA_Code/generate_directions.R")
source("SINATRA_Code/gp_inference.R")
source("SINATRA_Code/metric_curve_simulation.R")
source("SINATRA_Code/plotting_functions.R")
source("SINATRA_Code/roc_curve_simulation.R")
source("SINATRA_Code/shape_reconstruction.R")
source("SINATRA_Code/simulated_data_generation.R")
source("SINATRA_Code/RATEv2.R")
sourceCpp("SINATRA_Code/BAKRGibbs.cpp")

### Set the parameters for the analysis ###
set.seed(4913, kind = "L'Ecuyer-CMRG")
n.simulations <- 50

# take this input from command line
num_causal_region <- 1
num_shared_region <- 2
causal_points <- 10
shared_points <- 10


######################################################################################
######################################################################################
######################################################################################

### Setup DoParallel ###
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores, type="FORK")  
registerDoParallel(cl)  

lengths <- 10*(1:10)

### Run the analysis in Parallel ###
simulation_results <- foreach(i=1:n.simulations, .combine = 'rbind', .noexport = c('GaussKernel')) %:% 
  foreach(j=1:10, .combine = 'rbind', .noexport = c('GaussKernel')) %dopar% {
    
    set.seed(10*i+j)
    
    res <- tryCatch( generate_ROC_with_coned_directions(nsim = 25, curve_length = lengths[j], grid_size = 25, distance_to_causal_point = 0.1, 
                                                        causal_points = causal_points,shared_points = shared_points, desired_num_cones = 20, eta = 0.1, 
                                                        truncated = 300, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                        min_points = 3,directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                        mode = 'sphere', fpr = 0.05, start = 1, cusps = 50,
                                                        subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region),
                     error = function(x) {
                       return(matrix(nrow = 0,ncol = 3))
                     }
    )
    ### Label the results for each trial and directions ###
    rdf <- cbind(res, rep(lengths[j], dim(res)[1]) )
    rdf <- data.frame(cbind(rdf,rep(i,dim(rdf)[1])))
    rdf <- plyr::rename(rdf,c("X1" = "FPR","X2" = "TPR","X3" = "Class","X4" = "Index","X5" = "Lengths","X6" = "Trial"))
  }

stopCluster(cl)

######################################################################################
######################################################################################
######################################################################################
### Aggregate results ###
rdfmeans <- aggregate(simulation_results[c("FPR","TPR")],
                      by = list("Lengths" = simulation_results$Lengths,
                                "Index" = simulation_results$Index,
                                "Class" = simulation_results$Class), mean)

rdfmeans$Lengths <- as.factor(rdfmeans$Lengths)

## Plot results ###
ROC_curve_plt <- ggplot(data <- rdfmeans[rdfmeans$Class == 1,],aes(x = FPR, y = TPR, color = Lengths)) +
 geom_line(stat = "identity") +
 labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "Lengths") +
# ggtitle(sprintf("Curve Length Sensitivity - Hard Case")) +
 geom_abline(intercept = 0, slope = 1) +
 coord_equal(ratio=1) +
 theme_bw() +
 theme(plot.title = element_text(hjust = 0.5))
print(ROC_curve_plt)
######################################################################################
######################################################################################
######################################################################################
### save the results ###

# save the raw results
sim_results_file = sprintf("~/data/SINATRA/new/data_ROC_curve_length_causal%d_shared%d_region_size%s.RData",
                           num_causal_region,num_shared_region,causal_points)
save(simulation_results, file = sim_results_file)

# save the dataframe
df_results_file = sprintf("~/data/SINATRA/new/df_ROC_curvelength_causal%d_shared%d_%s.RData",
                          num_causal_region,num_shared_region,shared_points)
save(rdfmeans, file = df_results_file)