set.seed(4913, kind = "L'Ecuyer-CMRG")
library(sinatra)
library(FNN)
library(rgl)
library(Rvcg)
library(plyr)
library(pdist)
library(gglasso)
library(numbers)
library(data.table)
setwd('~/Documents')

#Parameters for the Analysis

cap_radius = 0.15
num_cones = 5
directions_per_cone = 5
len = 75

num_causal_region = 2
num_shared_region = 1
causal_points = 10
shared_points = 10

subdivision = 3

nsim = 25

cusps = 50
causal_dirs = generate_equidistributed_points(cusps, cusps +1)
causal_regions_1 = sample(1:cusps,num_causal_region)
causal_regions_2 = sample((1:cusps)[-causal_regions_1],num_causal_region)
shared_regions = sample((1:cusps)[-c(causal_regions_1,causal_regions_2)],num_shared_region)

runs = 100

#### Fig 2.a. ####
g1 = generate_averaged_ROC_with_coned_directions(runs = runs, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                   mode = 'sphere_baseline',num_cusps = cusps,
                                   subdivision = subdivision,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = max, write = TRUE)


write.csv(g1,'~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_elastic_net_max2.csv')
g1 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_elastic_net_max2.csv')[,-1][,-4]
roc_curve_frame1.1 = data.frame(g1)
library(ggplot2)
roc_curve_frame1.1 = roc_curve_frame1.1[roc_curve_frame1.1[,3] == 1,]
roc_curve_frame1.1[,3] = rep('Baseline (EN Max)',dim(roc_curve_frame1.1)[1])
g1.2 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                 causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                 truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                 min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                 mode = 'sphere_baseline',num_cusps = cusps,
                                                 subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = mean, write = FALSE)
#
#
write.csv(g1.2,'~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_elastic_net_mean2.csv')
g1.2 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_elastic_net_mean2.csv')[,-1][,-4]
roc_curve_frame1.2 = data.frame(g1.2)
library(ggplot2)
roc_curve_frame1.2 = roc_curve_frame1.2[roc_curve_frame1.2[,3] == 1,]
roc_curve_frame1.2[,3] = rep('Baseline (EN Mean)',dim(roc_curve_frame1.2)[1])
g1.3 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                 causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                 truncated = 400, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                 min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                 mode = 'sphere_baseline',num_cusps = cusps,
                                                 subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = max, write = FALSE)

#
write.csv(g1.3,'~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_ridge_max2.csv')
g1.3 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_ridge_max2.csv')[,-1][,-4]
roc_curve_frame1.3 = data.frame(g1.3)
library(ggplot2)
roc_curve_frame1.3 = roc_curve_frame1.3[roc_curve_frame1.3[,3] == 1,]
roc_curve_frame1.3[,3] = rep('Baseline (RR Max)',dim(roc_curve_frame1.3)[1])
g1.4 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                 causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                 truncated = 400, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                 min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                 mode = 'sphere_baseline',num_cusps = cusps,
                                                 subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = mean, write = FALSE)


write.csv(g1.4,'~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_ridge_mean2.csv')
g1.4 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared_ridge_mean2.csv')[,-1][,-4]
roc_curve_frame1.4 = data.frame(g1.4)
library(ggplot2)
roc_curve_frame1.4 = roc_curve_frame1.4[roc_curve_frame1.4[,3] == 1,]
roc_curve_frame1.4[,3] = rep('Baseline (RR Mean)',dim(roc_curve_frame1.4)[1])
roc_curve_frame = rbind(roc_curve_frame1.1,roc_curve_frame1.2,roc_curve_frame1.3,roc_curve_frame1.4)
library(ggplot2)
load('~/gitrepos/SINATRA/Simulations/Results/new/df_ROC_causal1_shared2_10.RData')
rdfmeans = rdfmeans[rdfmeans$Class == 1,]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[rdfmeans$Num_Cones==60,]
rdfmeans = rdfmeans[c('FPR','TPR')]
rdfmeans[,3] = rep('SINATRA',dim(rdfmeans)[1])
names(rdfmeans) = names(roc_curve_frame)
roc_curve_frame = rbind(roc_curve_frame,rdfmeans)
ggplot(roc_curve_frame, aes(x = V1,y = V2,group = V3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(V3) )) +
  geom_line(stat = "identity") +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("2 Causal Regions, 1 Shared Regions, Size 10")) +
  geom_abline(intercept = 0, slope = 1) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
ggsave('~/Documents/SINATRA/Scripts/Data/baseline_2causal_1_shared2.pdf')

#### Fig 2.b. ####
num_causal_region = 6
num_shared_region = 3
causal_points = 10
shared_points = 10
g2 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                 causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                 truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                 min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                 mode = 'sphere_baseline',num_cusps = cusps,
                                                 subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = max, write = FALSE)


write.csv(g2,'~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_elastic_net_max2.csv')
g2 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_elastic_net_max2.csv')[,-1][,-4]
roc_curve_frame2.1 = data.frame(g2)
library(ggplot2)
roc_curve_frame2.1 = roc_curve_frame2.1[roc_curve_frame2.1[,3] == 1,]
roc_curve_frame2.1[,3] = rep('Baseline (EN Max)',dim(roc_curve_frame2.1)[1])
g2.2 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                   mode = 'sphere_baseline',num_cusps = cusps,
                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = mean, write = FALSE)


write.csv(g2.2,'~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_elastic_net_mean2.csv')
g2.2 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_elastic_net_mean2.csv')[,-1][,-4]
roc_curve_frame2.2 = data.frame(g2.2)
library(ggplot2)
roc_curve_frame2.2 = roc_curve_frame2.2[roc_curve_frame2.2[,3] == 1,]
roc_curve_frame2.2[,3] = rep('Baseline (EN Mean)',dim(roc_curve_frame2.2)[1])
g2.3 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                   mode = 'sphere_baseline',num_cusps = cusps,
                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = max, write = FALSE)


write.csv(g2.3,'~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_ridge_max2.csv')
g2.3 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_ridge_max2.csv')[,-1][,-4]
roc_curve_frame2.3 = data.frame(g2.3)
library(ggplot2)
roc_curve_frame2.3 = roc_curve_frame2.3[roc_curve_frame2.3[,3] == 1,]
roc_curve_frame2.3[,3] = rep('Baseline (RR Max)',dim(roc_curve_frame2.3)[1])
g2.4 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                   mode = 'sphere_baseline',num_cusps = cusps,
                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = mean, write = FALSE)


write.csv(g2.4,'~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_ridge_mean2.csv')
g2.4 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared_ridge_mean2.csv')[,-1][,-4]
roc_curve_frame2.4 = data.frame(g2.4)
library(ggplot2)
roc_curve_frame2.4 = roc_curve_frame2.4[roc_curve_frame2.4[,3] == 1,]
roc_curve_frame2.4[,3] = rep('Baseline (RR Mean)',dim(roc_curve_frame2.4)[1])
roc_curve_frame = rbind(roc_curve_frame2.1,roc_curve_frame2.2,roc_curve_frame2.3,roc_curve_frame2.4)
library(ggplot2)
load('~/gitrepos/SINATRA/Simulations/Results/new/df_ROC_causal3_shared6_10.RData')
rdfmeans = rdfmeans[rdfmeans$Class == 1,]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[rdfmeans$Num_Cones==60,]
rdfmeans = rdfmeans[c('FPR','TPR')]
rdfmeans[,3] = rep('SINATRA',dim(rdfmeans)[1])
names(rdfmeans) = names(roc_curve_frame)
roc_curve_frame = rbind(roc_curve_frame,rdfmeans)
ggplot(roc_curve_frame, aes(x = V1,y = V2,group = V3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(V3) )) +
  geom_line(stat = "identity") +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("6 Causal Regions, 3 Shared Regions, Size 10")) +
  geom_abline(intercept = 0, slope = 1) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
ggsave('~/Documents/SINATRA/Scripts/Data/baseline_6causal_3_shared2.pdf')
#### Fig 2.c. ####
num_causal_region = 10
num_shared_region = 5
causal_points = 10
shared_points = 10
g3 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                 causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                 truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                 min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                 mode = 'sphere_baseline',num_cusps = cusps,
                                                 subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = max, write = FALSE)


write.csv(g3,'~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_elastic_net_max2.csv')
g3 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_elastic_net_max2.csv')[,-1][,-4]
roc_curve_frame3.1 = data.frame(g3)
library(ggplot2)
roc_curve_frame3.1 = roc_curve_frame3.1[roc_curve_frame3.1[,3] == 1,]
roc_curve_frame3.1[,3] = rep('Baseline (EN Max)',dim(roc_curve_frame3.1)[1])
g3.2 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                   mode = 'sphere_baseline',num_cusps = cusps,
                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 0.5,reduce = mean, write = FALSE)


write.csv(g3.2,'~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_elastic_net_mean2.csv')
g3.2 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_elastic_net_mean2.csv')[,-1][,-4]
roc_curve_frame3.2 = data.frame(g3.2)
library(ggplot2)
roc_curve_frame3.2 = roc_curve_frame3.2[roc_curve_frame3.2[,3] == 1,]
roc_curve_frame3.2[,3] = rep('Baseline (EN Mean)',dim(roc_curve_frame3.2)[1])
g3.3 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                   mode = 'sphere_baseline',num_cusps = cusps,
                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = max, write = FALSE)


write.csv(g3.3,'~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_ridge_max2.csv')
g3.3 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_ridge_max2.csv')[,-1][,-4]
roc_curve_frame3.3 = data.frame(g3.3)
library(ggplot2)
roc_curve_frame3.3 = roc_curve_frame3.3[roc_curve_frame3.3[,3] == 1,]
roc_curve_frame3.3[,3] = rep('Baseline (RR Max)',dim(roc_curve_frame3.3)[1])
g3.4 = generate_averaged_ROC_with_coned_directions(runs = 100, nsim = 50, curve_length = 30, grid_size = 25, distance_to_causal_point = 0.1, 
                                                   causal_points = causal_points,shared_points = shared_points, num_cones = 20, eta = 0.1, 
                                                   truncated = 200, two_curves = TRUE, ball = TRUE, ball_radius = 1.5, type = 'vertex',
                                                   min_points = 3, directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                                   mode = 'sphere_baseline',num_cusps = cusps,
                                                   subdivision = 3,num_causal_region = num_causal_region, num_shared_region = num_shared_region,alpha = 1,reduce = mean, write = FALSE)


write.csv(g3.4,'~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_ridge_mean2.csv')
g3.4 = read.csv('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared_ridge_mean2.csv')[,-1][,-4]
roc_curve_frame3.4 = data.frame(g3.4)
library(ggplot2)
roc_curve_frame3.4 = roc_curve_frame3.4[roc_curve_frame3.4[,3] == 1,]
roc_curve_frame3.4[,3] = rep('Baseline (RR Mean)',dim(roc_curve_frame3.4)[1])
roc_curve_frame = rbind(roc_curve_frame3.1,roc_curve_frame3.2,roc_curve_frame3.3,roc_curve_frame3.4)
library(ggplot2)
load('~/gitrepos/SINATRA/Simulations/Results/new/df_ROC_causal5_shared10_10.RData')
rdfmeans = rdfmeans[rdfmeans$Class == 1,]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[,-2]
rdfmeans = rdfmeans[rdfmeans$Num_Cones==60,]
rdfmeans = rdfmeans[c('FPR','TPR')]
rdfmeans[,3] = rep('SINATRA',dim(rdfmeans)[1])
names(rdfmeans) = names(roc_curve_frame)
roc_curve_frame = rbind(roc_curve_frame,rdfmeans)
ggplot(roc_curve_frame, aes(x = V1,y = V2,group = V3)) + geom_line(alpha = 0.8, size = 2,aes(color = factor(V3) )) +
  geom_line(stat = "identity") +
  labs(x = "FPR (False Positive Rate)", y = "TPR (True Positive Rate)", color = "# Cones") +
  ggtitle(sprintf("10 Causal Regions, 5 Shared Regions, Size 10")) +
  geom_abline(intercept = 0, slope = 1) + 
  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) 
ggsave('~/Documents/SINATRA/Scripts/Data/baseline_10causal_5_shared2.pdf')

#### Test New Function ####
library(glmnet)
data2 = generate_data_sphere_simulation_new(nsim = nsim,dir = dirs, curve_length = len,noise_points = shared_points,
                                       causal_points = causal_points,ball_radius = ball_radius, subdivision = subdivision,
                                       cusps = cusps, causal_regions_1 = causal_regions_1, causal_regions_2 = causal_regions_2,
                                       shared_regions = shared_regions, ec_type = ec_type)
groups = rep(1:(dim(data2$data)[2]/3),each = 3)
#lasso = cv.gglasso(x = data2$data[,-1], y = data2$data[,1], group = groups, loss = "logit")
lasso = cv.glmnet(x = data2$data[,-1], y = data2$data[,1], alpha = 0.0,  family = "binomial")
coefs  = coef(lasso, s = 'lambda.min')

#Binary Option
nonzero = which(abs(coefs[-1])>0.8)
results = unique(groups[nonzero])


# Heatmap Option


g = cbind(groups, abs(coefs[-1]))

library(dplyr)
df = as.data.table(g)

new_df = aggregate(df[,2],list(df$groups),mean)

abs_coefs = new_df$V2
results = groups[which(abs_coefs > 0.5)]
cuts = cut(abs_coefs,50,labels = FALSE)


color1='blue'
color2='lightgreen'
color3='red'
#col_pal = c("#00007F", "blue", "#007FFF", "cyan",  "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000")
col_pal=c(color1,color1,color2,color2,color2,color3)
col_pal=c(color1,color2,color3)

colfunc <- colorRampPalette(col_pal)
heat_colors=colfunc(max(cuts) - min(cuts))[cuts - min(cuts)]

mesh = vcgSphere(subdivision = 3)
mesh$vb[1:3,] = t(data2$complex_points[[1]])
cols = rep('white', dim(mesh$vb)[2])
cols[data$causal_points1] = 'red'
cols[data$causal_points2] = 'red'
cols[data$noise] = 'blue'
cols[results] = 'green'

mfrow3d(1,2)
plot3d(mesh, col = cols)
plot3d(mesh, col = heat_colors)
#### Sanity Checks ####

j = cor(t(data_summary$data[,-1]))
heatmap(j)


generate_ROC_with_coned_directions(nsim = 10, curve_length = 50, grid_size = 25, distance_to_causal_point = 0.1, 
                                   causal_points = causal_points,shared_points = shared_points,  eta = 0.1, 
                                   truncated = -1, two_curves = TRUE, ball = TRUE, ball_radius = 2.5, type = 'vertex',
                                   min_points = 3,directions_per_cone = 5, cap_radius = 0.15, radius = 1,ec_type = 'ECT',
                                   mode = 'sphere_baseline', 
                                   subdivision = 3,num_causal_region = 2, num_shared_region = 1)


#### Comp 2 ####
indices = c(261,295,833,224,939,278,252,593,177,357,4845,4215,4619,162,244,4618,4225,31,4705,4524,219,4301,4898,262,4860,69,64,4478,4452,4224,102,179,4176)

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

col_pal2=c(color1,color1,color2,color2,color3)
colfunc2 <- colorRampPalette(col_pal2)

