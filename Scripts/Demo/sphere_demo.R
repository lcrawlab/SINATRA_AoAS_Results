rm(list=ls())

#### This is a demo of the SINATRA pipeline with perturbed spheres (as demonstrated in the paper) ####
library(sinatra)
#library(Rvcg)
#library(rgl)
#library(pracma)
set.seed(55)

#Specifying Parameters
num_causal_region <- 10
num_shared_region <- 5
causal_points <- 10
shared_points <- 10

# EC Curve Parameters
nsim = 25 
curve_length = 50
grid_size = 25
distance_to_causal_point = 0.1 
eta = 0.1 
truncated = 100
two_curves = FALSE
ball = TRUE
ball_radius = 1.5
type = 'vertex'
min_points = 3
directions_per_cone = 5
cap_radius = 0.25
radius = 1
ec_type = 'ECT'
mode = 'sphere'
fpr = 0.05
start = 1
cusps = 50
subdivision = 3
num_cones = 25



#### Generate Data ####
causal_dirs = generate_equidistributed_points(cusps, 
                                              cusps)
causal_regions_1 = sample(1:cusps, num_causal_region)
causal_regions_2 = sample((1:cusps)[-causal_regions_1], 
                          num_causal_region)
shared_regions = sample((1:cusps)[-c(causal_regions_1, 
                                     causal_regions_2)], num_shared_region)
directions <- generate_equidistributed_cones(num_cones, 
                                             cap_radius, directions_per_cone)

#Generating perturbed spheres & visualizing
data <- generate_data_sphere_simulation(nsim = nsim, 
                                        dir = directions, curve_length = curve_length, noise_points = shared_points, 
                                        causal_points = causal_points, ball_radius = ball_radius, 
                                        subdivision = subdivision, cusps = cusps, causal_regions_1 = causal_regions_1, 
                                        causal_regions_2 = causal_regions_2, shared_regions = shared_regions, 
                                        ec_type = ec_type, write = FALSE)

rate_values <- find_rate_variables_with_other_sampling_methods(data$data, 
                                                               bandwidth = 0.01, type = "ESS")[, 2]

dummy_sphere = vcgSphere()
data_points = data$complex_points
dummy_sphere$vb[1:3, ] =  t(data_points[[1]])
cols =  rep('white',dim(dummy_sphere$vb)[2])
cols[data$causal_points1] = 'red'
cols[data$noise] = 'blue'
plot3d(dummy_sphere, col = cols)

####
#### Plot ROC Curve, and visualize ####


remove = c()
counter = 0
truncated = 500

roc_curve1 = compute_roc_curve_vertex(data = data, 
                                      class_1_causal_points = data$causal_points1, 
                                      class_2_causal_points = data$causal_points2, 
                                      curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, 
                                      rate_values = rate_values, grid_size = grid_size, 
                                      eta = eta, directions_per_cone = directions_per_cone, 
                                      directions = directions, class = 1, truncated = truncated, 
                                      ball_radius = ball_radius, radius = radius, 
                                      mode = mode, subdivision = subdivision)


roc_curve1 = compute_roc_curve_vertex(data = data, 
                                      class_1_causal_points = data$causal_points1, 
                                      class_2_causal_points = data$causal_points2, 
                                      curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, 
                                      rate_values = rate_values, grid_size = grid_size, 
                                      eta = eta, directions_per_cone = directions_per_cone, 
                                      directions = directions, class = 1, truncated = truncated, 
                                      ball_radius = ball_radius, radius = radius, 
                                      mode = mode, subdivision = subdivision)
roc_curve2 = compute_roc_curve_vertex(data = data, 
                                      class_1_causal_points = data$causal_points1, 
                                      class_2_causal_points = data$causal_points2, 
                                      curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, 
                                      rate_values = rate_values, grid_size = grid_size, 
                                      eta = eta, directions_per_cone = directions_per_cone, 
                                      directions = directions, class = 2, truncated = truncated, 
                                      ball_radius = ball_radius, radius = radius, 
                                      mode = mode, subdivision = subdivision)
roc_curve1= cbind(roc_curve1, rep('Class 1', truncated))
roc_curve2= cbind(roc_curve2, rep('Class 2', truncated))
roc_curve = rbind(roc_curve1, roc_curve2)
roc_curve_frame = data.frame(roc_curve)
roc_curve_frame[,1] = as.numeric(as.character(roc_curve_frame[,1]))
roc_curve_frame[,2] = as.numeric(as.character(roc_curve_frame[,2]))
ggplot() + 
  geom_line(data = subset(roc_curve_frame, X3 == "Class 1"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5) +
  geom_line(data = subset(roc_curve_frame, X3 == "Class 2"), aes(x = X1,y = X2,group = X3,color = factor(X3)),alpha = 0.75,  size = 1.5) +
  # geom_line(stat = "identity",aes(color = factor(X3)), linetype = 'dotted') +
  labs(x = "False Positive Rate (FPR)", y = "True Positive Rate (TPR)", color = "Method") +
  ggtitle(sprintf("ROC Curve for 2 Causal cusps, 1 shared cusp")) +
  coord_cartesian(xlim= c(0,1.0))+
  #  geom_abline(intercept = 0, slope = 1, alpha = 0.5) + 
  #  coord_equal(ratio=1) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = 'bold'),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),axis.text=element_text(size=20),
        axis.title=element_text(size=20,face="bold"), legend.text = element_text(size = 18), legend.title = element_text(size=20,face="bold")) + scale_colour_hue(l=40)


#### Visualize the spheres ####

color1='blue'
color2='lightgreen'
color3='red'
col_pal=c(color1,color1,color2,color2,color2,color3)
colfunc <- colorRampPalette(col_pal)

mfrow3d(2,2)

dummy_sphere = vcgSphere()
dummy_sphere$vb[1:3, ] =  t(data_points[[1]])
dummy_sphere_mesh = convert_off_file(dummy_sphere)
cols =  rep('white',dim(dummy_sphere$vb)[2])
cols[data$causal_points1] = 'red'
cols[data$noise] = 'blue'
plot3d(dummy_sphere, col = cols, specular="black", axes = FALSE)

sphere_heat = reconstruct_vertices_on_shape(dir = directions,complex = dummy_sphere_mesh,
                                            rate_vals = rate_values,len = curve_length,cuts = 500,cone_size = directions_per_cone,ball_radius = ball_radius,ball = TRUE,radius = radius)
cols_reconstructed=colfunc(max(sphere_heat[,1]) - min(sphere_heat[,1]))[sphere_heat[,1] - min(sphere_heat[,1])]

plot3d(dummy_sphere, col = cols_reconstructed, specular="black", axes = FALSE)



dummy_sphere$vb[1:3, ] =  t(data_points[[2]])
dummy_sphere_mesh = convert_off_file(dummy_sphere)
cols =  rep('white',dim(dummy_sphere$vb)[2])
cols[data$causal_points2] = 'red'
cols[data$noise] = 'blue'
plot3d(dummy_sphere, col = cols, specular="black", axes = FALSE)

sphere_heat = reconstruct_vertices_on_shape(dir = directions,complex = dummy_sphere_mesh,
                                            rate_vals = rate_values,len = curve_length,cuts = 500,cone_size = directions_per_cone,ball_radius = ball_radius,ball = TRUE,radius = radius)
cols_reconstructed=colfunc(max(sphere_heat[,1]) - min(sphere_heat[,1]))[sphere_heat[,1] - min(sphere_heat[,1])]

plot3d(dummy_sphere, col = cols_reconstructed, specular="black", axes = FALSE)