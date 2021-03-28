
######################################################################################
######################################################################################
######################################################################################

### Metric Curve Code ###

generate_metric_curve <- function(nsim = 20, curve_length = 20, grid_size = 30, distance_to_causal_point = 0.1, causal_points = 8,
                                  shared_points = 5, num_cones = 15, directions_per_cone = 4, eta = 0.1, ball_radius = 2.5, type){


  print("generating directions")
  # generate directions, length num_directions
  initial_cones <- 150
  directions <- generate_equidistributed_cones(initial_cones,0.1,directions_per_cone)
  num_cones = dim(directions)[1]/(directions_per_cone)


  print("generating data")
  # generate data
  data <- create_data_normal_fixed(num_sim = nsim, dir = directions, curve_length = curve_length,shared_points = shared_points,
                                   causal_points = causal_points,grid_size = grid_size,eta = eta,ball_radius = ball_radius)

  # Prune directions
  print("pruning data")
  temp <- prune_directions_to_desired_number(data = data$data[,-1],directions, initial_cones,curve_length,directions_per_cone,num_cones)
  directions <- temp[[1]]
  ec_curve_data <- temp[[2]]
  num_cones <- dim(directions)[1]/directions_per_cone

  metric_curve <- matrix(0, nrow = num_cones, ncol = 2)

  #want to vary directions, add on a new cone at each step
  for (i in 1:num_cones){
    print(sprintf("Analyzing Direction %i", i))

    # Generate the Rate using the data; take only data corresponding to desired directions
    print("computing rate values")
    rate_values <- find_rate_variables_with_other_sampling_methods(ec_curve_data[,1:(i*curve_length*directions_per_cone)],bandwidth = 0.01,type = 'ESS')[,2]

    print("computing metrics")
    if (type == 'vertex'){
      metrics <- compute_metrics_vertex(data_points = data$complex_points, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                        curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                        eta = eta, directions_per_cone = directions_per_cone, directions = directions[1:(i*directions_per_cone),], ball_radius = ball_radius,
                                        ball = ball)

    }
    if (type == 'feature'){
      metrics <- compute_metrics_feature(data_points = data$complex_points, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                         curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                         eta = eta, directions_per_cone = directions_per_cone, dir = directions[1:(i*directions_per_cone),], ball = ball,ball_radius = ball_radius,
                                         min_points = min_points)
    }
    if (type == 'cone'){
      metrics <- compute_metrics_cone(data_points = data$complex_points, class_1_causal_points = data$causal_points1, class_2_causal_points = data$causal_points2,
                                      curve_length = curve_length, distance_to_causal_point = distance_to_causal_point, rate_values = rate_values, grid_size = grid_size,
                                      eta = eta, directions_per_cone = directions_per_cone, dir = directions[1:(i*directions_per_cone),],
                                      ball = ball, ball_radius = ball_radius,  min_points = min_points, radius = radius)
    }

    metric_curve[i,] <- metrics
  }

  return(metric_curve)
}


######################################################################################
######################################################################################
######################################################################################

# Computing Metrics
compute_metrics_vertex <- function(data_points,class_1_causal_points,class_2_causal_points,distance_to_causal_point = 0.1,
                                   rate_values,grid_size,eta = 0.1,directions_per_cone, curve_length,directions, ball_radius, ball = ball){

  num_vertices = grid_size^2
  #Initializing the aggregate ROC curve frame
  total_metric = matrix(0, nrow = length(data_points),ncol = 2)



  # go down the list of complexes?
  for (i in 1:length(data_points)){

    #Interpolating based on the causal and shared points in R^3 for each shape
    predictions=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=data_points[[i]],eta=eta)
    complex=matrix_to_simplicial_complex(predictions,grid_length=grid_size)

    #Starting to Compute the ROC curve for a given complex
    class_1_true_vertices = c()
    class_2_true_vertices = c()

    for (j in 1:num_vertices){
      #computes the 2D euclidean distance on the grid between the points
      dist1=apply(X = class_1_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:2])
      dist2=apply(X = class_2_causal_points[,1:2],MARGIN = 1,FUN = difference,y=complex$Vertices[j,1:2])

      if (min(dist1)< distance_to_causal_point) class_1_true_vertices=c(class_1_true_vertices,j)
      if (min(dist2)< distance_to_causal_point) class_2_true_vertices=c(class_2_true_vertices,j)
    }
    combined_true_vertices = union(class_1_true_vertices,class_2_true_vertices)

    class_1_false_vertices = setdiff(1:num_vertices, class_1_true_vertices)
    class_2_false_vertices = setdiff(1:num_vertices, class_2_true_vertices)
    combined_false_vertices = setdiff(1:num_vertices, combined_true_vertices)

    rate_ROC <- matrix(0,nrow = 1,ncol = 2)
    for (threshold in quantile(rate_values,probs = seq(1,0,length.out = floor(length(rate_values)/5)) ) ){

      rate_positive_vertices <- compute_selected_vertices_cones(dir = directions, complex = complex, rate_vals = rate_values,
                                                                len = curve_length, threshold = threshold,
                                                                cone_size = directions_per_cone, ball_radius = ball_radius)

      rate_negative_vertices <- setdiff(1:num_vertices,rate_positive_vertices)

      TPR_FPR <- calculate_TPR_FPR(rate_positive_vertices,rate_negative_vertices,
                                   class_1_true_vertices,class_1_false_vertices)
      rate_ROC <- rbind(rate_ROC, TPR_FPR)

      # if FPR > 0.05, break
      if(TPR_FPR[1] > 0.05) break();
    }

    metrics <- c(TPR_at_specified_FPR_metric(0.05,rate_ROC),
                 size_of_intersection_metric(combined_true_vertices, rate_positive_vertices))

    total_metric[i,] <- metrics
  }

  averaged_metrics <- colMeans(total_metric)

  return(averaged_metrics)
}


######################################################################################
######################################################################################
######################################################################################

### Helper Functions ###
#' TPR/FPR
#' @export
calculate_TPR_FPR <- function(positive_vertices, negative_vertices, true_vertices, false_vertices){
  TP <- length(intersect(positive_vertices, true_vertices))
  FP <- length(positive_vertices) - TP
  TN <- length(intersect(negative_vertices, false_vertices))
  FN <- length(negative_vertices) - TN

  TPR = TP/(TP + FN)
  FPR = FP/(FP + TN)

  return(c(FPR,TPR))
}

#Input an ROC curve: a matrix of size n x 2, where n is the threshold size.
# Find the (FPR,TPR) pair such that FPR < specified FPR.
TPR_at_specified_FPR_metric <- function(FPR, ROC){
  sorted_ROC <- ROC[order(ROC[,1]),]
  desired_index <- which(sorted_ROC[,1] >= FPR)[1]
  return(sorted_ROC[desired_index,2])
}

# This is extremely similar to a TPR; intersection over total union.
# Inputs are lists of vertex indices, 1 to n. Take causal points to be the union of class1,class2 causal points;
# selected points should be the output of the rate & cone recontruction idea.
size_of_intersection_metric <- function(causal_points,selected_points){
  int <- intersect(causal_points,selected_points)
  union <- union(causal_points,selected_points)
  return(length(int)/length(union))
}
