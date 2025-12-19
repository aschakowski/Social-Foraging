################################################################################
#
# Title: Compute raw features
#
# Description: Computes raw feature distances based on feature matrix
#
# Authors: Schakowski, A.
#
# Last updated: 13/02/2024 (DD/MM/YYYY)
#
################################################################################
extract_features <- function(gps_data, spot_selection_data){
############### HELPER FUNCTIONS ###############################################
# Extract social distances
extract_social <- function(sub, gps_data, time_option){
  
  # subset spot identifier
  day_tmp = sub[1, "day"]
  year_tmp = sub[1, "year"]
  id_tmp = sub[1, "camera_id"]
  coord_tmp = sub[, c("x", "y")]
  time_tmp = sub[1, "timestamp_start"]
  
  # subset coordinates of other participants
  coordinates = gps_data %>% 
    filter(day == day_tmp &
             year == year_tmp &
             camera_id != id_tmp &
             timestamp == time_tmp - time_option) %>% 
    drop_na(position_lat, position_long) %>% 
    ungroup() %>% 
    dplyr::select(x = position_lat,y = position_long, ID = camera_id)
  
  # one decision being made at the very end of competition, where there is no
  # GPS data of other participants. in this case, take gps positions from 
  # 1 min prior.
  if (nrow(coordinates) < 15){
    coordinates = gps_data %>% 
      filter(day == day_tmp &
               year == year_tmp &
               camera_id != id_tmp &
               timestamp == time_tmp - 1*60) %>% 
      drop_na(position_lat, position_long) %>% 
      ungroup() %>% 
      dplyr::select(x = position_lat,y = position_long, ID = camera_id)
  }
  
  coordinates_xy = LongLatToUTM(coordinates$y, coordinates$x, zone = 35)
  coordinates_xy$ID = coordinates$ID
  
  # compute distances to target
  tmp_distances = fastPdist2(as.matrix(coord_tmp), as.matrix(coordinates_xy[, 2:3]))
  
  return(tmp_distances)
  
}

# extract roughness (distance to curve)
extract_roughness <- function(sub, depth_profile_points, lake_key){
  
  sub = left_join(sub, lake_key, by = c("day", "year"))
  points = st_coordinates(depth_profile_points[[sub$lake_index[1]]])
  
  coord_tmp = sub[, c("x", "y")]
  depth_profile = depth_profiles[[sub$lake_index[1]]]
  
  # compute distances to target
  d = fastPdist2(as.matrix(coord_tmp), as.matrix(points[, 1:2]))
  
  # if very small or 0 sometimes na
  d[which(is.na(d))] <- .00001
  
  # keep only 50 closest features for storage
  d = t(apply(d, 1, sort))[, 1:50]
  return(d)
  
}

# additional eco features (REVISION)
extract_ecological_features <- function(sub, depth_profiles,lake_key){
  
  sub = left_join(sub, lake_key, by = c("day", "year"))
  coord_tmp = sub[, c("x", "y")]
  profile = depth_profiles[[sub$lake_index[1]]]
  
  pts = list()
  for (i in 1:nrow(profile)){
    pts[[i]] = st_line_sample(profile[i,], density = 1/10)
    pts[[i]] = pts[[i]] %>% st_cast("POINT")
    pts[[i]] = st_coordinates(pts[[i]])
    pts[[i]] = as.data.frame(pts[[i]])
    pts[[i]]$depth = pull(profile[i, "depth"] %>% st_drop_geometry())
  }
  pts = do.call(rbind, pts)
  colnames(pts) <- c("x", "y", "depth")
  
  pts = pts[!is.na(pts$x) & !is.na(pts$y) & !is.na(pts$depth),]
  
  # depth model linear interpolation: 
  gs <- gstat(formula = depth ~ 1, locations= ~x+y, data = pts, nmax = 50, set = list(idp = 0))
  
  # simulate N = 100 points from disk of radius r=100
  N_pts = 100
  # simulate points using polar coordinate transformation
  roughness_matrix = matrix(nrow = nrow(sub), ncol = 4)
  for (pt in 1:nrow(sub)){
    
    # simulate random angle and radius within disk
    angles = runif(N_pts, -pi, pi)
    r = runif(N_pts, 0, 100)
    
    # simulate N points from location
    eval_points_x = sub$x[pt] + r * cos(angles)
    eval_points_y = sub$y[pt] + r * sin(angles)
    eval_df = data.frame(x = eval_points_x, y = eval_points_y)  
    
    # depth
    pred = predict(gs, eval_df)
    
    # turn into proper sf dataframe
    sf_df = st_as_sf(data.frame(x = eval_points_x, y = eval_points_y), coords = c("x", "y"))
    st_crs(sf_df) <- crs(profile)
    
    # distances to isoclines
    dists = st_distance(sf_df, profile)
    dist_to_isocline = apply(dists, 1, min)
    
    # compute variance in depth
    roughness_matrix[pt,1] = sd(pred$var1.pred, na.rm = T)
    roughness_matrix[pt,2] = max(pred$var1.pred)-min(pred$var1.pred)
    roughness_matrix[pt,3] = mean(dist_to_isocline, na.rm = T)
    roughness_matrix[pt,4] = mean(pred$var1.pred, na.rm = T)
    
  }
  
  colnames(roughness_matrix) = c("var_depth", "max_min_depth", "dist_iso", "mean_depth")
  
  return(roughness_matrix)

}


# extract alternative roughness (variation depth within radius) (not used in paper)
extract_alternative_roughness <- function(sub, lake_key){
  
  sub = left_join(sub, lake_key, by = c("day", "year"))
  #points = st_coordinates(depth_profile_points[[sub$lake_index[1]]])
  
  coord_tmp = sub[, c("x", "y")]
  depth_profile = depth_profiles[[sub$lake_index[1]]]
  
  # simulate N = 100 points from disk of radius r
  N_pts = 100
  radii = seq(from = 0, to = 700, by = 100)
  roughness_matrix = matrix(nrow = nrow(sub), ncol = length(radii))
  c=1
  for (radius in radii){
    
    # simulate points using polar coordinate transformation
    for (pt in 1:nrow(sub)){
      
      # simulate random angle and radius within disk
      angles = runif(N_pts, -pi, pi)
      r = runif(N_pts, 0, radius)
      
      # simulate N points from location
      eval_points_x = sub$x[pt] + r * cos(angles)
      eval_points_y = sub$y[pt] + r * sin(angles)
      
      # turn into proper sf dataframe
      sf_df = st_as_sf(data.frame(x = eval_points_x, y = eval_points_y), coords = c("x", "y"))
      st_crs(sf_df) <- crs(depth_profile)
      
      # compute depth for each point
      depth = depth_profile[apply(st_distance(sf_df, depth_profile), 1, which.min), ]$depth
      
      # compute variance in depth
      roughness_matrix[pt, c] = sd(depth)
      
    }
    # count to next radius
    c=c+1
  }
  
  return(roughness_matrix)

}

# extract alternative roughness 2 (max-min depth within radius) (not used in paper)
extract_alternative_roughness2 <- function(sub, lake_key){
  
  sub = left_join(sub, lake_key, by = c("day", "year"))
  #points = st_coordinates(depth_profile_points[[sub$lake_index[1]]])
  
  coord_tmp = sub[, c("x", "y")]
  depth_profile = depth_profiles[[sub$lake_index[1]]]
  
  # simulate N = 100 points from disk of radius r
  N_pts = 100
  radii = seq(from = 0, to = 700, by = 100)
  roughness_matrix = matrix(nrow = nrow(sub), ncol = length(radii))
  c=1
  for (radius in radii){
    
    # simulate points using polar coordinate transformation
    for (pt in 1:nrow(sub)){
      
      # simulate random angle and radius within disk
      angles = runif(N_pts, -pi, pi)
      r = runif(N_pts, 0, radius)
      
      # simulate N points from location
      eval_points_x = sub$x[pt] + r * cos(angles)
      eval_points_y = sub$y[pt] + r * sin(angles)
      
      # turn into proper sf dataframe
      sf_df = st_as_sf(data.frame(x = eval_points_x, y = eval_points_y), coords = c("x", "y"))
      st_crs(sf_df) <- crs(depth_profile)
      
      # compute depth for each point
      depth = depth_profile[apply(st_distance(sf_df, depth_profile), 1, which.min), ]$depth
      
      # compute variance in depth
      roughness_matrix[pt, c] = max(depth) - min(depth)
      
    }
    # count to next radius
    c=c+1
  }
  
  return(roughness_matrix)
  
}

# extract alternative roughness 3 (Average dist. to curves within radius) (not used in paper)
extract_alternative_roughness3 <- function(sub, lake_key){
  
  sub = left_join(sub, lake_key, by = c("day", "year"))
  
  coord_tmp = sub[, c("x", "y")]
  depth_profile = depth_profiles[[sub$lake_index[1]]]
  
  # simulate N = 100 points from disk of radius r
  N_pts = 100
  radii = seq(from = 0, to = 700, by = 100)
  roughness_matrix = matrix(nrow = nrow(sub), ncol = length(radii))
  c=1
  for (radius in radii){
  
    # simulate points using polar coordinate transformation
    for (pt in 1:nrow(sub)){
      
      # simulate random angle and radius within disk
      angles = runif(N_pts, -pi, pi)
      r = runif(N_pts, 0, radius)
      
      # simulate N points from location
      eval_points_x = sub$x[pt] + r * cos(angles)
      eval_points_y = sub$y[pt] + r * sin(angles)
      
      # turn into proper sf dataframe
      sf_df = st_as_sf(data.frame(x = eval_points_x, y = eval_points_y), coords = c("x", "y"))
      st_crs(sf_df) <- crs(depth_profile)
      
      # compute depth for each point
      dist = apply(st_distance(sf_df, depth_profile), 1, min)
      
      # compute variance in depth
      roughness_matrix[pt, c] = mean(dist)
      
    }
    # count to next radius
    c=c+1
  }
  
  return(roughness_matrix)
  
}

# extract depth (not used in paper)
extract_depth <- function(sub, lake_key){
  
  sub = left_join(sub, lake_key, by = c("day", "year"))
  
  coord_tmp = sub[, c("x", "y")]
  depth_profile = depth_profiles[[sub$lake_index[1]]]
  
  # simulate N = 100 points from disk of radius r
  N_pts = 100
  radii = seq(from = 0, to = 700, by = 100)
  depth_matrix = matrix(nrow = nrow(sub), ncol = length(radii))
  c=1
  for (radius in radii){
    
    # simulate points using polar coordinate transformation
    for (pt in 1:nrow(sub)){
      
      # simulate random angle and radius within disk
      angles = runif(N_pts, -pi, pi)
      r = runif(N_pts, 0, radius)
      
      # simulate N points from location
      eval_points_x = sub$x[pt] + r * cos(angles)
      eval_points_y = sub$y[pt] + r * sin(angles)
      
      # turn into proper sf dataframe
      sf_df = st_as_sf(data.frame(x = eval_points_x, y = eval_points_y), coords = c("x", "y"))
      st_crs(sf_df) <- crs(depth_profile)
      
      # compute depth for each point
      depth = depth_profile[apply(st_distance(sf_df, depth_profile), 1, which.min), ]$depth
      
      # compute variance in depth
      depth_matrix[pt, c] = mean(depth)
      
    }
    # count to next radius
    c=c+1
  }
  
  return(depth_matrix)
  
}

# compute distance to successes and losses
extract_spot_distance <- function(sub, spot_selection, option){ #option either 0 or 1
  
  # define cpp function for cluster export
  #sourceCpp("functions/fast_dist.cpp")
  
  # subset spot identifier
  day_tmp = sub[1, "day"]
  year_tmp = sub[1, "year"]
  id_tmp = sub[1, "camera_id"]
  coord_tmp = sub[, c("x", "y")]
  spot_id_tmp = sub[1, "spot_id"]
  
  # get coordinates of last spots
  spots = spot_selection %>% 
    ungroup() %>% 
    filter(day == day_tmp & year == year_tmp & camera_id == id_tmp & spot_id<spot_id_tmp & Return == option) %>% 
    dplyr::select(x, y, spot_id)
  
  if(nrow(spots) == 0){
    result = rep(-99, nrow(sub)) # Important: choices and non-choices get same 
    # value if no info available and will thus not affect estimates 
    # (subjected to softmax, each spot gets same choice prob).
  } else {
    result = fastPdist2(as.matrix(coord_tmp), as.matrix(spots[, c("x", "y")]))
    
    # sometimes missing, if distances very close or 0, happens very rarely (<.01%).
    result[is.na(result)]<-.0001
    
    # keep only 50 closest features for storage
    if (dim(result)[2] == 1){
      result = result
    } else {
      result = t(apply(result, 1, sort))[, 1:min(ncol(result), 50)] 
    }
    
  }
  
  return(result)
  
}

# extract distance to last spot
extract_locality_feature <- function(sub, spot_selection_data){
  
  # define cpp function for cluster export
  #sourceCpp("functions/fast_dist.cpp")
  
  # subset spot identifier
  day_tmp = sub[1, "day"]
  year_tmp = sub[1, "year"]
  id_tmp = sub[1, "camera_id"]
  coord_tmp = sub[, c("x", "y")]
  spot_id_tmp = sub[1, "spot_id"]
  
  # get coordinates of last spots
  spots = spot_selection_data %>% 
    ungroup() %>% 
    filter(day == day_tmp & year == year_tmp & camera_id == id_tmp & spot_id<spot_id_tmp) %>% 
    dplyr::select(x, y, spot_id, Return)
  
  if(nrow(spots) == 0){
    result = rep(-99, nrow(sub)) 
  } else {
    result = fastPdist2(as.matrix(coord_tmp), as.matrix(spots[nrow(spots), c("x", "y")]))
  }
  
  # sometimes missing, if distances very close (rarely, <.01%)
  result[is.na(result)]<-0.0001
  
  return(result)
  
}

# extract angle to last spot (not used in paper)
extract_angle_feature <- function(sub, spot_selection_data){
  
  # define cpp function for cluster export
  #sourceCpp("functions/fast_dist.cpp")
  # subset spot identifier
  day_tmp = sub[1, "day"]
  year_tmp = sub[1, "year"]
  id_tmp = sub[1, "camera_id"]
  coord_tmp = sub[, c("x", "y")]
  spot_id_tmp = sub[1, "spot_id"]
  
  # get coordinates of last spots
  spots = spot_selection_data %>% 
    ungroup() %>% 
    filter(day == day_tmp & year == year_tmp & camera_id == id_tmp & spot_id<spot_id_tmp) %>% 
    dplyr::select(x, y, spot_id, Return)
  
  if(nrow(spots) <= 1){
    result = rep(-99, nrow(sub))
  } else {
    result = vector()
    for (row in 1:nrow(coord_tmp)){
      result[row] = momentuHMM:::turnAngle(as.vector(as.matrix(spots[nrow(spots)-1, c("x", "y")])),
                             as.vector(as.matrix(spots[nrow(spots), c("x", "y")])), 
                             as.vector(as.matrix(coord_tmp[row,])))
    }
  }
  
  # sometimes missing, if distances very close
  result[is.na(result)]<-0
  
  return(result)
  
}

################################################################################
# Load feature matrix
Feature_matrix = readRDS(file = "utils/data/processed_data/feature_matrix_spot_selection.rds")

# load identifiers
identifiers = readRDS("utils/data/processed_data/identifiers_spot_selection.rds")

# define cpp function for cluster export
sourceCpp("utils/library/fast_dist.cpp")

################################################################################
# extract social distances
################################################################################
# at time of joining (time option = 0, results also hold for negative offsets, i.e., 
# x minutes before joining)
if (file.exists("utils/data/processed_data/social_distance_matrix.rds")){
  
} else {
  social_distance_matrix = map(Feature_matrix, function(x) extract_social(x, gps_data, time_option = 0), .progress = TRUE)
  
  # save social distances
  saveRDS(social_distance_matrix, file = "utils/data/processed_data/social_distance_matrix.rds")
  
  rm(social_distance_matrix)
  gc()
  gc()
  
}

################################################################################
# Extract Roughness
################################################################################
# Import Lake Data
lake_key = data.frame(day = c(1:5, 1:3, 5,4), 
                      year = rep(c(2023, 2022), each = 5),
                      lake_index = 1:10)

lake_data = readRDS("utils/data/raw_data/lake_data/depth_profiles.rds")
competition_area = readRDS("utils/data/raw_data/lake_data/study_lakes.rds")

# get in same order
depth_profiles = list()
depth_profiles[[1]] <- lake_data[[7]]
depth_profiles[[2]] <- lake_data[[10]]
depth_profiles[[3]] <- lake_data[[6]]
depth_profiles[[4]] <- lake_data[[9]]
depth_profiles[[5]] <- lake_data[[8]]
depth_profiles[[6]] <- lake_data[[1]]
depth_profiles[[7]] <- lake_data[[2]]
depth_profiles[[8]] <- lake_data[[3]]
depth_profiles[[9]] <- lake_data[[5]]
depth_profiles[[10]] <- lake_data[[4]]

# depth profile points
depth_profile_points = list()
for (i in 1:10){
  tmp = st_line_sample(depth_profiles[[i]] %>% filter(depth>0), density = .1, type = "regular")
  depth_profile_points[[i]] = tmp[!st_is_empty(tmp)] %>% st_cast("POINT")
}


if (file.exists("utils/data/processed_data/roughness_distance_matrix.rds")){
  
} else {
  # extract feature
  roughness_distance_matrix = map(Feature_matrix, function(x) extract_roughness(x, depth_profile_points, lake_key), .progress = TRUE)
  
  # save result
  saveRDS(roughness_distance_matrix, file = "utils/data/processed_data/roughness_distance_matrix.rds")
  rm(roughness_distance_matrix)
  gc()
}

# ALTERNATIVE WAYS TO QUANTIFY ROUGHNESS AND DEPTH (REVISION):
if (file.exists("utils/data/processed_data/alternative_roughness_distance_matrix.rds")){
  
} else {
  
  cl <- makeClusterPSOCK(35)
  future::plan(future::cluster, workers = cl, gc = T)
  alternative_roughness_distance_matrix = future_map(Feature_matrix, function(x) extract_ecological_features(x, depth_profiles, lake_key),
                         .progress = TRUE,
                         .options = furrr_options(seed = TRUE))
  parallel::stopCluster(cl)
  
  # save result
  saveRDS(alternative_roughness_distance_matrix, file = "utils/data/processed_data/alternative_roughness_distance_matrix.rds")
  rm(alternative_roughness_distance_matrix)
  gc()
}
################################################################################
# Extract spot distances
################################################################################
# extract feature
if (file.exists("utils/data/processed_data/success_distance_matrix.rds")){
  
} else {
  
  success_distance_matrix = map(Feature_matrix, function(x) extract_spot_distance(x, spot_selection_data, option = 1), .progress = TRUE)
  saveRDS(success_distance_matrix, file = "utils/data/processed_data/success_distance_matrix.rds")
}

if (file.exists("utils/data/processed_data/loss_distance_matrix.rds")){
  
} else {
  
  loss_distance_matrix = map(Feature_matrix, function(x) extract_spot_distance(x, spot_selection_data, option = 0), .progress = TRUE)
  saveRDS(loss_distance_matrix, file = "utils/data/processed_data/loss_distance_matrix.rds")
}

################################################################################
if (file.exists("utils/data/processed_data/locality_matrix.rds")){
  
} else {
  # extract locality feature
  locality_matrix = list()
  for (i in 1:length(Feature_matrix)){
    locality_matrix[[i]] <- extract_locality_feature(Feature_matrix[[i]], spot_selection_data)
  }
  
  saveRDS(locality_matrix, file = "utils/data/processed_data/locality_matrix.rds")
  
  rm(locality_matrix)
  gc()
  
}

if (file.exists("utils/data/processed_data/angle_matrix.rds")){
  
} else {
  # extract angle feature (not used in paper)
  angle_matrix = list()
  for (i in 1:length(Feature_matrix)){
    angle_matrix[[i]] <- extract_angle_feature(Feature_matrix[[i]], spot_selection_data)
  }
  
  saveRDS(angle_matrix, file = "utils/data/processed_data/angle_matrix.rds")
  
  rm(angle_matrix)
  gc()

}

################################################################################
# END
################################################################################
}