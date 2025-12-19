################################################################################
#
# Title: Simulate spot selection
#
# Description: Simulates trajectories for each id based on posterior estimates
#
# Authors: Schakowski, A.
#
# Last updated: 02/04/2024 (DD/MM/YYYY)
#
################################################################################
simulate_spot_selection <- function(spot_selection_data, gps_data){

################################################################################
# compute gaussian kernel for kde2d
kern <- function(d, h, option, ncol){
  d = as.matrix(d[, 1:ncol])
  if(option == "social"){
    d = exp(-d^2 / (2 * h^2)) * (1.0 / (1 + exp(-20 * (d - 5))))
  } else {
    d = exp(-d^2 / (2 * h^2))
  }
  
  d = d %*% rep(1, ncol)
  
  return(d)
}

# softmax function to compute features
softmax <- function(par){
  n.par <- length(par)
  par1 <- sort(par, decreasing = TRUE)
  Lk <- par1[1]
  for (k in 1:(n.par-1)) {
    Lk <- max(par1[k+1], Lk) + log1p(exp(-abs(par1[k+1] - Lk))) 
  }
  val <- exp(par - Lk)
  return(val)
}

# function to compute heading angles
heading_angle <- function(x, y){
  v <- c(y[1]-x[1],y[2]-x[2])
  angle <- atan2(v[2],v[1])
  
  while(angle<=(-pi))
    angle <- angle + 2*pi
  while(angle>pi)
    angle <- angle - 2*pi
  
  return(angle)
}

# library
library(tidybayes)
library(Rcpp)

# source helper functions
sourceCpp("utils/library/fast_dist.cpp")
sourceCpp("utils/library/ptinpoly.cpp")

################################################################################
# load posterior estimates model
spot_selection_fit = readRDS(file = "utils/data/processed_data/model_fit/spot_selection_models/f_conditional.rds")

# spot selection data
spot_selection_data = fread("utils/data/raw_data/spot_selection_data.csv")
spot_selection_data = spot_selection_data %>% filter(!exclude)

# data
stan_data = readRDS(file = "utils/data/processed_data/stan_data_spot_selection.rds")

# compute individual level parameters
library(dplyr)
main_effects = spot_selection_fit %>% 
  spread_draws(weights[feature]) %>% 
  group_by(feature) %>% 
  filter(feature <=8) %>% 
  summarize(weights = mean(weights))
trip_effects = spot_selection_fit %>% 
  spread_draws(v_trip[trip, feature]) %>% 
  group_by(trip, feature) %>% 
  filter(feature <=8) %>% 
  summarize(weights = mean(v_trip))
lake_effects = spot_selection_fit %>% 
  spread_draws(v_lake[lake, feature]) %>% 
  group_by(lake, feature) %>% 
  filter(feature <=8) %>% 
  summarize(weights = mean(v_lake))
participant_effects = spot_selection_fit %>% 
  spread_draws(v_id[id, feature]) %>% 
  group_by(id, feature) %>% 
  filter(feature <=8) %>% 
  summarize(weights = mean(v_id))

# compute kernel parameters
bandwidths = spot_selection_fit %>% 
  spread_draws(lambdas[feature]) %>% 
  group_by(feature) %>% 
  summarize(h = mean(lambdas))
bandwidths_offset = spot_selection_fit %>% 
  spread_draws(v_lake[lake, feature]) %>% 
  filter(feature > 8 & feature < 11) %>% 
  group_by(lake, feature) %>% 
  summarize(offset = mean(v_lake)) %>% 
  mutate(feature = feature - 8)
bandwidths = left_join(bandwidths, bandwidths_offset)
bandwidths = bandwidths %>% 
  mutate(bandwidth = exp(h + offset))

# load ids
identifiers = readRDS(file = "utils/data/processed_data/identifiers_model_key.rds")

# convert time
identifiers$time_since_start = as.numeric(identifiers$time_since_start) / max(as.numeric(identifiers$time_since_start))
identifiers$day = as.numeric(identifiers$day)
identifiers$year = as.numeric(identifiers$year)
identifiers$spot_id = as.numeric(identifiers$spot_id)
identifiers$participant_id =as.numeric(identifiers$participant_id)
identifiers$track_id = as.numeric(identifiers$track_id)
identifiers$unique_spot_id = as.numeric(identifiers$unique_spot_id)
identifiers$lake_id = as.numeric(identifiers$lake_id)

# add to model data
spot_selection_data = left_join(spot_selection_data, identifiers)

# impute indices for first spot and remove entirely missing trips
spot_selection_data = spot_selection_data %>% 
  group_by(day, year, camera_id) %>% 
  mutate(missing_trip = ifelse(sum(is.na(trip_index)) == n(), 1, 0)) %>% 
  filter(!missing_trip) %>% 
  mutate(lake_index = ifelse(is.na(lake_index) & spot_id == 1, lake_index[2], lake_index),
         trip_index = ifelse(is.na(trip_index) & spot_id == 1, trip_index[2], trip_index),
         id_index = ifelse(is.na(id_index) & spot_id == 1, id_index[2], id_index),
         lake_index = ifelse(is.na(lake_index) & spot_id == 1, lake_index[2], lake_index))
spot_selection_data = spot_selection_data %>% arrange(trip_index, spot_id)

# add time since start
spot_selection_data = left_join(spot_selection_data, identifiers %>% dplyr::select(day, year, spot_id, camera_id, time_since_start, lake_index, trip_index, id_index))
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

# polygon lakes for faster execution time
point_lakes = list()
for (i in 1:10){
  point_lakes[[i]] = competition_area[[i]] %>% 
    st_buffer(-10) %>% 
    st_cast("POINT") %>% 
    as.data.frame()
  point_lakes[[i]] = point_lakes[[i]] %>% 
    mutate(x = unlist(map(point_lakes[[i]]$polygons,1)),
           y = unlist(map(point_lakes[[i]]$polygons,2))) %>% 
    dplyr::select(x,y)
}

################################################################################
# compute competition area
check = vector()
c = 1
competition_area_new=list()
for (i in 1:10){
  
  # subset lake area
  area = point_lakes[[i]]
  
  # coordinates of spot choices
  coords = spot_selection_data %>% filter(day == lake_key[i, "day"] & year == lake_key[i, "year"]) %>% dplyr::select(x, y) %>% filter(!is.na(x))
  
  # convex hull around spot choices
  hull = coords[coords %>% chull(),c("x", "y")]
  colnames(hull) <- c("x", "y")
  hull2 = rbind(area, hull)
  
  # hull encompassing lake area and spot choices
  hull3 = hull2[hull2 %>% chull(),c("x", "y")]
  
  # check whether any spots are outside competition area
  for (j in 1:nrow(coords)){
    check[c] = pnpoly(cbind(coords$x[j], coords$y[j]), as.matrix(hull3)) 
    c = c+1
  }
  
  # save hull
  competition_area_new[[i]] <- hull3
  
}

# 
depth_profile_points = list()
for (i in 1:10){
  tmp = st_line_sample(depth_profiles[[i]] %>% filter(depth>0), density = .1, type = "regular")
  depth_profile_points[[i]] = tmp[!st_is_empty(tmp)] %>% st_cast("POINT")
}

################################################################################
# run simulation
################################################################################
# number of alternative spots (same as in model)
N_spots = 50

# load identifiers
identifiers = readRDS("utils/data/processed_data/identifiers_model_key.rds")
identifiers = identifiers %>% 
  group_by(day, year, camera_id) %>% 
  mutate(lake_index = ifelse(is.na(lake_index), lake_index[2], lake_index),
         trip_index = ifelse(is.na(trip_index), trip_index[2], trip_index),
         id_index = ifelse(is.na(id_index), id_index[2], id_index))
################################################################################
simulate_trajectory <- function(iter=1,
                                spot_selection_data = spot_selection_data,
                                identifiers = identifiers,
                                competition_area_new=competition_area_new,
                                depth_profiles=depth_profiles,
                                gps_data=gps_data,
                                main_effects=main_effects,
                                trip_effects=trip_effects,
                                lake_effects=lake_effects,
                                participant_effects=participant_effects,
                                bandwidths=bandwidths,
                                N_spots=N_spots,
                                option = "brw"){
  
  # load c++ functions (needs to be here for future_map export)
  sourceCpp("utils/library/fast_dist.cpp")
  sourceCpp("utils/library/ptinpoly.cpp")
  
  # movement constraints
  steps = spot_selection_data$step[spot_selection_data$step > -99 & spot_selection_data$angle>-99]
  turning_angles = spot_selection_data$angle[spot_selection_data$step > -99 & spot_selection_data$angle>-99]
  
  c=1
  Result = list()
  
  for (trip_tmp in unique(spot_selection_data$trip_index)){
    
    # subset spots for this participant
    spot_sub = spot_selection_data %>% 
      filter(trip_index == trip_tmp)
    
    # select angling times
    times = spot_sub$timestamp_start
    
    # select successful spots
    return = spot_sub$Return
    
    # resample returns
    if (option != "random"){ #keep sequence for random to match to data later. otherwise only keep number
      return = sample(return, length(return), replace = F)
    }
    
    # time since start
    time_since_start = spot_sub$time_since_start
    
    # select coordinates of all spots
    coordinates = as.matrix(spot_sub[, c("x", "y")])
    
    # identifiers
    day_tmp = spot_sub$day[1]
    year_tmp = spot_sub$year[1]
    camera_id_tmp = spot_sub$camera_id[1]
    participant_tmp = pull(identifiers[identifiers$day == day_tmp & 
                                    identifiers$year == year_tmp & 
                                    identifiers$camera_id == camera_id_tmp,"id_index"])[1]
    lake_tmp = pull(identifiers[identifiers$day == day_tmp & 
                                    identifiers$year == year_tmp & 
                                    identifiers$camera_id == camera_id_tmp,"lake_index"])[1]
    trip_tmp = pull(identifiers[identifiers$day == day_tmp & 
                             identifiers$year == year_tmp & 
                             identifiers$camera_id == camera_id_tmp,"trip_index"])[1]
    
    # lake surface as points
    lake_points = competition_area_new[[lake_key[lake_key$day == day_tmp & lake_key$year == year_tmp, "lake_index"]]]
    
    # depth profile
    depth_profile = depth_profiles[[lake_key[lake_key$day == day_tmp & lake_key$year == year_tmp, "lake_index"]]]
    
    # social coordinates
    coordinates_social = list()
    counter = 1
    for (j in times){
      coordinates_social_tmp = gps_data %>% 
        filter(day == day_tmp &
                 year == year_tmp &
                 camera_id != camera_id_tmp &
                 timestamp == j) %>% 
        drop_na(position_lat, position_long) %>% 
        ungroup() %>% 
        dplyr::select(x = position_lat,y = position_long, camera_id, timestamp_start = timestamp)
      
      # sometimes if decisions close to end, no gps available (n = 1)
      if (nrow(coordinates_social_tmp) == 0){
        coordinates_social_tmp = gps_data %>% 
          filter(day == day_tmp &
                   year == year_tmp &
                   camera_id != camera_id_tmp &
                   timestamp == j-1*60) %>% 
          drop_na(position_lat, position_long) %>% 
          ungroup() %>% 
          dplyr::select(x = position_lat,y = position_long, camera_id, timestamp_start = timestamp)
      }
      coordinates_social[[counter]] = coordinates_social_tmp
      counter = counter + 1
    }
    coordinates_social = do.call(rbind, coordinates_social)
    
    ################################################################################
    # preallocate spot choices: start simulating from initial spot
    spot_choices = matrix(nrow = length(times),
                          ncol = 2)
    spot_choices[1,1:2] <- c(coordinates[1,1],coordinates[1,2])
    
    # preallocate feature vector
    social_features_choice = vector()
    success_features_choice = vector()
    roughness_features_choice = vector()
    locality_features_choice = vector()
    loss_features_choice = vector()
    time_since_last_success_choice = vector()
    check = vector()
    chosen_angle = vector()
    chosen_radius = vector()
    
    # simulate alternative spots
    if (nrow(spot_sub)<2){
      
    } else {
      
     for (n in 2:nrow(spot_sub)){
     
      # compute heading angle
      if (n >= 3){
        heading = heading_angle(x = spot_choices[n-2,], y = spot_choices[n-1,])
      } else {
        heading = sample(turning_angles, 1, replace = T)
      }
      
       # sample turning angle and step size
      angle = sample(turning_angles, N_spots)
      radius = sample(steps, N_spots)
      
      # add to coordinates of last spot
      simulated_x = spot_choices[n-1,1] + radius * cos(heading + angle)
      simulated_y = spot_choices[n-1,2] + radius * sin(heading + angle)
      
      # check if points are on lake surface, if not, resample  
      for (pt in 1:length(simulated_x)){
        counter = 1
        while(!pnpoly(c(simulated_x[pt], simulated_y[pt]), as.matrix(lake_points)) & counter < 200){
          
          sampled_index = sample(1:length(turning_angles), 1)
          angle[pt] = turning_angles[sampled_index]
          radius[pt] = steps[sampled_index]
          
          simulated_x[pt] =  spot_choices[n-1,1] + radius[pt] * cos(heading + angle[pt])
          simulated_y[pt] =  spot_choices[n-1,2] + radius[pt] * sin(heading + angle[pt])
          counter = counter + 1
        }
        
      }
      
      # alternative coordinates
      coord_tmp = cbind(simulated_x, simulated_y)
      
      # compute features for each alternative
      ####### Social Feature ########
      coordinates_social_sub = coordinates_social[coordinates_social$timestamp_start == times[n],]
      if (nrow(coordinates_social_sub) == 0){
        coordinates_social_sub = coordinates_social[coordinates_social$timestamp_start == (times[n]-1*60),]
      }
      
      coordinates_xy = LongLatToUTM(coordinates_social_sub$y, coordinates_social_sub$x, zone = 35)
      social_distances  = fastPdist2(as.matrix(coord_tmp), as.matrix(coordinates_xy[,2:3]))
      social_distances[is.na(social_distances)] <- sample(social_distances, 1)
      social_distances = t(apply(social_distances, 1, sort))[, 1:min(nrow(coordinates_xy), 15)]
      social_feature = kern(social_distances, h = pull(bandwidths[bandwidths$lake == lake_tmp & bandwidths$feature == 1, "bandwidth"]), option = "social", ncol = min(nrow(coordinates_xy), 15))
        
      ####### Roughness Feature ######
      d = fastPdist2(as.matrix(coord_tmp), as.matrix(lake_points[, 1:2]))
      d[which(is.na(d))] <- .00001
      d = t(apply(d, 1, sort))
      roughness_feature = d[,1]/1000
      
      ######## Success Feature #######
      if (is.null(dim(spot_choices[1:(n-1),]))){
        successful_spots = cbind(t(as.matrix(spot_choices[1:(n-1),])), return[1:(n-1)])
        successful_spots = successful_spots[successful_spots[, 3] == 1, 1:2]
      } else {
        successful_spots = cbind(spot_choices[1:(n-1),], return[1:(n-1)])
        successful_spots = successful_spots[successful_spots[, 3] == 1, 1:2]
      }
    
      time_since_last_success = ifelse(length(successful_spots)>0, n-max(which(return[1:(n-1)]==1)), NA)
      
      if(length(successful_spots) == 0){
        success_feature = rep(0, (N_spots))
      } else {
        
        if (is.null(dim(successful_spots))){
          result = fastPdist2(as.matrix(coord_tmp), t(as.matrix(successful_spots)))
        } else {
          result = fastPdist2(as.matrix(coord_tmp), as.matrix(successful_spots))
        }
        
        # sometimes missing, if distances very close
        result[is.na(result)]<-0
        
        # keep only 50 closest features for storage
        if (dim(result)[2] == 1){
          result = result
        } else {
          result = t(apply(result, 1, sort))[, 1:min(ncol(result), 15)] 
        }
        
        success_feature = kern(result, h = pull(bandwidths[bandwidths$lake == lake_tmp & bandwidths$feature == 2, "bandwidth"]), option = "success", ncol = min(ncol(result), 15))
        #success_feature = kern(result, h = 8, option = "success", ncol = min(ncol(result), 15))
        
      }
      
      ######## Loss Feature #######
      if (is.null(dim(spot_choices[1:(n-1),]))){
        unsuccessful_spots = cbind(t(as.matrix(spot_choices[1:(n-1),])), return[1:(n-1)])
        unsuccessful_spots = unsuccessful_spots[unsuccessful_spots[, 3] == 0, 1:2]
      } else {
        unsuccessful_spots = cbind(spot_choices[1:(n-1),], return[1:(n-1)])
        unsuccessful_spots = unsuccessful_spots[unsuccessful_spots[, 3] == 0, 1:2]
      }
      
      if(length(unsuccessful_spots) == 0){
        loss_feature = rep(0, N_spots)
      } else {
        
        if (is.null(dim(unsuccessful_spots))){
          result = fastPdist2(as.matrix(coord_tmp), t(as.matrix(unsuccessful_spots)))
        } else {
          result = fastPdist2(as.matrix(coord_tmp), as.matrix(unsuccessful_spots))
        }
        
        # sometimes missing, if distances very close
        result[is.na(result)]<-0
        
        # keep only 50 closest features for storage
        if (dim(result)[2] == 1){
          result = result
        } else {
          result = t(apply(result, 1, sort))[, 1:min(ncol(result), 15)] 
        }
        loss_feature = kern(result, h = pull(bandwidths[bandwidths$lake == lake_tmp & bandwidths$feature == 2, "bandwidth"]), option = "loss", ncol = min(ncol(result), 15))
        
      }
      
      
      # locality
      locality_feature = radius/1000
      
      # matrix
      f = cbind(social_feature,
                roughness_feature,
                success_feature,
                loss_feature,
                locality_feature,
                locality_feature*time_since_start[n],
                loss_feature * social_feature,
                success_feature * social_feature)
      
      w =  main_effects$weights + trip_effects$weights[trip_effects$trip == trip_tmp] + lake_effects$weights[lake_effects$lake==lake_tmp] + participant_effects$weights[participant_effects$id == participant_tmp]
      
      if (option == "random"){
        w[1] <- 0
        w[2] <- 0
        w[3] <- 0
        w[4] <- 0
        w[5] <- 0
        w[6] <- 0
        w[7] <- 0
        w[8] <- 0
        
      } else if (option == "brw"){
      
      } else if (option == "non-social"){
        w[1] <- 0
        f[, 1] <- 0
        w[7] <- 0
        w[8] <- 0
      } else if (option == "non-success"){
        w[3] <- 0
        f[,3]<-0
        w[8] <- 0
      } else if (option == "no-interaction"){
        w[7]<-0
        w[8]<-0
      } else if (option == "no-loss"){
        w[4] <- 0
        f[,4]<-0
        w[7] <- 0
      } else if (option == "no-roughness"){
        w[2] <- 0
      }
      
      # compute softmax
      s = softmax(f%*%w)
      
      # select spot acc to softmax
      if (sum(is.na(s))>0){
        ind_fw = sample(1:N_spots, 1, prob = rep(1/N_spots, N_spots))
      } else {
        ind_fw = sample(1:N_spots, 1, prob = s)
      }
      
      # save results
      spot_choices[n,] <- coord_tmp[ind_fw,]
      
      # features
      social_features_choice[n] = social_feature[ind_fw]
      success_features_choice[n] = success_feature[ind_fw]
      loss_features_choice[n] = loss_feature[ind_fw]
      roughness_features_choice[n] = roughness_feature[ind_fw]
      locality_features_choice[n] = locality_feature[ind_fw]
      time_since_last_success_choice[n] = time_since_last_success
      chosen_angle[n] = angle[ind_fw]
      chosen_radius[n] = radius[ind_fw]
      check[n] = sum(is.na(s))
      
    }
    
    # progress
    
    # save results tidy
    if (nrow(spot_choices)>2){
      res <- moveHMM::prepData(as.data.frame(spot_choices), type = 'UTM', coordNames = c('V1', "V2")) 
    } else {
      res = data.frame(x = spot_choices[,1], y = spot_choices[, 2], step = NA, angle = NA)
    }
    
    #hist(res$angle, breaks = 20)
    if ("ID" %in% colnames(res)){
      res = res %>% dplyr::select(-ID)
    } else {
      
    }
    
    res$camera_id = camera_id_tmp
    res$id_index = participant_tmp
    res$day = day_tmp
    res$year = year_tmp
    res$spot_id = 1:n
    res$trip_index = trip_tmp
    res$lake_index = lake_tmp
    res$social_feature = social_features_choice
    res$success_feature = success_features_choice
    res$loss_feature = loss_features_choice
    res$roughness_feature = roughness_features_choice
    res$locality_feature = locality_features_choice
    res$time_since_start = time_since_start
    res$Return = return
    res$time_since_last_success = time_since_last_success_choice
    #res$time_since_last_success = ifelse(res$Return == 1, 0, res$time_since_last_success)
    res$iter = iter
    res$check = check
    res$chosen_angle = chosen_angle
    res$chosen_radius = chosen_radius
  
    Result[[c]] <- res
    c=c+1
    
  }

  }
  
  return(do.call(rbind, Result)) 
}

################################################################################
library(future)
library(furrr)
library(sf)
library(sp)
session_info()
cl <- makeClusterPSOCK(25)
future::plan(future::cluster, workers = cl, gc = T)
options(future.globals.maxSize = 6000*1024^2)
model_sim = future_map(1:50, function(x) simulate_trajectory(x,
                                                            spot_selection_data,
                                                            identifiers,
                                                            competition_area_new,
                                                            depth_profiles,
                                                            gps_data,
                                                            main_effects,
                                                            trip_effects,
                                                            lake_effects,
                                                            participant_effects,
                                                            bandwidths,
                                                            N_spots,
                                                            option = "brw"),
                       .progress = TRUE,
                       .options = furrr_options(seed = TRUE))
options(future.globals.maxSize = 500*1024^2)
predictions = do.call(rbind, model_sim)
parallel::stopCluster(cl)

# random walk
cl <- makeClusterPSOCK(25)
future::plan(cluster, workers = cl, gc = T)
options(future.globals.maxSize = 6000*1024^2)
random_sim = future_map(1:50, function(x) simulate_trajectory(x,
                                                              spot_selection_data,
                                                              identifiers,
                                                              competition_area_new,
                                                              depth_profiles,
                                                              gps_data,
                                                              main_effects,
                                                              trip_effects,
                                                              lake_effects,
                                                              participant_effects,
                                                              bandwidths,
                                                              N_spots,
                                                              option = "random"),
                        .progress = TRUE,
                        .options = furrr_options(seed = TRUE))
options(future.globals.maxSize = 500*1024^2)
parallel::stopCluster(cl)
random = do.call(rbind, random_sim)

# non success walk
cl <- makeClusterPSOCK(25)
future::plan(cluster, workers = cl, gc = T)
options(future.globals.maxSize = 6000*1024^2)
nonsuccess_sim = future_map(1:50, function(x) simulate_trajectory(x,
                                                             spot_selection_data,
                                                             identifiers,
                                                             competition_area_new,
                                                             depth_profiles,
                                                             gps_data,
                                                             main_effects,
                                                             trip_effects,
                                                             lake_effects,
                                                             participant_effects,
                                                             bandwidths,
                                                             N_spots,
                                                             option = "non-success"),
                        .progress = TRUE,
                        .options = furrr_options(seed = TRUE))
options(future.globals.maxSize = 500*1024^2)
parallel::stopCluster(cl)
nonsuccess = do.call(rbind, nonsuccess_sim)

# non social walk
cl <- makeClusterPSOCK(25)
future::plan(cluster, workers = cl, gc = T)
options(future.globals.maxSize = 6000*1024^2)
nonsocial_sim = future_map(1:50, function(x) simulate_trajectory(x,
                                                                 spot_selection_data,
                                                                 identifiers,
                                                                 competition_area_new,
                                                                 depth_profiles,
                                                                 gps_data,
                                                                 main_effects,
                                                                 trip_effects,
                                                                 lake_effects,
                                                                 participant_effects,
                                                                 bandwidths,
                                                                 N_spots,
                                                                 option = "non-social"),
                           .progress = TRUE,
                           .options = furrr_options(seed = TRUE))
options(future.globals.maxSize = 500*1024^2)
parallel::stopCluster(cl)
nonsocial = do.call(rbind, nonsocial_sim)

# non loss walk
cl <- makeClusterPSOCK(25)
future::plan(cluster, workers = cl, gc = T)
options(future.globals.maxSize = 6000*1024^2)
nonloss_sim = future_map(1:50, function(x) simulate_trajectory(x,
                                                                 spot_selection_data,
                                                                 identifiers,
                                                                 competition_area_new,
                                                                 depth_profiles,
                                                                 gps_data,
                                                                 main_effects,
                                                                 trip_effects,
                                                                 lake_effects,
                                                                 participant_effects,
                                                                 bandwidths,
                                                                 N_spots,
                                                                 option = "no-loss"),
                           .progress = TRUE,
                           .options = furrr_options(seed = TRUE))
options(future.globals.maxSize = 500*1024^2)
parallel::stopCluster(cl)
nonloss = do.call(rbind, nonloss_sim)

################################################################################
saveRDS(random, file = "utils/data/processed_data/spot_selection_simulations/sim_random.rds")
saveRDS(predictions, file = "utils/data/processed_data/spot_selection_simulations/sim_predictions.rds")
saveRDS(nonsocial, file = "utils/data/processed_data/spot_selection_simulations/sim_non_social.rds")
saveRDS(nonloss, file = "utils/data/processed_data/spot_selection_simulations/sim_non_loss.rds")
saveRDS(nonsuccess, file = "utils/data/processed_data/spot_selection_simulations/sim_non_success.rds")
################################################################################
}
################################################################################
# END
################################################################################
