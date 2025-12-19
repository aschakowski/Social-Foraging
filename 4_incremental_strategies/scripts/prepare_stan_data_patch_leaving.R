################################################################################
#
# Title: Run patch-leaving models
#
# Authors: Schakowski, A.
#
# Last updated: 27/07/2024 (DD/MM/YYYY)
#
################################################################################
prepare_stan_data_patch_leaving <- function(){

################################################################################
# load angling data
patch_leaving_data <- read_csv("utils/data/raw_data/patch_leaving_data.csv")
patch_leaving_data = patch_leaving_data %>% filter(!exclude)
patch_leaving_data = patch_leaving_data %>% group_by(day, year, camera_id) %>% mutate(n_spots = length(unique(spot_id)))
first_spots = patch_leaving_data %>% filter(spot_id == 1 & n_spots >1)
patch_leaving_data = patch_leaving_data %>% filter(spot_id>1)

# compute unique trip id
patch_leaving_data$unique_trip_id = as.numeric(as.factor(paste(patch_leaving_data$day, patch_leaving_data$year, patch_leaving_data$camera_id)))

# compute unique participant id
patch_leaving_data$unique_participant_id = as.numeric(as.factor(paste(patch_leaving_data$participant_id)))

# compute unique lake id
patch_leaving_data$unique_lake_id = as.numeric(as.factor(paste(patch_leaving_data$day, patch_leaving_data$year)))

# how many unique angling spots (spot ids not in numerical order, but variable not used in models.)
patch_leaving_data$unique_spot_id = as.numeric(as.factor(paste(patch_leaving_data$unique_trip_id, patch_leaving_data$spot_id)))

# compute cumulative catch at spots
patch_leaving_data = patch_leaving_data %>% 
  group_by(unique_spot_id) %>% 
  mutate(cumulative_catch = cumsum(reward))

# compute time since event
patch_leaving_data = patch_leaving_data %>% 
  arrange(day, year, camera_id, spot_id, bin)

time_since_event = vector()
time_since_event[1] = 0
for (i in 2:nrow(patch_leaving_data)){
  
  if (patch_leaving_data$reward[i] == 1 || patch_leaving_data$unique_spot_id[i] != patch_leaving_data$unique_spot_id[i-1]){
    
    time_since_event[i] = 0
    
  } else {
    time_since_event[i] = 1 + time_since_event[i-1]
  }
  
}

patch_leaving_data$time_since_event = time_since_event

# compute leaving state (if last spot and not left at end, censored event)
patch_leaving_data = 
  patch_leaving_data %>% 
  group_by(unique_spot_id) %>% 
  mutate(state = ifelse(bin == max(bin), 1, 0)) %>% 
  mutate(state = ifelse(last_spot == 1, 0, state))

# compute angling time
patch_leaving_data = 
  patch_leaving_data %>% 
  group_by(unique_spot_id) %>% 
  mutate(angling_time = 1:n())

# compute catch
patch_leaving_data$catch = patch_leaving_data$reward

# for each trip compute start and end row (first arrange in correct order)
patch_leaving_data = patch_leaving_data %>% 
  arrange(day, year, camera_id, spot_id, bin)
patch_leaving_data$row_index = 1:nrow(patch_leaving_data)

trip_start = patch_leaving_data %>% group_by(unique_trip_id) %>% summarize(t = min(row_index))
trip_end = patch_leaving_data %>% group_by(unique_trip_id) %>% summarize(t = max(row_index))
trip_start = trip_start %>% arrange(unique_trip_id)
trip_end = trip_end %>% arrange(unique_trip_id)
trip_start = as.data.frame(trip_start)
trip_end = as.data.frame(trip_end)

# compute initial values
initial_values = first_spots %>% 
  group_by(camera_id, day, year) %>%
  mutate(cumulative_catch = cumsum(reward)) %>% 
  filter(bin == max(bin)) %>% 
  summarize(initial_value = ifelse(cumulative_catch == 0, bin, bin / (cumulative_catch+1)))
hist(initial_values$initial_value)

# left join unique trip id to initial values
initial_values = left_join(initial_values, patch_leaving_data %>% group_by(unique_trip_id) %>% slice(1) %>%  dplyr::select(camera_id, day, year, unique_trip_id))
initial_values = initial_values %>% arrange(unique_trip_id)
initial_values = as.data.frame(initial_values %>% ungroup() %>%  dplyr::select(unique_trip_id, initial_value))

# add spatial features (social feature will be recomputed to obtain bandwidth within model)
spatial_features = readRDS(file = "utils/data/processed_data/fitted_features.rds")
spatial_features = spatial_features %>% group_by(day, year, camera_id, spot_id) %>% slice(1) %>% dplyr::select(success_feature, loss_feature, social_feature, roughness_feature, day, year, camera_id, spot_id)
spatial_features$day = as.numeric(spatial_features$day)
spatial_features$year = as.numeric(spatial_features$year)
patch_leaving_data = left_join(patch_leaving_data, spatial_features)

# compute contiuous social feature across whole angling time
# load gps data
gps_data <- fread("utils/data/raw_data/gps_data.csv")

# fast distance computation
sourceCpp("utils/library/fast_dist.cpp")

# bandwidth estimates (here, social feature is computed based on bandwidth estimates from spot selection model, yields same results as 
# reported for additional bandwidth estimate reported in the paper)
spot_selection_fit = readRDS(file = "utils/data/processed_data/model_fit/spot_selection_models/f_conditional.rds")
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

# project coordinates
coords_xy = LongLatToUTM(gps_data$position_long, gps_data$position_lat, zone = 35)
gps_data$position_x = coords_xy[, 2]
gps_data$position_y = coords_xy[, 3]

# compute time (reconvert time bin to time in seconds at start of time bin)
patch_leaving_data$time_sec = patch_leaving_data$time_start + ((patch_leaving_data$bin - 1) * 10)
patch_leaving_data$timestamp_sec = patch_leaving_data$timestamp_start + ((patch_leaving_data$bin - 1) * 10)

# lake key
lake_key = patch_leaving_data %>% group_by(day, year) %>% slice(1) %>% dplyr::select(day, year, lake_index = unique_lake_id)

# write as function
compute_social_feature <- function(timestamp_sec,
                                   unique_lake_id,
                                   camera_id,
                                   lake_key,
                                   gps_data,
                                   bandwidths, 
                                   option){
  
  # define cpp function for cluster export
  sourceCpp("utils/library/fast_dist.cpp")
  
  # get day and year
  tmp = lake_key[lake_key$lake_index == unique_lake_id, ]
  day_tmp = tmp$day
  year_tmp = tmp$year
  
  # tmp id
  id_tmp = camera_id
  
  target_coordinates = gps_data %>% 
    filter(camera_id == id_tmp & year == year_tmp & day == day_tmp & timestamp == timestamp_sec) %>% 
    dplyr::select(position_x, position_y)
  
  # sometimes missing, if at end of competition (take location 5 minutes prior)
  if (nrow(target_coordinates) == 0){
    target_coordinates = gps_data %>% 
      filter(camera_id == id_tmp & year == year_tmp & day == day_tmp & timestamp == timestamp_sec - 5*60) %>% 
      dplyr::select(position_x, position_y)
    
  }
  
  # extract coordinates social
  social_coordinates = gps_data %>% 
    filter(camera_id != id_tmp & year == year_tmp & day == day_tmp & timestamp == timestamp_sec) %>% 
    dplyr::select(position_x, position_y)
  
  # sometimes missing if at end of competition, take 5 minutes prior (rare, 0.01% of points)
  if(nrow(social_coordinates) < 15){
    social_coordinates = gps_data %>% 
      filter(camera_id != id_tmp & year == year_tmp & day == day_tmp & timestamp == timestamp_sec - 5*60) %>% 
      dplyr::select(position_x, position_y)
  }
  
  # compute distances
  distances = fastPdist2(as.matrix(target_coordinates), as.matrix(social_coordinates))
  
  if (option == "selection"){
    # extract bandwidth 
    bandwidth_tmp = bandwidths$bandwidth[bandwidths$lake == unique_lake_id & bandwidths$feature == 1]
    
    # compute kde
    return(kern(d = sort(distances)[1:15], h = bandwidth_tmp, option = "success", ncol = 15))
    
  } else {
    return(sort(distances)[1:15])
  }
  
}

# check if continous social feature already computed, if not recompute (not used in paper)
if(!file.exists("utils/data/processed_data/social_feature_continuous.rds")){

# parallelize
cl <- makeClusterPSOCK(40)
future::plan(cluster, workers = cl, gc = T)
options(future.globals.maxSize= 6000*1024^2)

# compute distances
start_time= Sys.time()
social_feature_angling = future_pmap(list(patch_leaving_data$timestamp_sec, patch_leaving_data$unique_lake_id, patch_leaving_data$camera_id),function(x, y, z) compute_social_feature(x, y, z, lake_key, gps_data, bandwidths, option = "selection"), .progress = TRUE)
end_time = Sys.time()

# stop cluster
parallel::stopCluster(cl)

# combine
social_feature_angling_vector = unlist(social_feature_angling)

# save
saveRDS(social_feature_angling_vector, file = "utils/data/processed_data/social_feature_continuous.rds")
}

################################################################################
# REVISION: recompute social features
# check if continous social feature already computed, if not recompute (used in paper)
if(!file.exists("utils/data/processed_data/social_distances_continuous.rds")){
  
  # parallelize
  cl <- makeClusterPSOCK(40)
  future::plan(cluster, workers = cl, gc = T)
  options(future.globals.maxSize= 6000*1024^2)
  
  # compute distances
  start_time= Sys.time()
  social_feature_angling = future_pmap(list(patch_leaving_data$timestamp_sec, patch_leaving_data$unique_lake_id, patch_leaving_data$camera_id),function(x, y, z) compute_social_feature(x, y, z, lake_key, gps_data, bandwidths, option = "leaving"), .progress = TRUE)
  end_time = Sys.time()
  
  # stop cluster
  parallel::stopCluster(cl)
  
  # combine
  social_feature_angling_matrix = do.call(rbind, social_feature_angling)
  
  # save
  saveRDS(social_feature_angling_matrix, file = "utils/data/processed_data/social_distances_continuous.rds")
}

################################################################################

# read continuous social feature
social_feature_angling_vector = readRDS(file = "utils/data/processed_data/social_feature_continuous.rds")

# read continous social distances matrix
social_feature_angling_matrix = readRDS("utils/data/processed_data/social_distances_continuous.rds")

# add to data
patch_leaving_data$social_feature_cont = social_feature_angling_vector

# lake and id 
lakes = patch_leaving_data$unique_lake_id
ids = patch_leaving_data$unique_participant_id

# test
plot(patch_leaving_data$social_feature_cont[500:1000], patch_leaving_data$social_feature[500:1000])

# check that spots within trips are in right chronological order. They are, if function is monotonous
tid = sample(1:length(unique(patch_leaving_data$unique_trip_id)),1)
patch_leaving_data %>% 
  filter(unique_trip_id == tid) %>%
  ungroup() %>% 
  mutate(r = 1:n()) %>% 
  ggplot(aes(x = r, y = sequence)) + 
  geom_point()

# check that trip indices are correct (trip id should change)
tid = sample(1:length(unique(patch_leaving_data$unique_trip_id)),1)
patch_leaving_data[(trip_start[tid,"t"]-1):(trip_end[tid,"t"]+1),] %>% 
  ungroup() %>% 
  mutate(r = 1:n()) %>% 
  ggplot(aes(x = r, y = unique_trip_id)) + 
  geom_point()

# make stan data
stan_data <- list(T = nrow(patch_leaving_data),
                  state = patch_leaving_data$state,
                  catch = patch_leaving_data$catch,
                  time_since_event = patch_leaving_data$time_since_event,
                  cumulative_catch = patch_leaving_data$cumulative_catch,
                  angling_time = patch_leaving_data$angling_time,
                  success_feature = patch_leaving_data$success_feature,
                  loss_feature = patch_leaving_data$loss_feature,
                  roughness_feature = patch_leaving_data$roughness_feature,
                  social_feature = patch_leaving_data$social_feature_cont,
                  social_matrix = social_feature_angling_matrix,
                  grainsize = 1,
                  N_trips = nrow(initial_values),
                  trips = 1:nrow(initial_values),
                  initial_values = initial_values$initial_value,
                  trip_start = trip_start$t,
                  trip_end = trip_end$t,
                  lakes = lakes, 
                  ids = ids,
                  N_ids = length(unique(ids)),
                  N_lakes = length(unique(lakes)))
saveRDS(stan_data, file = "utils/data/processed_data/stan_data_patch_leaving.rds")

# make identifiers as for spatial model
# need day, year, camera id, time bin, ids
identifiers = data.frame(camera_id = patch_leaving_data$camera_id,
                         day = patch_leaving_data$day,
                         spot_id = patch_leaving_data$spot_id,
                         year = patch_leaving_data$year,
                         bin = patch_leaving_data$bin,
                         participant_id = patch_leaving_data$participant_id,
                         track_id = patch_leaving_data$track_id,
                         lake_id = patch_leaving_data$lake_id,
                         model_participant_id = ids,
                         model_trip_id = patch_leaving_data$unique_trip_id,
                         model_lake_id = lakes)

# save
saveRDS(identifiers, file = "utils/data/processed_data/identifiers_patch_leaving.rds")

}
################################################################################
# END
################################################################################