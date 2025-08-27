################################################################################
#
# Title: Run patch-leaving models including group level predictors
#
# Authors: Schakowski, A.
#
# Last updated: 27/07/2024 (DD/MM/YYYY)
#
################################################################################
prepare_stan_data_patch_leaving_grp_predictors <- function(){
  
  ################################################################################
  # load angling data
  patch_leaving_data <- read_csv("utils/data/raw_data/patch_leaving_data.csv")
  patch_leaving_data = patch_leaving_data %>% filter(!exclude)
  patch_leaving_data = patch_leaving_data %>% group_by(day, year, camera_id) %>% mutate(n_spots = length(unique(spot_id)))
  first_spots = patch_leaving_data %>% filter(spot_id == 1 & n_spots >1)
  patch_leaving_data = patch_leaving_data %>% filter(spot_id>1)
  
  ##############################################################################
  # join id level predictors to data
  catch_data <- read_csv("utils/data/raw_data/catch_data.csv")
  variables = catch_data %>% 
    group_by(participant_id) %>% 
    slice(1) %>% 
    dplyr::select(sex, year_of_birth)
  variables$age = 2023-variables$year_of_birth
  variables$sex = 2-as.numeric(as.factor(variables$sex))
  patch_leaving_data = left_join(patch_leaving_data, variables)
  
  survey_data = read_csv("utils/data/raw_data/survey_data.csv")
  skill = survey_data %>% dplyr::select(participant_id, angling_skill)
  patch_leaving_data = left_join(patch_leaving_data, skill)
  
  # which ids are missing?
  missing_ids = unique(patch_leaving_data$participant_id[which(is.na(patch_leaving_data$angling_skill))])
  
  # filter
  patch_leaving_data = patch_leaving_data %>% filter(!participant_id %in% missing_ids)
  
  # avg catch by lake
  catches = catch_data %>% group_by(day,year) %>% summarize(avg_catch = mean(catch))
  patch_leaving_data = left_join(patch_leaving_data, catches)
  
  successes = patch_leaving_data %>% 
    group_by(camera_id, spot_id, lake_id) %>% 
    summarize(return = sum(reward)>0) %>% 
    group_by(lake_id) %>% 
    summarize(mean_return = mean(return))
  patch_leaving_data = left_join(patch_leaving_data, successes)
  
  ##############################################################################
  
  # compute unique trip id
  patch_leaving_data$unique_trip_id = as.numeric(as.factor(paste(patch_leaving_data$day, patch_leaving_data$year, patch_leaving_data$camera_id)))
  
  # compute unique participant id
  patch_leaving_data$unique_participant_id = as.numeric(as.factor(paste(patch_leaving_data$participant_id)))
  
  # compute unique lake id
  patch_leaving_data$unique_lake_id = as.numeric(as.factor(paste(patch_leaving_data$day, patch_leaving_data$year)))
  
  # how many unique angling spots
  patch_leaving_data$unique_spot_id = as.numeric(as.factor(paste(patch_leaving_data$unique_trip_id, patch_leaving_data$spot_id)))
  
  # how many unique angling spots
  patch_leaving_data$unique_trip_id = as.numeric(as.factor(paste(patch_leaving_data$unique_trip_id)))
  
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
  
  # compute leaving state
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
  
  # for each trip compute start and end row
  patch_leaving_data = patch_leaving_data %>% 
    arrange(day, year, camera_id, spot_id, bin)
  patch_leaving_data$row_index = 1:nrow(patch_leaving_data)
  
  trip_start = patch_leaving_data %>% group_by(unique_trip_id) %>% summarize(t = min(row_index))
  trip_end = patch_leaving_data %>% group_by(unique_trip_id) %>% summarize(t = max(row_index))
  trip_start = trip_start %>% arrange(unique_trip_id)
  trip_end = trip_end %>% arrange(unique_trip_id)
  trip_start = as.data.frame(trip_start)
  trip_end = as.data.frame(trip_end)
  
  # compute time
  patch_leaving_data$time_sec = patch_leaving_data$time_start + (patch_leaving_data$bin - 1 * 10)
  patch_leaving_data$timestamp_sec = patch_leaving_data$timestamp_start + (patch_leaving_data$bin - 1 * 10)
  
  # lake key
  lake_key = patch_leaving_data %>% group_by(day, year) %>% slice(1) %>% dplyr::select(day, year, lake_index = unique_lake_id)
  
  # lake and id 
  lakes = patch_leaving_data$unique_lake_id
  ids = patch_leaving_data$unique_participant_id
  
  # grp level predictors
  # id level predictors
  id_level_predictors = patch_leaving_data %>% 
    group_by(unique_participant_id) %>% 
    slice(1) %>% 
    arrange(unique_participant_id) %>% 
    dplyr::select(age, sex, angling_skill)
  id_level_predictors = as.data.frame(id_level_predictors)[2:4]
  
  # lake level preditors
  lake_level_predictors = patch_leaving_data %>% 
    group_by(unique_lake_id) %>% 
    slice(1) %>% 
    arrange(unique_lake_id) %>% 
    dplyr::select(avg_catch, mean_return)
  lake_level_predictors = as.data.frame(lake_level_predictors)[2:3]
  
  # make stan data
  stan_data <- list(T = nrow(patch_leaving_data),
                    state = patch_leaving_data$state,
                    catch = patch_leaving_data$catch,
                    time_since_event = patch_leaving_data$time_since_event,
                    cumulative_catch = patch_leaving_data$cumulative_catch,
                    angling_time = patch_leaving_data$angling_time,
                    grainsize = 1,
                    N_trips = length(unique(patch_leaving_data$unique_trip_id)),
                    trips = 1:length(unique(patch_leaving_data$unique_trip_id)),
                    trip_start = trip_start$t,
                    trip_end = trip_end$t,
                    lakes = lakes, 
                    ids = ids,
                    N_ids = length(unique(patch_leaving_data$unique_participant_id)),
                    N_lakes = length(unique(patch_leaving_data$unique_lake_id)),
                    id_level_predictors = id_level_predictors,
                    lake_level_predictors = lake_level_predictors)
  saveRDS(stan_data, file = "utils/data/processed_data/stan_data_patch_leaving_grp.rds")
  
  # make identifiers as for spatial model
  # need day, year, camera id, time bin, ids
  identifiers = patch_leaving_data
  
  # save
  saveRDS(identifiers, file = "utils/data/processed_data/identifiers_patch_leaving_grp.rds")
  
}
################################################################################
# END
################################################################################