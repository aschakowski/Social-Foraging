################################################################################
#
# Title: Prepare Stan Data
#
# Description: Prepare stan data for model fit
#
# Authors: Schakowski, A.
#
# Last updated: 13/02/2024 (DD/MM/YYYY)
#
################################################################################
prepare_stan_data_spot_selection <- function(){
  
  ################################################################################
  # Helper Functions
  ################################################################################
  pad_matrix <- function(matrix, new_rows, new_columns, filler){
    
    if ( new_rows > nrow(matrix) & new_columns > ncol(matrix)){
      
      # add columns
      column_diff = new_columns - ncol(matrix)
      columns = replicate(n = column_diff, expr = rep(filler, nrow(matrix)))
      matrix = cbind(matrix, columns)
      
      # add rows
      row_diff = new_rows - nrow(matrix)
      rows = t(replicate(n = row_diff, expr = rep(filler, ncol(matrix))))
      matrix = rbind(matrix, rows)
      
    } else if (new_rows > nrow(matrix) & new_columns == ncol(matrix)){
      
      # add rows
      row_diff = new_rows - nrow(matrix)
      rows = t(replicate(n = row_diff, expr = rep(filler, ncol(matrix))))
      matrix = rbind(matrix, rows)
      
    } else if (new_rows == nrow(matrix) & new_columns > ncol(matrix) & new_rows > 1){
      
      # add columns
      column_diff = new_columns - ncol(matrix)
      columns = replicate(n = column_diff, expr = rep(filler, nrow(matrix)))
      matrix = cbind(matrix, columns)
      
    } else if (new_rows == 1){
      
      ## add columns
      column_diff = new_columns - ncol(matrix)
      columns = replicate(n = column_diff, expr = rep(filler, 1))
      matrix = c(matrix, columns)
      
    } else {
      # do nothing
      #print("New dimensions are old dimensions.")
    }
    
    return(matrix)
  }
  
  prepare_stan_data <- function(social_distances,
                                roughness_distances,
                                ecological_features,
                                success_distances, 
                                loss_distances, 
                                locality,
                                time_since_start,
                                grainsize, 
                                identifiers, 
                                ncolumns,
                                option = "full",
                                missing_ids = 0){
    
    # subset ids
    if (option == "id_level"){
      sub_ids = identifiers[identifiers$spot_id != 1 & !identifiers$participant_id %in% missing_ids, ]
    } else {
      sub_ids = identifiers[identifiers$spot_id != 1, ]
    }
    
    # lake index
    sub_ids$lake_index = as.numeric(as.factor(sub_ids$lake_id))
    lake = sub_ids$lake_index
    N_lakes = length(unique(sub_ids$lake_index))
    
    # trip index
    sub_ids$trip_index = as.numeric(as.factor(sub_ids$track_id))
    trip = sub_ids$trip_index
    N_trips = length(unique(sub_ids$trip_index))
    
    # participant index
    sub_ids$id_index = as.numeric(as.factor(sub_ids$participant_id))
    id = sub_ids$id_index
    N_ids = length(unique(sub_ids$id_index))
    
    # join back
    identifiers = left_join(identifiers, sub_ids)
    if (option == "id_level"){
      saveRDS(identifiers, file = "utils/data/processed_data/identifiers_model_key_grp_level_predictors.rds")
    } else {
      saveRDS(identifiers, file = "utils/data/processed_data/identifiers_model_key.rds")
    }
    
    # id level predictors
    id_level_predictors = sub_ids %>% 
      group_by(id_index) %>% 
      slice(1) %>% 
      arrange(id_index) %>% 
      dplyr::select(age, sex, angling_skill)
    id_level_predictors = as.data.frame(id_level_predictors)[2:4]
    
    # lake level preditors
    lake_level_predictors = sub_ids %>% 
      group_by(lake_index) %>% 
      slice(1) %>% 
      arrange(lake_index) %>% 
      dplyr::select(catch, mean_return)
    lake_level_predictors = as.data.frame(lake_level_predictors)[2:3]
    
    # keep key to link new ids to old ones
    key_participant_id = data.frame(participant_id = sub_ids$participant_id,
                                    model_participant_id = id)
    key_track_id = data.frame(track_id = sub_ids$track_id,
                              model_track_id = trip)
    key_lake_id = data.frame(lake_id = sub_ids$lake_id,
                             model_lake_id = lake)
    keys = list(key_participant_id,
                key_track_id,
                key_lake_id)
    
    if (option == "id_level"){
      saveRDS(key, file = "utils/data/processed_data/model_key_grp_level_predictors.rds")
    } else {
      saveRDS(key, file = "utils/data/processed_data/model_key.rds")
    }
    
    
    # compute success columns
    success_columns = vector(mode = "numeric", length = length(success_distances))
    for (i in 1:length(success_distances)){
      
      if (is.null(dim(success_distances[[i]]))){
        success_columns[i] <- 1
      } else {
        success_columns[i] <- ncol(success_distances[[i]])
        success_columns[i] <- min(ncolumns, success_columns[[i]]) 
      }
    }
    
    
    # compute loss columns
    loss_columns = vector(mode = "numeric", length = length(loss_distances))
    for (i in 1:length(loss_distances)){
      
      if (is.null(dim(loss_distances[[i]]))){
        loss_columns[i] <- 1
      } else {
        loss_columns[i] <- ncol(loss_distances[[i]])
        loss_columns[i] <- min(ncolumns, loss_columns[[i]]) 
      }
    }
    
    # select closest ncol spots
    for (i in 1:length(social_distances)){
      
      # social
      social_distances[[i]][is.na(social_distances[[i]])] <- 0.0001
      social_distances[[i]] = t(apply(social_distances[[i]], 1, sort))
      social_distances[[i]] <- social_distances[[i]][, 1:ncolumns]
      
      # roughness
      roughness_distances[[i]][is.na(roughness_distances[[i]])] <- 0.0001
      roughness_distances[[i]] = t(apply(roughness_distances[[i]], 1, sort))
      roughness_distances[[i]] <- as.matrix(roughness_distances[[i]])[, 1:ncolumns]
      
      
      # select last ncol spots
      if (success_columns[[i]] == 1){
        success_distances[[i]] = success_distances[[i]]
      } else {
        success_distances[[i]] = success_distances[[i]][, 1:success_columns[i]]
      }
      
      if (loss_columns[[i]] == 1){
        loss_distances[[i]] = loss_distances[[i]]
      } else {
        loss_distances[[i]] = loss_distances[[i]][, 1:loss_columns[i]] 
      }
      
    }
    
    # pad matrix to max columns
    for (i in 1:length(social_distances)){
      
      social_distances[[i]] <- pad_matrix(social_distances[[i]],
                                          new_rows = nrow(social_distances[[i]]),
                                          new_columns = ncolumns,
                                          filler = -99)
      
      success_distances[[i]] <- pad_matrix(as.matrix(success_distances[[i]]),
                                           new_rows = nrow(as.matrix(success_distances[[i]])),
                                           new_columns = ncolumns,
                                           filler = -99)
      
      roughness_distances[[i]] <- pad_matrix(roughness_distances[[i]],
                                             new_rows = nrow(roughness_distances[[i]]),
                                             new_columns = ncolumns,
                                             filler = -99)
      
      loss_distances[[i]] <- pad_matrix(as.matrix(loss_distances[[i]]),
                                        new_rows = nrow(as.matrix(loss_distances[[i]])),
                                        new_columns = ncolumns,
                                        filler = -99)
      
    }
    
    # convert list of matrices to stan conform array
    social_feature <- abind(social_distances, along=3)
    social_feature <- aperm(social_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
    
    success_feature <- abind(success_distances, along=3)
    success_feature <- aperm(success_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
    
    loss_feature <- abind(loss_distances, along=3)
    loss_feature <- aperm(loss_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
    
    roughness_feature <- abind(roughness_distances, along=3)
    roughness_feature <- aperm(roughness_feature, c(3,1,2)) #dim 1: time, dim2: alternative spots, dim3: distances to 55 players
    
    ecological_features <- abind(ecological_features, along = 3)
    ecological_features <- aperm(ecological_features, c(3, 1, 2))
    
    # locality matrix: columns: all spot choices, rows: alternative spots
    locality_feature = do.call(cbind, locality)
    
    # number of alternative spots
    N_spots = nrow(social_distances[[1]])
    
    # choice vector
    choices = 1:length(social_distances)
    
    # grainsize for within chain parallelization
    grainsize = grainsize
    
    # return padded distance matrix and columns
    if (option == "id_level"){
      return(list(Distances_social = social_feature,
                  Distances_roughness = roughness_feature, 
                  Ecological_features = ecological_features,
                  Distances_successes = success_feature,
                  Distances_losses = loss_feature,
                  Locality = locality_feature,
                  Time_since_start = as.numeric(time_since_start)/max(as.numeric(time_since_start)),
                  success_columns = success_columns, 
                  loss_columns = loss_columns,
                  T = length(social_distances),
                  N_spots = N_spots, 
                  choices = choices, 
                  grainsize = grainsize,
                  columns = ncolumns,
                  lake = lake,
                  N_lakes = N_lakes,
                  trip = trip,
                  N_trips = N_trips, 
                  id = id,
                  N_ids = N_ids,
                  id_level_predictors = id_level_predictors,
                  lake_level_predictors = lake_level_predictors)) 
    } else {
      return(list(Distances_social = social_feature,
                  Distances_roughness = roughness_feature, 
                  Ecological_features = ecological_features,
                  Distances_successes = success_feature,
                  Distances_losses = loss_feature,
                  Locality = locality_feature,
                  Time_since_start = as.numeric(time_since_start)/max(as.numeric(time_since_start)),
                  success_columns = success_columns, 
                  loss_columns = loss_columns,
                  T = length(social_distances),
                  N_spots = N_spots, 
                  choices = choices, 
                  grainsize = grainsize,
                  columns = ncolumns,
                  lake = lake,
                  N_lakes = N_lakes,
                  trip = trip,
                  N_trips = N_trips, 
                  id = id,
                  N_ids = N_ids))
    }
    
  }
  ################################################################################
  # spot ids
  identifiers = readRDS(file = "utils/data/processed_data/identifiers_spot_selection.rds")
  
  # join id level predictors to identifiers
  catch_data <- read_csv("utils/data/raw_data/catch_data.csv")
  variables = catch_data %>% 
    group_by(participant_id) %>% 
    slice(1) %>% 
    dplyr::select(sex, year_of_birth)
  variables$age = 2023-variables$year_of_birth
  variables$sex = 2-as.numeric(as.factor(variables$sex))
  identifiers = left_join(identifiers, variables)
  
  survey_data = read_csv("utils/data/raw_data/survey_data.csv")
  skill = survey_data %>% dplyr::select(participant_id, angling_skill)
  identifiers = left_join(identifiers, skill)
  
  # join lake level predictors
  lake_variables = catch_data %>% group_by(day, year) %>% summarize(catch = mean(catch))
  ret = spot_selection_data %>% group_by(day, year) %>% summarize(mean_return = mean(Return))
  lake_variables = left_join(lake_variables, ret)
  identifiers = left_join(identifiers, lake_variables)
  
  # load raw features
  social_distances = readRDS(file = "utils/data/processed_data/social_distance_matrix.rds")
  roughness_distances = readRDS(file = "utils/data/processed_data/roughness_distance_matrix.rds")
  ecological_features = readRDS(file = "utils/data/processed_data/alternative_roughness_distance_matrix.rds")
  success_distances = readRDS(file = "utils/data/processed_data/success_distance_matrix.rds")
  loss_distances = readRDS(file = "utils/data/processed_data/loss_distance_matrix.rds")
  locality = readRDS(file = "utils/data/processed_data/locality_matrix.rds")
  
  # time since start as feature
  time_since_start = identifiers$time_since_start
  
  # drop first choices
  social_distances = social_distances[(identifiers$spot_id!= 1)]
  roughness_distances = roughness_distances[(identifiers$spot_id!= 1)]
  ecological_features = ecological_features[(identifiers$spot_id!= 1)]
  success_distances = success_distances[(identifiers$spot_id!= 1)]
  loss_distances = loss_distances[(identifiers$spot_id!= 1)]
  locality = locality[(identifiers$spot_id!= 1)]
  time_since_start = time_since_start[identifiers$spot_id!=1]
  
  # compute feature matrix
  data_list = prepare_stan_data(social_distances = social_distances,
                                roughness_distances = roughness_distances, 
                                ecological_features = ecological_features,
                                success_distances = success_distances,
                                loss_distances = loss_distances,
                                locality = locality, 
                                time_since_start = time_since_start,
                                grainsize = 1, identifiers = identifiers, ncolumns = 15)
  saveRDS(data_list, file = "utils/data/processed_data/stan_data_spot_selection.rds")
  
  ################################################################################
  # REVISION: Create stan data for id level predictions (remove missings first)
  ################################################################################
  # which ids are missing?
  missing_ids = unique(identifiers$participant_id[which(is.na(identifiers$angling_skill))])
  
  # compute feature matrix
  data_list = prepare_stan_data(social_distances = social_distances,
                                roughness_distances = roughness_distances, 
                                ecological_features = ecological_features,
                                success_distances = success_distances,
                                loss_distances = loss_distances,
                                locality = locality, 
                                time_since_start = time_since_start,
                                grainsize = 1, identifiers = identifiers, ncolumns = 15,
                                option = "id_level",
                                missing_ids = missing_ids)
  saveRDS(data_list, file = "utils/data/processed_data/stan_data_spot_selection_grp_level_predictors.rds")
  
}

################################################################################
# END
################################################################################