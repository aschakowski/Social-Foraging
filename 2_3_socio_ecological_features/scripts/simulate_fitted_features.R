#################################################################################
#
# Title: Compute fitted features
#
# Description: Simulate spatial features from posterior model estimates.
#
# Authors: Schakowski, A.
#
# Last updated: 13/02/2024 (DD/MM/YYYY)
#
################################################################################
simulate_fitted_features<-function(){
  
  # load model
  spot_selection_fit = readRDS(file = "utils/data/processed_data/model_fit/spot_selection_models/f_conditional.rds")
  
  # load data
  stan_data = readRDS(file = "utils/data/processed_data/stan_data_spot_selection.rds")
  
  # avg bandwidth parameters for each lake and each feature
  main_lambdas = spot_selection_fit %>% 
    spread_draws(lambdas[feature]) %>% 
    group_by(feature) %>% summarize(log_lambda = mean(lambdas))
  lake_lambdas = spot_selection_fit %>% 
    spread_draws(v_lake[lake,feature]) %>% 
    filter(feature %in% c(9, 10)) %>% 
    mutate(feature = feature - 8) %>% 
    group_by(feature, lake) %>% 
    summarize(log_lambda_offset = mean(v_lake))
  
  # join to one dataframe
  lambdas = left_join(main_lambdas, lake_lambdas)
  lambdas$lake_lambda = exp(lambdas$log_lambda_offset + lambdas$log_lambda)
  
  # loop through all choices and compute feature values
  fitted_features = list()
  for (t in 1:stan_data$T){#stan_data$T){
    
    # preallocate results
    success_feature = vector()
    loss_feature = vector()
    social_feature = vector()
    roughness_feature = vector()
    lake = vector()
    id = vector()
    trip = vector()
    choice = vector()
    spot_id = vector()
    chosen = vector()
    
    # success feature
    if (stan_data$Distances_successes[t,1,1] < 0){
      success_feature = rep(0, stan_data$N_spots)
    } else {
      if (stan_data$success_columns[t] == 1){
        success_feature = exp( - (stan_data$Distances_successes[t,, 1:stan_data$success_columns[t]] ^ 2) /  (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 2])))
      } else {
        success_feature = (exp( - (stan_data$Distances_successes[t,, 1:stan_data$success_columns[t]] ^ 2) /  (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 2]))) %*% rep(1.0, stan_data$success_columns[t]));  
      }
    }
    
    # loss feature
    if (stan_data$Distances_losses[t,1,1] < 0){
      loss_feature = rep(0, stan_data$N_spots);
    } else {
      if (stan_data$loss_columns[t] == 1){
        loss_feature = exp( - (stan_data$Distances_losses[t,, 1:stan_data$loss_columns[t]] ^ 2) /  (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 2])));  
      } else {
        loss_feature = (exp( - (stan_data$Distances_losses[t,, 1:stan_data$loss_columns[t]] ^ 2) /  (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 2]))) %*% rep(1.0, stan_data$loss_columns[t])); 
      }
    }
    
    # roughness feature
    roughness_feature = stan_data$Distances_roughness[t,, 1] / 1000;
    
    # social feature
    social_feature = rowSums(exp( - (stan_data$Distances_social[t,,] ^ 2) / (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 1]))) * (1.0 / (1 + exp(-20 * (stan_data$Distances_social[t,,] - 5)))));
    
    # locality feature
    locality_feature = stan_data$Locality[, t]/1000;
    
    # identifiers
    lake = stan_data$lake[t]
    id = stan_data$id[t]
    trip = stan_data$trip[t]
    choice = stan_data$choices[t]
    chosen = c(1, rep(0, stan_data$N_spots-1))
    
    # compute spot_id
    if(t == 1){
      spot_id = rep(2, stan_data$N_spots)
    } else if (stan_data$trip[t] == stan_data$trip[t-1]){
      spot_id = fitted_features[[t-1]]$spot_id + 1
    } else if (stan_data$trip[t] != stan_data$trip[t-1]){
      spot_id = rep(2, stan_data$N_spots)
    }
    
    
    # make tidy dataframe
    fitted_features_df = data.frame(success_feature = success_feature,
                                    loss_feature = loss_feature,
                                    social_feature = social_feature,
                                    roughness_feature = roughness_feature,
                                    locality_feature = locality_feature,
                                    lake = lake,
                                    id = id,
                                    trip = trip,
                                    choice = choice,
                                    spot_id = spot_id,
                                    chosen = chosen)
    
    # save to list
    fitted_features[[t]]<-fitted_features_df
  }
  
  # make df from list
  fitted_features = do.call(rbind, fitted_features)
  
  # load ids to match back to dataset
  identifiers = readRDS(file = "utils/data/processed_data/identifiers_spot_selection.rds")
  identifiers = identifiers[identifiers$spot_id != 1,]
  identifiers$lake = as.numeric(as.factor(identifiers$lake_id))
  identifiers$trip = as.numeric(as.factor(identifiers$track_id))
  identifiers$id = as.numeric(as.factor(identifiers$participant_id))
  identifiers$spot_id = as.numeric(identifiers$spot_id)
  fitted_features = left_join(fitted_features, identifiers)
  saveRDS(fitted_features, file = "utils/data/processed_data/fitted_features.rds")
  
}

################################################################################
# END
################################################################################