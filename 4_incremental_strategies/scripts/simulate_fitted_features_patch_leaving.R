################################################################################
#
# Title: Simulate fitted features patch leaving
#
# Description: 
#
# Authors: Schakowski, A.
#
# Last updated: 02/04/2024 (DD/MM/YYYY)
#
################################################################################
simulate_fitted_features_patch_leaving<-function(){
  
  # load model
  patch_leaving_model_fit = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/e_patch_leaving_model_bandwidth_fit.rds")
  stan_data = readRDS(file = "utils/data/processed_data/stan_data_patch_leaving.rds")
  
  # avg bandwidth parameters for each lake and each feature
  main_lambda = patch_leaving_model_fit %>% 
    spread_draws(lambda) %>% 
    summarize(log_lambda = mean(lambda))
  lake_lambdas = patch_leaving_model_fit %>% 
    spread_draws(v_lake[lake,feature]) %>% 
    filter(feature %in% c(12)) %>% 
    mutate(feature = feature - 11) %>% 
    group_by(feature, lake) %>% 
    summarize(log_lambda_offset = mean(v_lake))
  main_lambda$feature = 1
  
  # join to one dataframe
  lambdas = left_join(main_lambda, lake_lambdas)
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
    choice = vector()
    spot_id = vector()
    chosen = vector()
    
    # success feature
    success_feature = stan_data$success_feature[t];  
    
    # loss feature
    loss_feature = stan_data$loss_feature[t]
    
    # roughness feature
    roughness_feature = stan_data$roughness_feature[t]
    
    # social feature
    social_feature = sum(exp( - (stan_data$social_matrix[t,] ^ 2) / (2*square(lambdas$lake_lambda[lambdas$lake == stan_data$lake[t] & lambdas$feature == 1]))));
    
    # identifiers
    lake = stan_data$lake[t]
    id = stan_data$id[t]
    
    # make tidy dataframe
    fitted_features_df = data.frame(success_feature = success_feature,
                                    loss_feature = loss_feature,
                                    social_feature = social_feature,
                                    roughness_feature = roughness_feature,
                                    lake = lake,
                                    id = id,
                                    t = t)
    # save to list
    fitted_features[[t]]<-fitted_features_df
  }
  
  # make df from list
  fitted_features = fitted_features %>% dplyr::bind_rows()
  fitted_features$unique_trip = as.numeric(as.factor(paste(fitted_features$lake, fitted_features$id)))
  # load ids to match back to dataset
  saveRDS(fitted_features, file = "utils/data/processed_data/fitted_features_patch_leaving.rds")
  
}

################################################################################
# END
################################################################################