################################################################################
#
# Title: Run and compare depletion models
#
# Authors: Schakowski, A.
#
# Last updated: 27/07/2024 (DD/MM/YYYY)
#
################################################################################
run_model_comparison_patch_leaving <- function(){

################################################################################
# helper function

# fit models
model_fit <- function(model_location, save_location, data_list, warmup, iterations, chains){
  # load model
  stan_model = cmdstan_model(model_location, 
                             cpp_options = list(stan_threads = TRUE))
  
  # fit models
  stan_model_fit <- stan_model$sample(
    data = data_list,
    chains = chains,
    parallel_chains = chains,
    refresh = 5,
    iter_sampling = iterations,
    iter_warmup = warmup,
    threads_per_chain = round(40/chains)
  )
  
  # save
  stan_model_fit$save_object(file = save_location)
  
  # empty return
  return()
  
}

# extract log lik
compute_waic <- function(fit_location, save_location){
  
  # load model fit
  model_fit = readRDS(file = fit_location)
  
  ll = model_fit$draws(variables = "log_lik")
  waic = waic(ll)
  
  # save loo
  saveRDS(waic, file = save_location)
  
  # keep workspace clear
  gc()
  
  # empty return
  return()
  
}

# extract log lik per lake
compute_waic_lakes <- function(fit_location, save_location, stan_data){
  
  # load model fit
  model_fit = readRDS(file = fit_location)
  
  # extract log likelihood matrix
  ll = model_fit$draws(variables = "log_lik")
  
  # repeat for all lakes
  waic_list = list()
  for (i in 1:10){
    # subset lake indices
    indices = which(stan_data$lake == i)
    ll_sub = ll[,,indices]
    waic_list[[i]] = waic(ll_sub)
  }
  
  # save loo
  saveRDS(waic_list, file = save_location)
  
  # keep workspace clear
  gc()
  
  # empty return
  return()
  
}

################################################################################
# load stan data
stan_data = readRDS(file = "utils/data/processed_data/stan_data_patch_leaving.rds")

# step 1: Fit different strategies
# fit models
model_locations = c("utils/stan/patch_leaving_models/a_patch_leaving_model_baseline.stan",
                    "utils/stan/patch_leaving_models/b_patch_leaving_model_time_wo_catch.stan",#gut
                    "utils/stan/patch_leaving_models/b_patch_leaving_model_incremental.stan",
                    "utils/stan/patch_leaving_models/b_patch_leaving_model_fixed_n.stan",
                    "utils/stan/patch_leaving_models/b_patch_leaving_model_fixed_time.stan")
save_locations = c("utils/data/processed_data/model_fit/patch_leaving_models/a_patch_leaving_model_baseline_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_time_wo_catch_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_incremental_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_fixed_n_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_fixed_time_fit.rds")

for (i in 4:length(model_locations)){
  model_fit(model_location = model_locations[i],
            save_location = save_locations[i],
            data_list = stan_data,
            warmup = 1000,
            iterations = 2000,
            chains = 4)
  gc()
  gc()
  gc()
}

# run model comparison
fit_locations = save_locations
save_locations = c("utils/data/processed_data/model_fit/patch_leaving_models/comparison/a_loo_patch_leaving_model_baseline.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_patch_leaving_model_time_wo_catch.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_patch_leaving_model_incremental.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_patch_leaving_model_fixed_n.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_patch_leaving_model_fixed_time.rds")

for (i in 1:length(fit_locations)){
  compute_waic(fit_location = fit_locations[i],
               save_location = save_locations[i])
  gc()
  gc()
}

# conduct model comparison
loo_list = list(readRDS(save_locations[1]), 
                readRDS(save_locations[2]), 
                readRDS(save_locations[3]), 
                readRDS(save_locations[4]), 
                readRDS(save_locations[5]))
comparison = loo_compare(loo_list)

# save results
saveRDS(comparison, file = "utils/data/processed_data/model_fit/patch_leaving_models/comparison/comparison_basic_strategies.rds")
############################################
# fit models
model_locations = c("utils/stan/patch_leaving_models/c_patch_leaving_model_patch_discovery.stan",
                    "utils/stan/patch_leaving_models/d_patch_leaving_model_globallocal_updating.stan",
                    "utils/stan/patch_leaving_models/e_patch_leaving_model_bandwidth.stan")
save_locations = c("utils/data/processed_data/model_fit/patch_leaving_models/c_patch_leaving_model_patch_discovery_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/d_patch_leaving_model_globallocal_updating_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/e_patch_leaving_model_bandwidth_fit.rds")

for (i in 1:length(model_locations)){
  model_fit(model_location = model_locations[i],
            save_location = save_locations[i],
            data_list = stan_data,
            warmup = 1000,
            iterations = 2000,
            chains = 4)
  gc()
  gc()
  gc()
}

# run model comparison
fit_locations = save_locations
save_locations = c("utils/data/processed_data/model_fit/patch_leaving_models/comparison/c_loo_patch_discovery.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/comparison/d_loo_updating.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/comparison/e_loo_bandwidth.rds")

for (i in 1:length(fit_locations)){
  compute_waic(fit_location = fit_locations[i],
               save_location = save_locations[i])
  gc()
  gc()
  gc()
}

# conduct model comparison
loo_list = list(readRDS("utils/data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_patch_leaving_model_time_wo_catch.rds"),
                readRDS(save_locations[1]), 
                readRDS(save_locations[2]), 
                readRDS(save_locations[3]))
comparison = loo_compare(loo_list)

# save results
saveRDS(comparison, file = "utils/data/processed_data/model_fit/patch_leaving_models/comparison/comparison_extended_strategies.rds")

################################################################################
# REVISION: Comparisons by lake
fit_locations = c("utils/data/processed_data/model_fit/patch_leaving_models/a_patch_leaving_model_baseline_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_time_wo_catch_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_incremental_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_fixed_n_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_fixed_time_fit.rds")
save_locations = c("utils/data/processed_data/model_fit/patch_leaving_models/comparison/a_loo_patch_leaving_model_baseline_lakes.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_patch_leaving_model_time_wo_catch_lakes.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_patch_leaving_model_incremental_lakes.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_patch_leaving_model_fixed_n_lakes.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_patch_leaving_model_fixed_time_lakes.rds")

for (i in 1:length(fit_locations)){
  compute_waic_lakes(fit_location = fit_locations[i],
               save_location = save_locations[i],
               stan_data = stan_data)
  gc()
  gc()
}

# conduct model comparison across lakes
comp_list = list()
for (i in 1:10){
  loos = list(readRDS(save_locations[[1]])[[i]],
              readRDS(save_locations[[2]])[[i]],
              readRDS(save_locations[[3]])[[i]],
              readRDS(save_locations[[4]])[[i]],
              readRDS(save_locations[[5]])[[i]])
  comp_list[[i]] = loo_compare(loos)
  comp_list[[i]] = as.data.frame(comp_list[[i]])
  comp_list[[i]]$lake = i
  comp_list[[i]]$model_id = rownames(comp_list[[i]])
}
comparison = do.call(rbind, comp_list)
comparison$model_name = ifelse(comparison$model_id == "model1", "Baseline", 
                               ifelse(comparison$model_id == "model2", "GUT",
                                      ifelse(comparison$model_id == "model3", "Incremental",
                                             ifelse(comparison$model_id == "model4", "Fixed-n", "Fixed-time"))))

# save results
saveRDS(comparison, file = "utils/data/processed_data/model_fit/patch_leaving_models/comparison/comparison_basic_strategies_lakes.rds")

################################################################################
# Extended strategies by lake
fit_locations =  c("utils/data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_time_wo_catch_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/c_patch_leaving_model_patch_discovery_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/d_patch_leaving_model_globallocal_updating_fit.rds",
                   "utils/data/processed_data/model_fit/patch_leaving_models/e_patch_leaving_model_bandwidth_fit.rds")
save_locations =  c("utils/data/processed_data/model_fit/patch_leaving_models/comparison/b_loo_patch_leaving_model_time_wo_catch_lake.rds",
                    "utils/data/processed_data/model_fit/patch_leaving_models/comparison/c_loo_patch_leaving_model_patch_discovery_lake.rds",
                    "utils/data/processed_data/model_fit/patch_leaving_models/comparison/d_loo_patch_leaving_model_globallocal_updating_lake.rds",
                    "utils/data/processed_data/model_fit/patch_leaving_models/comparison/e_loo_patch_leaving_model_bandwidth_lake.rds")

for (i in 1:length(fit_locations)){
  compute_waic_lakes(fit_location = fit_locations[i],
                     save_location = save_locations[i],
                     stan_data = stan_data)
  gc()
  gc()
}

# conduct model comparison across lakes
comp_list = list()
for (i in 1:10){
  loos = list(readRDS(save_locations[[1]])[[i]],
              readRDS(save_locations[[2]])[[i]],
              readRDS(save_locations[[3]])[[i]],
              readRDS(save_locations[[4]])[[i]])
  comp_list[[i]] = loo_compare(loos)
  comp_list[[i]] = as.data.frame(comp_list[[i]])
  comp_list[[i]]$lake = i
  comp_list[[i]]$model_id = rownames(comp_list[[i]])
}
comparison = do.call(rbind, comp_list)
comparison$model_name = ifelse(comparison$model_id == "model1", "GUT", 
                               ifelse(comparison$model_id == "model2", "+Fish\nDiscovery",
                                      ifelse(comparison$model_id == "model3", "+Updating","+Spatial\nFeatures")))

# save results
saveRDS(comparison, file = "utils/data/processed_data/model_fit/patch_leaving_models/comparison/comparison_extended_strategies_lakes.rds")

}
################################################################################
# END
################################################################################
