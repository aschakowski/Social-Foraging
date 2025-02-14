################################################################################
#
# Title: Run and compare spot selection models
#
# Authors: Schakowski, A.
#
# Last updated: 14/02/2025 (DD/MM/YYYY)
#
################################################################################
run_model_comparison_spot_selection <- function(){
  
################################################################################
# helper function

# fit models
model_fit <- function(model_location, save_location, data_list, warmup, iterations, chains, fixed_param){
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
    threads_per_chain = round(32/chains),
    fixed_param = fixed_param
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
  
  # compute loo
  #loo <- model_fit$loo(cores = 40)
  
  ll = model_fit$draws(variables = "log_lik")
  waic = waic(ll)
  
  # save loo
  saveRDS(waic, file = save_location)
  
  # keep workspace clear
  gc()
  
  # empty return
  return()
  
}

################################################################################
# load stan data
stan_data = readRDS(file = "utils/data/raw_data/stan_data_spot_selection.rds")

# standardized time since start from 0 to 1
stan_data$Time_since_start = as.numeric(stan_data$Time_since_start)
stan_data$Time_since_start = stan_data$Time_since_start / max(stan_data$Time_since_start)

# any 0s in dataset? If so, replace with small value
stan_data$Distances_social[stan_data$Distances_social == 0] <- .01
stan_data$Distances_roughness[stan_data$Distances_roughness == 0] <- .01
stan_data$Distances_success[stan_data$Distances_success == 0] <- .01
stan_data$Distances_loss[stan_data$Distances_loss == 0] <- .01
stan_data$Locality[stan_data$Locality == 0] <- .01
min(stan_data$Time_since_start)
max(stan_data$Time_since_start)
sum(is.na(stan_data$Time_since_start))

# any missings in dataset?
stan_data$Distances_social[is.na(stan_data$Distances_social)]
stan_data$Distances_roughness[is.na(stan_data$Distances_roughness)]
stan_data$Distances_success[is.na(stan_data$Distances_success)]
stan_data$Distances_loss[is.na(stan_data$Distances_loss)]
stan_data$Locality[is.na(stan_data$Locality)]

# fit models
model_locations = c("utils/stan/spot_selection_models/a_random_selection.stan",
                    "utils/stan/spot_selection_models/b_no_roughness.stan",
                    "utils/stan/spot_selection_models/c_no_personal.stan",
                    "utils/stan/spot_selection_models/d_no_social.stan",
                    "utils/stan/spot_selection_models/e_all_features.stan",
                    "utils/stan/spot_selection_models/f_conditional.stan")
save_locations = c("utils/data/processed_data/model_fit/spot_selection_models/a_random_selection.rds",
                   "utils/data/processed_data/model_fit/spot_selection_models/b_no_roughness.rds",
                   "utils/data/processed_data/model_fit/spot_selection_models/c_no_personal.rds",
                   "utils/data/processed_data/model_fit/spot_selection_models/d_no_social.rds",
                   "utils/data/processed_data/model_fit/spot_selection_models/e_all_features.rds",
                   "utils/data/processed_data/model_fit/spot_selection_models/f_conditional.rds")

for (i in 1:length(model_locations)){
  
  if (i == 1){
    model_fit(model_location = model_locations[i],
              save_location = save_locations[i],
              data_list = stan_data,
              warmup = 1000,
              iterations = 2000,
              chains = 4,
              fixed_param = T)
  } else {
    model_fit(model_location = model_locations[i],
              save_location = save_locations[i],
              data_list = stan_data,
              warmup = 1000,
              iterations = 2000,
              chains = 4,
              fixed_param = F)
  }
  
}

# run model comparison
fit_locations = save_locations
save_locations = c("utils/data/processed_data/model_fit/spot_selection_models/comparison/a_random_selection.rds",
                   "utils/data/processed_data/model_fit/spot_selection_models/comparison/b_no_roughness.rds",
                   "utils/data/processed_data/model_fit/spot_selection_models/comparison/c_no_personal.rds",
                   "utils/data/processed_data/model_fit/spot_selection_models/comparison/d_no_social.rds",
                   "utils/data/processed_data/model_fit/spot_selection_models/comparison/e_all_features.rds",
                   "utils/data/processed_data/model_fit/spot_selection_models/comparison/f_conditional.rds")

for (i in 1:length(fit_locations)){
  compute_waic(fit_location = fit_locations[i],
               save_location = save_locations[i])
}

# conduct model comparison
loo_list = list(readRDS(save_locations[1]), 
                readRDS(save_locations[2]), 
                readRDS(save_locations[3]), 
                readRDS(save_locations[4]), 
                readRDS(save_locations[5]),
                readRDS(save_locations[6]))
comparison = loo_compare(loo_list)

# save results
saveRDS(comparison, file = "utils/data/processed_data/model_fit/spot_selection_models/comparison/comparison.rds")

}
################################################################################
# END
################################################################################
