################################################################################
#
# Title: Run depletion models
#
# Authors: Schakowski, A.
#
# Last updated: 27/07/2024 (DD/MM/YYYY)
#
################################################################################
run_depletion_analyses <- function(){

# load patch leaving data
patch_leaving_data = fread(file = "utils/data/raw_data/patch_leaving_data.csv")
patch_leaving_data = patch_leaving_data %>% filter(!exclude) %>% filter(!exclude2)

# add unique participant id
catch_data <- read_delim("utils/data/raw_data/catch_data.csv", 
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
catch_data = catch_data %>% dplyr::select(camera_id, day, year, participant_id)
patch_leaving_data = left_join(patch_leaving_data, catch_data)

# depletion dat
depletion_dat = patch_leaving_data %>% 
  group_by(day, year, camera_id, spot_id) %>%
  mutate(cumulative_r = cumsum(reward)) %>% 
  mutate(s = sum(reward)) %>% 
  filter(s > 1) %>% 
  group_by(day, year, camera_id, spot_id) %>% 
  mutate(test = n()) %>% 
  group_by(day, year, camera_id, spot_id) %>% 
  arrange(day, year, camera_id, spot_id, bin)

# unique spot id
depletion_dat$unique_spot_id = as.numeric(as.factor(paste(depletion_dat$day, depletion_dat$year, depletion_dat$camera_id, depletion_dat$spot_id)))

# test
depletion_dat %>%
  filter(unique_spot_id == 230) %>% 
  ggplot(aes(x = bin, y = cumulative_r)) + 
  geom_point() + 
  geom_path()

# add time for within chain parallelization
depletion_dat$time = 1:nrow(depletion_dat)

# hierarchical no pooling
stan_data_depletion = list(T = nrow(depletion_dat),
                           cumulative_R = depletion_dat$cumulative_r,
                           bin = depletion_dat$bin-1,
                           spot_id = depletion_dat$unique_spot_id,
                           N_spots = length(unique(depletion_dat$unique_spot_id)),
                           time = depletion_dat$time,
                           grainsize = 20)

################################################################################
# Path depletion models
exponential_model_hierarchical = cmdstan_model("utils/stan/depletion_models/exponential_depletion_nopooling.stan",
                                               stanc_options = list("O1"),cpp_options = list(stan_threads = TRUE))
michaelis_menten_model_hierarchical = cmdstan_model("utils/stan/depletion_models/michaelis_menten_depletion_nopooling.stan",
                                                    stanc_options = list("O1"),cpp_options = list(stan_threads = TRUE))
linear_model_hierarchical = cmdstan_model("utils/stan/depletion_models/linear_depletion_nopooling.stan",
                                          stanc_options = list("O1"),cpp_options = list(stan_threads = TRUE))
holling_model_hierarchical = cmdstan_model("utils/stan/depletion_models/holling_depletion_nopooling.stan",
                                           stanc_options = list("O1"),cpp_options = list(stan_threads = TRUE))
################################################################################
# fit depletion curves
holling <- holling_model_hierarchical$sample(
  data = stan_data_depletion,
  chains = 4,
  parallel_chains = 4,
  refresh = 2,
  iter_sampling = 1000,
  iter_warmup = 1000,
  threads_per_chain = 8
)

exponential <- exponential_model_hierarchical$sample(
  data = stan_data_depletion,
  chains = 4,
  parallel_chains = 4,
  refresh = 2,
  iter_sampling = 1000,
  iter_warmup = 1000,
  threads_per_chain = 8
)

michaelis_menten <- michaelis_menten_model_hierarchical$sample(
  data = stan_data_depletion,
  chains = 4,
  parallel_chains = 4,
  refresh = 2,
  iter_sampling = 1000,
  iter_warmup = 1000,
  threads_per_chain = 8
)

linear <- linear_model_hierarchical$sample(
  data = stan_data_depletion,
  chains = 4,
  parallel_chains = 4,
  refresh = 2,
  iter_sampling = 1000,
  iter_warmup = 1000,
  threads_per_chain = 8
)

################################################################################
# save models
holling$save_object(file = "utils/data/processed_data/model_fit/depletion_models/holling_depletion.rds")
exponential$save_object(file = "utils/data/processed_data/model_fit/depletion_models/exponential_depletion.rds")
linear$save_object(file = "utils/data/processed_data/model_fit/depletion_models/linear_depletion.rds")
michaelis_menten$save_object(file = "utils/data/processed_data/model_fit/depletion_models/michaelis_menten_depletion.rds")

rm(list = ls())
gc()

################################################################################
# Model Comparison
fit_locations = c("utils/data/processed_data/model_fit/depletion_models/holling_depletion.rds",
                  "utils/data/processed_data//model_fit/depletion_models/exponential_depletion.rds",
                  "utils/data/processed_data//model_fit/depletion_models/linear_depletion.rds",
                  "utils/data/processed_data//model_fit/depletion_models/michaelis_menten_depletion.rds")

# preallocate waic list
loo_list = list()
for (i in 1:length(fit_locations)){
  fit = readRDS(fit_locations[i])
  loo_list[[i]] = waic(fit$draws(variables = "log_lik"))
}

# make comparison
comparison = loo_compare(loo_list)

# save comparison result
saveRDS(comparison, file = "utils/data/processed_data/model_fit/depletion_models/comparison/depletion_comparison.rds")

################################################################################
# convergence diagnostics
draws_exponential = exponential$draws(format = "df")
rhats = linear$summary()
rhats = rhats[, "rhat"]
max(rhats)
################################################################################

################################################################################
# END
################################################################################
}