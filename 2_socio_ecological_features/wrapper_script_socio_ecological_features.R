################################################################################
#
# Title: Descriptive Analyses and Figures
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 29/11/2024 (DD/MM/YYYY)
#
################################################################################
run_2_socio_ecological_features <- function(wd){
# set working directory
if (getwd() != wd){
  setwd(wd)
}

# check
getwd()

# load library and functions
source("utils/library/library.R")

################################################################################
# load data
################################################################################

# data contains information for all foraging locations (ID equivalent to camera_id in other datasets)
spot_selection_data = fread(file = "utils/data/spot_selection_data.csv")
spot_selection_data = spot_selection_data %>% filter(!exclude)

# Import competition areas
competition_area = readRDS("utils/data/study_lakes.rds")

# polygon lakes for faster execution time
point_lakes = list()
for (i in 1:10){
  point_lakes[[i]] = competition_area[[i]] %>% 
    st_buffer(10) %>% 
    st_cast("POINT") %>% 
    as.data.frame()
  point_lakes[[i]] = point_lakes[[i]] %>% 
    mutate(x = unlist(map(point_lakes[[i]]$polygons,1)),
           y = unlist(map(point_lakes[[i]]$polygons,2))) %>% 
    dplyr::select(x,y)
}

# Import GPS Data
gps_data <- fread("utils/data/gps_data.csv")

# add unique ID as track identifier
spot_selection_data$unique_ID = paste(spot_selection_data$day, spot_selection_data$year, spot_selection_data$ID, paste = "")

# remove IDs with only one spot left
spot_selection_data = spot_selection_data %>% group_by(day, year, ID) %>% mutate(n = max(spot_id)) %>% filter(n > 1)

# compute spot id
spot_selection_data = spot_selection_data %>% 
  group_by(year, day, ID) %>% 
  mutate(spot_id = 1:n())

# arrange data chronologically
spot_selection_data = spot_selection_data %>% 
  arrange(year, day, ID, spot_id)

# how many alternative spots to simulate
N_spots = 50

# call simulation script
source("utils/library/simulate_available_spots.r")
simulate_available_spots(spot_selection_data, N_spots, gps_data, point_lakes, competition_area)

# extract socio-ecological features of simulated spots
source("utils/library/compute_raw_features.r")
compute_raw_features(gps_data, spot_selection_data)

# create stan data list from extracted features
source("utils/library/prepare_stan_data.r")
prepare_stan_data_spot_selection()

# load model
spot_selection_model <- cmdstan_model("utils/stan/spot_selection_model.stan",
                                      cpp_options = list(stan_threads = TRUE),
                                      stanc_options = list("O1"))

# load data
stan_data = readRDS(file = "utils/data/stan_data_spot_selection_model.rds")

# run model
spot_selection_model_fit <- spot_selection_model$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4, 
  refresh = 20,
  iter_sampling = 1000,
  iter_warmup = 2000,
  threads_per_chain = 8
)

# save output
spot_selection_model_fit$save_object(file = "utils/data/stan/spot_selection_model_fit.rds")

# simulate fitted spatial features
source("utils/library/simulate_fitted_features.r")
simulate_fitted_features()

# visualize results and store figures in output folder
source("utils/library/visualize_model_results.r")
visualize_spot_selection_model()

# run spot selection simulation
source("utils/library/simulate_spot_selection.r")
simulate_spot_selection()

}
