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
run_2_3_socio_ecological_features <- function(wd){

################################################################################
# load data
################################################################################

# data contains information for all foraging locations (ID equivalent to camera_id in other datasets)
spot_selection_data = fread(file = "utils/data/raw_data/spot_selection_data.csv")
spot_selection_data = spot_selection_data %>% filter(!exclude)

# Import competition areas
competition_area = readRDS("utils/data/raw_data/lake_data/study_lakes.rds")

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
gps_data <- fread("utils/data/raw_data/gps_data.csv")

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
source("2_3_socio_ecological_features/scripts/simulate_available_spots.r")
simulate_available_spots(spot_selection_data, N_spots, gps_data, point_lakes, competition_area)

# extract socio-ecological features of simulated spots
source("2_3_socio_ecological_features/scripts/compute_raw_features.r")
extract_features(gps_data, spot_selection_data)

# create stan data list from extracted features
source("2_3_socio_ecological_features/scripts/prepare_stan_data_spot_selection.r")
prepare_stan_data_spot_selection()

# run model comparison spot selection models
source("2_3_socio_ecological_features/scripts/run_model_comparison_spot_selection.r")
run_model_comparison_spot_selection()

# simulate fitted spatial features from full model
source("2_3_socio_ecological_features/scripts/simulate_fitted_features.r")
simulate_fitted_features()

# visualize results and store figures in output folder
source("2_3_socio_ecological_features/scripts/visualize_model_results.r")
visualize_spot_selection_model()

# run spot selection simulation
source("2_3_socio_ecological_features/scripts/simulate_spot_selection.r")
simulate_spot_selection(spot_selection_data, gps_data)

# visualize spot selection simulation
source("2_3_socio_ecological_features/scripts/visualize_spot_selection_simulation.r")
visualize_spot_selection_simulation()


}
################################################################################
# END
################################################################################
