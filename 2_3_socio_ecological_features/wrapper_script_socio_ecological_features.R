################################################################################
#
# Title: Wrapper spatial model section 2
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 19/12/2025 (DD/MM/YYYY)
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
  
  # arrange data chronologically
  spot_selection_data = spot_selection_data %>% 
    arrange(year, day, camera_id, spot_id)
  
  # how many alternative spots to simulate
  N_spots = 50
  
  # call simulation script. Note that the feature matrix included in the utils file
  # is not identical to the one used in the paper. The code has been rerun to 
  # ensure reproducibility and new files were created. This means that the model
  # estimates and simulated locations in the spot-selection model can be diverging
  # slightly from those reported in the paper. We reran this script multiple times
  # to ensure that these divergences do not affect the direction, size or 
  # credibility of any of the effects reported in the paper.
  if (file.exists("utils/data/processed_data/feature_matrix_spot_selection.rds")){
    #
  } else {
    source("2_3_socio_ecological_features/scripts/simulate_available_spots.r")
    simulate_available_spots(spot_selection_data, N_spots, gps_data, point_lakes, competition_area)
  }
  
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
