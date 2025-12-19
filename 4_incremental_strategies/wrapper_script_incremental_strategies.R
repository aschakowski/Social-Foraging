################################################################################
#
# Title: Wrapper script patch leaving analyses
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 19/12/2025 (DD/MM/YYYY)
#
################################################################################
run_4_incremental_strategies <- function(wd){
  
  # load library and functions
  source("utils/library/library.R")
  
  ##############################################################################
  # run depletion models
  ##############################################################################
  source("4_incremental_strategies/scripts/run_depletion_analysis.r")
  run_depletion_analyses()
  
  source("4_incremental_strategies/scripts/visualize_depletion_results.r")
  visualize_depletion_results()
  
  ##############################################################################
  # run patch leaving models
  ##############################################################################
  
  # prepare stan data
  source("4_incremental_strategies/scripts/prepare_stan_data_patch_leaving.r")
  prepare_stan_data_patch_leaving()
  
  # prepare stan data
  source("4_incremental_strategies/scripts/prepare_stan_data_patch_leaving_grp_predictors.r")
  prepare_stan_data_patch_leaving_grp_predictors()
  
  # run models
  source("4_incremental_strategies/scripts/run_model_comparison_patch_leaving.r")
  run_model_comparison_patch_leaving()
  
  ##############################################################################
  # visualize patch leaving results
  ##############################################################################
  source("4_incremental_strategies/scripts/visualize_patch_leaving_model.r")
  visualize_patch_leaving()
  
  ##############################################################################
  # gut simulations supplement
  ##############################################################################
  #run simulations
  source("4_incremental_strategies/scripts/run_gut_simulations.r")
  run_gut_simulations()
  
  # visualize
  source("4_incremental_strategies/scripts/visualize_gut_simulations.r")
  visualize_gut_simulations()
  
  ##############################################################################
  # END
  ##############################################################################
}
