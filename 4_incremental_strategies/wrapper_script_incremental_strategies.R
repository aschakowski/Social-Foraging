################################################################################
#
# Title: Wrapper script patch leaving analyses
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 14/02/2025 (DD/MM/YYYY)
#
################################################################################
run_4_incremental_strategies <- function(){
  
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
  # visualize
  source("4_incremental_strategies/scripts/visualize_gut_simulations.r")
  visualize_gut_simulations()
  
  ##############################################################################
  # END
  ##############################################################################
}
