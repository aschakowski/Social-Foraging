################################################################################
#
# Title: Wrapper script interindividual differences
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 29/11/2024 (DD/MM/YYYY)
#
################################################################################
run_5_inter_individual_differences <- function(wd){
  # set working directory
  if (getwd() != wd){
    setwd(wd)
  }
  
  # check
  getwd()
  
  # load library and functions
  source("utils/library/library.R")
  
  ##############################################################################
  # run analyses and make visualizations
  ##############################################################################
  source("5_inter_individual_differences/scripts/run_int_diff_analyses.r")
  analyses_interindividual_differences()
  
  ##############################################################################
  # END
  ##############################################################################
}
