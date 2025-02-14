################################################################################
#
# Title: Wrapper Script Social Foraging Dynamics Icefishing
#
# Description: 
#
# Authors: Schakowski, A.
#
# Last updated: 14/02/2025 (DD/MM/YYYY)
#
################################################################################
# set working directory
if (Sys.info()[1] == "Linux"){
  wd="add your working directory here; this was run on a fileserver with different wd depending on system"
} else {
  wd="add your working directory here"
}

if (getwd() != wd){
  setwd(wd)
}

# check
getwd()

# load library and functions
source("utils/library/library.R")

# brms option
options(buildtools.check = function(action) TRUE)
################################################################################
# set seed for whole analysis
set.seed(28012024)

################################################################################
# 1. Ice-fishing competitions to study human social foraging
################################################################################

# The script loads data for descriptive analyses and outputs figures
source("1_icefishing_social_foraging/wrapper_script_icefishing_social_foraging.R")

# This function runs the analyses and creates output. To reproduce analyses,
# script can be run on its own.
run_1_icefishing_social_foraging()

################################################################################
# 2. Socio-ecological features drive adaptive spatial search
################################################################################
################################################################################
# 3. Social information use modulates area--restricted search
################################################################################

# The script runs spot selection model, model simulation and visualizes results.
source("2_3_socio_ecological_features/wrapper_script_socio_ecological_features.R")

# To see results of analyses, run script on its own.
run_2_3_socio_ecological_features()

################################################################################
# 4. Non-social heuristics drive patch leaving decisions
################################################################################

# The script runs depletion models and patch--leaving model.
# To see intermediary results, run script stand-alone.
source("4_incremental_strategies/wrapper_script_incremental_strategies.R")
run_4_incremental_strategies()

################################################################################
# 5. Inter--individual differences in foraging strategies
################################################################################

# The scripts investigates inter--individual differences and produces figures.
# To see results, run script on its own.
source("5_inter_individual_differences/wrapper_script_inter_ind_differences.R")
run_5_inter_individual_differences()

################################################################################
# END
################################################################################

