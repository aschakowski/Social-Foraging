################################################################################
#
# Title: Wrapper Script Social Foraging Dynamics Icefishing
#
# Description: 
#
# Authors: Schakowski, A.
#
# Last updated: 29/11/2024 (DD/MM/YYYY)
#
################################################################################
# set working directory
wd="P:/IceFishing/icefishing_git"

if (getwd() != wd){
  setwd(wd)
}

# check
getwd()
################################################################################

################################################################################
# 1. Ice-fishing competitions to study human social foraging
################################################################################

# The script loads data for descriptive analyses and outputs figures
source("1_icefishing_social_foraging/scripts/wrapper_script_icefishing_social_foraging.R")

# This function runs the analyses and creates output. To reproduce analyses,
# script can be run on its own.
run_1_icefishing_social_foraging(wd)

################################################################################
# 2. Socio-ecological features drive adaptive spatial search
################################################################################
################################################################################
# 3. Social information use modulates area--restricted search
################################################################################

# The script runs spot selection model, model simulation and visualizes results.
source("2_socio_ecological_features/scripts/wrapper_script_socio_ecological_features.R")

# To see results of analyses, run script on its own.
run_2_socio_ecological_features(wd)

################################################################################
# 4. Non-social heuristics drive patch leaving decisions
################################################################################

# The script runs depletion models and patch--leaving model.
# To see intermediary results, run script stand-alone.

################################################################################
# 5. Inter--individual differences in foraging strategies
################################################################################

# The scripts investigates inter--individual differences and produces figures.
# To see results, run script on its own.

################################################################################
# END
################################################################################

