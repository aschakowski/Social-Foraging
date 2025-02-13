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

# run model comparison spot selection models
source("2_3_socio_ecological_features/scripts/run_model_comparison_spot_selection.r")
run_model_comparison_spot_selection()

# visualize results and store figures in output folder
source("2_3_socio_ecological_features/scripts/visualize_model_results.r")
visualize_spot_selection_model()

# visualize spot selection simulation
source("2_3_socio_ecological_features/scripts/visualize_spot_selection_simulation.r")
visualize_spot_selection_simulation()


}
################################################################################
# END
################################################################################
