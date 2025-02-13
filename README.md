# Social-Foraging

This repository contains the data and scripts necessary to reproduce all analyses and figures for:

Schakowski, A., Deffner, D., Kortet, R., Niemelä, P.T., Kavelaars, M.M., Monk, C.T., Pykälä, M., & Kurvers, R.H.J.M (2025) Socioecology drives adaptive social foraging dynamics in the wild.

Paper here: link to preprint

The raw GPS and video data is available upon request. 

To reproduce analyses, an additional folder "utils" has to be downloaded from OSF (), and placed into the same working directory as the remaining folders. Before running the analyses, an additional empty folder named "output" has to be placed in each subdirectory of the 4 main folder (1_..., 2_3_..., 4_..., and 5_...).

"wrapper_script_social_foraging_dynamics.r" sources all helper scripts and runs analyses. 
Scripts source raw data stored in "utils/data/raw_data" and store processed data in the folder "utils/data/processed_data".
"utils" directory needs to be unzipped before running analyses.
"utils/library/" contains helper functions and loads necessary packages.
Helper scripts for analyses are stored in the folder for the respective sub-paragraph of the manuscript. 
Figures for each sub-paragraph are stored in the "output" folder in the respective sub-directory.

Software requirements

The analysis code was written in R (v4.0.3.) Statistical models are fit using the Stan MCMC engine (v2.32.6) and cmdstanr (v0.8.1) packages, which require a C++ compiler. Installation instructions are available at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started and https://mc-stan.org/cmdstanr/articles/cmdstanr.html. See also the Stan user guide at https://mc-stan.org/users/documentation. All other required R packages can be found in "utils/library/library.R"

Detailed package versions can be found in "session_info.txt".
