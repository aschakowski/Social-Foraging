# Social-Foraging

This repository contains the full data and scripts to reproduce all analyses and figures for:

Schakowski, A., Deffner, D., Kortet, R., Niemelä, P.T., Kavelaars, M.M., Monk, C.T., Pykälä, M., & Kurvers, R.H.J.M (2025) Socioecology drives adaptive social foraging dynamics in the wild.

Paper here: link to preprint

"wrapper_script_social_foraging_dynamics.r" sources all helper scripts and runs analyses. 
Scripts source raw data stored in "data/raw_data" and store processed data in the folder "data/processed_data".
Helper scripts for analyses are stored in the folder for the respective sub-paragraph of the manuscript. 
Figures for each sub-paragraph are stored in the "output" folder in the respective sub-directory.

Software requirements

The analysis code was written in R 4.0.3. Statistical models are fit using the Stan MCMC engine via the rstan (2.21.2) and cmdstanr (0.5.3) packages, which require a C++ compiler. Installation instructions are available at https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started and https://mc-stan.org/cmdstanr/articles/cmdstanr.html. See also the Stan user guide at https://mc-stan.org/users/documentation. The rethinking package (2.12) is required to process fitted model outputs (installation instructions at http://xcelab.net/rm/software/).
