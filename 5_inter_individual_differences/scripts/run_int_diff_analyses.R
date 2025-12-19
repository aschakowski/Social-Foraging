################################################################################
#
# Title: Run and compare depletion models
#
# Authors: Schakowski, A.
#
# Last updated: 27/07/2024 (DD/MM/YYYY)
#
################################################################################
analyses_interindividual_differences <- function(){

library("ggsci")
detach("package:weights", unload=TRUE)
################################################################################
# Part 0: Load Demographic Data and Catch Results
################################################################################
# load participant data and demographics
catch_data <- fread("utils/data/raw_data/catch_data.csv")
  
# survey data
survey_data <- fread("utils/data/raw_data/survey_data.csv")

################################################################################
# Part 1: Analyses Spatial Model
################################################################################
# load spatial model fit
spatial_model = readRDS(file = "utils/data/processed_data/model_fit/spot_selection_models/f_conditional.rds")
identifiers = readRDS(file = "utils/data/processed_data/identifiers_model_key.rds")

################################################################################
# 1a. Variance Decomposition Spatial Strategies
################################################################################
variance_decomposition_spatial = spatial_model %>% 
  spread_draws(sigma_id[parameter],
               sigma_trip[parameter],
               sigma_lake[parameter]) %>% 
  filter(parameter < 5) %>% 
  mutate(Id = sigma_id^2 / (sigma_trip^2 + sigma_lake ^2 + sigma_id ^2),
         Lake = sigma_lake^2 / (sigma_trip^2 + sigma_lake ^2 + sigma_id ^2)) %>% 
  pivot_longer(cols = c(Id, Lake), names_to = "Component", values_to = "est") %>% 
  mutate(feature = parameter %>% factor(levels = c("1", "2", "3", "4")) %>% fct_recode("Social" = "1","Roughness"="2", "Successful" = "3", "Unsuccessful" = "4")) %>% 
  mutate(name = factor(feature, levels=c("Roughness", "Social", "Unsuccessful", "Successful"))) %>%
  ggplot(aes(x = name, y = est, color = Component, group = Component)) + 
  stat_pointinterval(position = position_dodge(.2),size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  #scale_x_discrete(labels = rev(c("Social", "Roughness", "Success", "Loss"))) + 
  coord_flip(ylim = c(0, 1))+ 
  ylab("Posterior Estimate (95% CrI)\nProp. Expl. Variance~Component")+
  labs(tag = "A", color = "Component  ")+
  ggtitle("Variance Decomposition")+
  xlab("Feature")+
  scale_color_manual(values = c("green4", "lightblue3"), breaks = c("Lake", "Id"))+
  #scale_color_npg()+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, angle = 51.5, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans"))
variance_decomposition_spatial

################################################################################
# 1b. Correlations parameters (loo analysis)
################################################################################
# load id, trip, lake offsets
id_offsets = spatial_model %>% 
  spread_draws(v_id[id_index, parameter])
trip_offsets = spatial_model %>% 
  spread_draws(v_trip[trip_index, parameter])
lake_offsets = spatial_model %>% 
  spread_draws(v_lake[lake_index, parameter])

# join together
trip_offsets = left_join(trip_offsets, identifiers %>% group_by(trip_index) %>% slice(1) %>% dplyr::select(day, year, participant_id, track_id, lake_id, camera_id, lake_index, trip_index, id_index))
trip_offsets = left_join(trip_offsets, lake_offsets)
trip_offsets = left_join(trip_offsets, id_offsets)

# mutate
trip_offsets$trip_estimate = trip_offsets$v_lake + trip_offsets$v_id + trip_offsets$v_trip

# loop through all lakes
if (!file.exists("utils/data/processed_data/cors_spatial.rds")){
Cors_res = list()
for (i in 1:10){
  
  # make two sets
  sub1 = trip_offsets %>% 
    filter(lake_index == i) %>% 
    group_by(id_index, parameter, .draw) %>% 
    dplyr::summarize(estimate = mean(trip_estimate))
  
  sub2 = trip_offsets %>% 
    filter(lake_index != i) %>% 
    group_by(id_index, parameter, .draw) %>% 
    dplyr::summarize(estimate = mean(trip_estimate))
  
  # which ids in first set
  ids = unique(sub1$id_index)
  
  # filter ids  in set 2
  sub2 = sub2 %>% filter(id_index %in% ids)
  
  # filter in other direction
  sub1 = sub1 %>% filter(id_index %in% unique(sub2$id_index))
  
  # for each parameter and each draw, compute correlation
  max_dr = max(sub2$.draw)
  cors = matrix(nrow = max_dr, ncol = 8)
  for (par in 1:8){
    
    for (dr in 1:max_dr){
      
      cors[dr, par] = cor(sub1[sub1$.draw == dr & sub1$parameter == par,"estimate"],
                          sub2[sub2$.draw == dr & sub2$parameter == par,"estimate"])
      
    }
    
  }
  
  # Save tidy dataframe
  cors = as.data.frame(cors)
  colnames(cors) = paste("parameter_", 1:8, sep = "")
  cors$draw = 1:max_dr
  cors$N1 = length(unique(sub2$id_index))
  cors$N2 = length(unique(sub1$id_index))
  cors$lake = i
  
  Cors_res[[i]] <- cors
}
cors_spatial = do.call(rbind, Cors_res)

saveRDS(cors_spatial, file = "utils/data/processed_data/cors_spatial.rds")
}

# plot
cors_spatial = readRDS("utils/data/processed_data/cors_spatial.rds")
cors_spatial %>% 
  group_by(lake) %>% 
  dplyr::summarize(N = mean(N1)) %>% 
  dplyr::summarize(N = mean(N))
cors_spatial_plot = cors_spatial %>%
  pivot_longer(1:8) %>% 
  filter(name %in% c("parameter_1",
                     "parameter_2",
                     "parameter_3",
                     "parameter_4")) %>% 
  mutate(feature = name %>% factor(levels = c("parameter_2", "parameter_1", "parameter_4", "parameter_3")) %>% fct_recode("Social" = "parameter_1","Roughness"="parameter_2", "Successful" = "parameter_3", "Unsuccessful" = "parameter_4")) %>% 
  ggplot(aes(x = feature, y = value, group = lake)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = .6)+
  stat_pointinterval(position = position_dodge(.9), size = .6, point_size = 1, alpha = .3) + 
  scale_size_continuous(range = c(3, 10))+
  #scale_x_discrete(labels = rev(c("Social", "Roughness", "Success", "Loss"))) + 
  coord_flip(ylim = c(-1, 1))+ 
  ylab("Posterior Estimate (95% CrI)\nCorrelation Ind. Est. Across Lakes")+
  #ggtitle("Individual Consistency")+
  xlab("Feature")+
  #scale_color_manual(values = c("green4", "lightblue3"), breaks = c("Lake", "Id"))+
  #scale_color_npg()+
  labs(tag = "A")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, angle = 50.5, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.6, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_blank(),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans"))
cors_spatial_plot

################################################################################
# 1c. Predicting Ind. level strategies from person variables
################################################################################
id_effects = spatial_model %>% 
  spread_draws(v_id[id_index, parameter])

# join participant ids
id_effects = left_join(id_effects, identifiers %>% group_by(participant_id) %>% filter(spot_id >1) %>% slice(1) %>% dplyr::select(participant_id, id_index))

# join id level predictors to identifiers
catch_data <- read_csv("utils/data/raw_data/catch_data.csv")
variables = catch_data %>% 
  group_by(participant_id) %>% 
  slice(1) %>% 
  dplyr::select(sex, year_of_birth)
variables$age = 2023-variables$year_of_birth #note. for participants who participated in 2022, age will be overestimated. 
#But, to assign one age variable to each person (people participated across years, so some people might in theory have two age variables) we use 2023 as ref.
variables$sex = as.numeric(as.factor(variables$sex))-1
variables$participant_id = as.character(variables$participant_id)
identifiers = left_join(identifiers, variables)

survey_data = read_csv("utils/data/raw_data/survey_data.csv")
skill = survey_data %>% dplyr::select(participant_id, angling_skill)
skill$participant_id = as.character(skill$participant_id)
dem_dat = left_join(skill, variables) %>% dplyr::select(participant_id, angling_skill, sex, age)

# join demographic data to original dataset, only keep participants with full info
id_effects = left_join(id_effects, dem_dat)
id_effects = id_effects[!is.na(id_effects$age),]
id_effects = id_effects[!is.na(id_effects$angling_skill),]
id_effects = id_effects[!is.na(id_effects$sex),]

# predict individual level effects from demographic data
if (!file.exists("utils/data/processed_data/ind_regression_spatial.rds")){
start_time = Sys.time()
regression_list = list()
for (draw in unique(id_effects$.draw)){
  
  b = matrix(nrow = 4, ncol = 4)
  for (par_tmp in 1:4){
    
    # subset draw
    sub = id_effects[id_effects$.draw == draw & id_effects$parameter == par_tmp,]
    
    # compute X and y
    X = cbind(1, scale(sub$age)[,1], sub$sex, scale(sub$angling_skill)[,1])
    y = scale(sub$v_id)[,1]
    
    # compute normal equation
    b[par_tmp,] <- t(solve(t(X)%*%X)%*%t(X)%*%y)
    
  }
  
  # tidy results
  colnames(b) = c("Intercept", "age", "sex (ref. = m)", "skill")
  b = as.data.frame(b)
  b$parameter = 1:4
  b$draw = draw
  b$N = nrow(sub)
  
  # save in list
  regression_list[[draw]] = b
  
}
end_time = Sys.time()
end_time - start_time
regression_list = do.call(rbind, regression_list)
parameter_key = data.frame(parameter = 1:4, name = c("Social", "Roughness", "Successful", "Unsuccessful"))
regression_list_plot = left_join(regression_list, parameter_key)
regression_list_plot = regression_list_plot %>% pivot_longer(2:4, names_to = "Predictor", values_to = "Estimate")
# mislabeled reference category in computation loop. relabel here. (0 is female, 1 is male)
regression_list_plot = regression_list_plot %>% mutate(Predictor = ifelse(Predictor == "sex (ref. = m)", "gender\n(ref. = f)", Predictor))

# save results
saveRDS(regression_list_plot, "utils/data/processed_data/ind_regression_spatial.rds")
}
regression_list_plot = readRDS("utils/data/processed_data/ind_regression_spatial.rds")


# note: exact numerical parameter estimates not reproducible, because of randomness in simulation 
# of alternative spots used to fit spot selection model. 
# Model was fit multiple times to ensure that the direction, size and credibility
# of effects is not affected by this randomness.
# table
regression_list_plot %>% 
  group_by(name, Predictor) %>% 
  dplyr::summarize(mean_est = mean(Estimate),
                   lower = quantile(Estimate, .025),
                   upper = quantile(Estimate, .975))

# regression results
library(ggsci)
regression_plot_spatial = regression_list_plot %>% 
  filter(name %in% c("Social", "Roughness", "Successful", "Unsuccessful")) %>% 
  mutate(name = factor(name, levels = rev(c("Successful", "Unsuccessful", "Social", "Roughness")))) %>% 
  ggplot(aes(x = name, y = Estimate, group = Predictor, color = Predictor)) + 
  stat_pointinterval(position = position_dodge(.4),size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  #scale_x_discrete(labels = rev(c("Social", "Roughness", "Success", "Loss"))) + 
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  scale_color_npg(breaks = c("skill", "gender\n(ref. = f)", "age")) + 
  coord_flip(ylim = c(-1.5, 1.5)) +
  scale_y_continuous(breaks = seq(from = -1.5, to = 1.5, by = .5))+
  geom_hline(yintercept = 0, linetype = "dotted") + 
  #xlab("Feature") + 
  ylab("Posterior Estimate (95% CrI)\nInd. Feature Weight ~ Predictor") + 
  labs(tag = "C", color = "Predictor")+
  ggtitle("Individual Variation")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans"),
        legend.box.margin = margin(10, 10, 10, 10))
regression_plot_spatial

################################################################################
# 1d. Predicting Success From individual level strategies
################################################################################
# individual level parameters
ind_effects = spatial_model %>% 
  spread_draws(v_id[model_id, parameter])

# average individual level rank
catch_data$participant_id = as.character(catch_data$participant_id)
catch_data = left_join(catch_data, identifiers %>% group_by(participant_id) %>% filter(spot_id >1) %>% slice(1) %>% dplyr::select(participant_id, model_id = id_index))
rank_data = catch_data %>% 
  group_by(day, year) %>% 
  mutate(rank = rank(catch)/n(),
         catch_s = scale(catch)[, 1]) %>% 
  group_by(model_id) %>% 
  dplyr::summarize(rank = mean(rank),
                   catch = mean(catch_s))
rank_data$scaled_rank = scale(rank_data$rank)[,1]
rank_data = rank_data %>% filter(!is.na(model_id))
ind_effects = left_join(ind_effects, rank_data)
# compute posterior regression results predicting success in kg from strategies
# same results for rank (not shown in paper)
start_time = Sys.time()
b = matrix(nrow = max(ind_effects$.draw), ncol = 5)
for (draw in unique(ind_effects$.draw)){
  
  # subset draw
  sub = ind_effects %>% 
    filter(.draw == draw) %>% 
    pivot_wider(names_from = parameter, values_from = v_id)
  
  # compute X and y
  X = cbind(1, scale(sub[, 8:11])[, 1:4])
  y = scale(sub$catch)[,1]
  
  # compute normal equation
  b[draw,] <- t(solve(t(X)%*%X)%*%t(X)%*%y)
  
}

# tidy results
colnames(b) = c("Intercept", "Social", "Roughness", "Successful", "Unsuccessful")
b = as.data.frame(b)
b$draw = 1:draw
end_time = Sys.time()
end_time - start_time

success_prediction_spatial = b %>% 
  tidyr::pivot_longer(2:5) %>% 
  mutate(name = factor(name, levels = c("Social", "Roughness", "Unsuccessful", "Successful"))) %>% 
  ggplot(aes(x = name, y = value)) + 
  stat_pointinterval(size = .6, point_size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  coord_flip() + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ylab("Posterior Estimate (95% CrI)\nForaging Success ~ Ind. Feature Weight") + 
  xlab("Feature") + 
  labs(tag = "A")+
  #ggtitle("Foraging Success")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, angle = 55, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.2, .5),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans"),
        legend.box.margin = margin(0, 0, 0, 0))
success_prediction_spatial

################################################################################
# 1f. ICC catch success
################################################################################
# icc is computed as relative variance in success explained by id and lake intercepts.
# to account for skewed distribution, alternative is to use lognormal + zero inflation, 
# (hurdle logn) this yields similar results (not shown in paper).
mod_catch = brm(data = catch_data, 
                formula = catch/1000 ~ (1|participant_id) + (1|lake_id))
mod_catch_ln = brm(data = catch_data, 
                   formula = bf(catch/1000 ~ (1|participant_id) + (1|lake_id),
                                hu ~ (1|participant_id) + (1|lake_id)),
                   family = hurdle_lognormal())

# icc is computed as relative variance divided by total var. + res.
icc_catch = mod_catch %>% 
  spread_draws(sd_lake_id__Intercept,
               sd_participant_id__Intercept,
               sigma) %>% 
  mutate(var_lake = sd_lake_id__Intercept^2 / (sd_lake_id__Intercept^2 + sd_participant_id__Intercept^2 + sigma^2),
         var_id = sd_participant_id__Intercept^2 / (sd_lake_id__Intercept^2 + sd_participant_id__Intercept^2 + sigma^2)) %>% 
  pivot_longer(7:8) %>% 
  ggplot(aes(x = name, y = value)) + 
  stat_pointinterval(size = .6, point_size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = c("Individual", "Lake")) + 
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  #scale_color_npg(breaks = c("skill", "gender (ref. = f)", "age")) + 
  coord_flip() + 
  labs(tag = "C")+
  ylab("Posterior Estimate (95% CrI)\nExpl. Variance Foraging Success ~ Component") + 
  xlab("Cmoponent") + 
  ggtitle("Repeatability Foraging Success")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, angle = 55, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.2, .5),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_blank(),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans"),
        legend.box.margin = margin(10, 10, 10, 10))
icc_catch

################################################################################
# 1g. Predicting Lake--level variation by lake--characteristics
################################################################################
# lake characteristics (need model data)
spot_selection_data = read.csv("utils/data/raw_data/spot_selection_data.csv")
spot_selection_data = spot_selection_data %>% 
  filter(!exclude)
identifiers = readRDS(file = "utils/data/processed_data/identifiers_model_key.rds")
key = identifiers %>% group_by(lake_index) %>% slice(1) %>% dplyr::select(day, year, lake_index)
key = key %>% filter(!is.na(lake_index)) #missing lake index for first spots - filter out

# lake estimates
lake_estimates = spatial_model %>% 
  spread_draws(v_lake[lake_index, parameter])
lake_estimates = left_join(lake_estimates, key)

# 1. ratio of successful angling spots [%]
ratio_successful = spot_selection_data %>% group_by(day, year) %>% dplyr::summarize(ratio_successful = mean(Return))
lake_estimates$day = as.numeric(lake_estimates$day)
lake_estimates$year = as.numeric(lake_estimates$year)
ratio_successful = left_join(lake_estimates, ratio_successful)

# 2. average amount of fish caught [kg]
avg_catch = catch_data %>% group_by(day, year) %>% dplyr::summarize(avg_catch = mean(catch))
avg_catch = left_join(ratio_successful, avg_catch)
test_data = avg_catch

# 3. posterior estimate, v_lake from lake characteristics
if (!file.exists("utils/data/processed_data/eco_regression_spatial.rds")){
results = list()
for (j in 1:4){
  
  b = matrix(ncol = 2, nrow = max(test_data$.draw))
  
  for (i in 1:max(test_data$.draw)){
    # subset draw
    sub = test_data %>% 
      filter(.draw == i & parameter == j)
    
    for (k in 1:2){
      
      # compute X and y
      X = cbind(1, scale(sub$avg_catch)[, 1],
                scale(sub$ratio_successful)[, 1])
      
      X = X[, c(1, (k+1))]
      y = scale(sub$v_lake)[,1]
      
      # compute normal equation
      b[i,k] <- (t(solve(t(X)%*%X)%*%t(X)%*%y))[, 2]
      
    }
    
  }
  
  b = as.data.frame(b)
  colnames(b) <- c("% Successful", "Avg. Catch (kg)     ")
  b$parameter = j
  b$draw = 1:max(test_data$.draw)
  
  results[[j]] <- b
  
}
results = do.call(rbind, results)

# make figure
parameter_key = data.frame(parameter = 1:8, name = c("Social", "Roughness", "Successful", "Unsuccessful"))
regression_list_plot_eco = left_join(results, parameter_key)
regression_list_plot_eco = regression_list_plot_eco %>% pivot_longer(1:2, names_to = "Predictor", values_to = "Estimate")

# save
saveRDS(regression_list_plot_eco, file = "utils/data/processed_data/eco_regression_spatial.rds")
} else {
  regression_list_plot_eco = readRDS(file = "utils/data/processed_data/eco_regression_spatial.rds")
}
################################################################################
# note. exact numerical estimates not reproducible because of random simulation of
# alternative locations. Spots were simulated multiple times to ensure that
# the reported effects are not affected in size, direction or credibility.
regression_list_plot_eco %>% 
  group_by(name, Predictor) %>% 
  dplyr::summarize(mean = mean(Estimate),
                   lower = quantile(Estimate, .025),
                   upper = quantile(Estimate, .975))

# regression results
regression_plot_spatial_eco = regression_list_plot_eco %>% 
  filter(name %in% c("Social", "Roughness", "Successful", "Unsuccessful")) %>% 
  mutate(name = factor(name, levels = (c("Roughness", "Social", "Unsuccessful", "Successful")))) %>% 
  ggplot(aes(x = name, y = Estimate, group = Predictor, color = Predictor)) + 
  stat_pointinterval(position = position_dodge(.4),size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  scale_color_nejm(breaks = c("% Successful", "Avg. Catch (kg)     ")) + 
  coord_flip(ylim = c(-1, 1)) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  xlab("Feature") + 
  ylab("Posterior Estimate (95% CrI)\nLake Feature Weight ~ Predictor") + 
  labs(tag = "E", color = "Predictor")+
  ggtitle("Ecological Variation")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans"),
        legend.box.margin = margin(10, 10, 10, 10))
regression_plot_spatial_eco

################################################################################
# Part 2: Analyses Patch Leaving Model
################################################################################
# load patch leaving model fit
tmp_model = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/e_patch_leaving_model_bandwidth_fit.rds")
identifiers = readRDS(file = "utils/data/processed_data/identifiers_patch_leaving.rds")

################################################################################
# 2a. Variance Decomposition tmp Strategies
################################################################################
variance_decomposition_tmp = tmp_model %>% 
  spread_draws(sigma_id[parameter],
               sigma_trip[parameter],
               sigma_lake[parameter]) %>% 
  filter(parameter %in% c(2,3)) %>% 
  mutate(Id = sigma_id^2 / (sigma_trip^2 + sigma_lake ^2 + sigma_id ^2),
         Lake = sigma_lake^2 / (sigma_trip^2 + sigma_lake ^2 + sigma_id ^2)) %>% 
  pivot_longer(cols = c(Id, Lake), names_to = "Component", values_to = "est") %>% 
  ggplot(aes(x = as.factor(-parameter), y = est, color = Component, group = Component)) + 
  stat_pointinterval(position = position_dodge(.2),size = .6, point_size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = rev(c("Fish\nDiscovery", "Time w/o\nCatch"))) + 
  coord_flip(ylim = c(0, 1))+ 
  ylab("Posterior Estimate (95% CrI)\nProp. Expl. Variance ~ Component")+
  labs(tag = "B", color = "Component:")+
  ggtitle("Variance Decomposition Spot-Leaving")+
  xlab("Parameter")+
  scale_color_manual(values = c("green4", "lightblue3"), breaks = c("Lake", "Id"))+
  #scale_color_npg()+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, angle = 52, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        #legend.position = "none",
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_blank(),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.position = "bottom",
        legend.key.height = unit(.1, "cm"))
variance_decomposition_tmp

################################################################################
# 2b. Correlations parameters (loo analysis)
################################################################################
# loop through all lakes
if (!file.exists("utils/data/processed_data/cors_tmp.rds")){
  Cors_res = list()
  for (i in 1:10){
    
    # make two sets
    sub1 = trip_offsets %>% 
      filter(model_lake == i) %>% 
      group_by(model_id, parameter, .draw) %>% 
      dplyr::summarize(estimate = mean(trip_estimate))
    
    sub2 = trip_offsets %>% 
      filter(model_lake != i) %>% 
      group_by(model_id, parameter, .draw) %>% 
      dplyr::summarize(estimate = mean(trip_estimate))
    
    # which ids in first set
    ids = unique(sub1$model_id)
    
    # filter ids  in set 2
    sub2 = sub2 %>% filter(model_id %in% ids)
    
    # filter in other direction
    sub1 = sub1 %>% filter(model_id %in% unique(sub2$model_id))
    
    # for each parameter and each draw, compute correlation
    max_dr = max(sub2$.draw)
    cors = matrix(nrow = max_dr, ncol = 11)
    for (par in 1:11){
      
      for (dr in 1:max_dr){
        
        cors[dr, par] = cor(sub1[sub1$.draw == dr & sub1$parameter == par,"estimate"],
                            sub2[sub2$.draw == dr & sub2$parameter == par,"estimate"])
        
      }
      
    }
    
    # Save tidy dataframe
    cors = as.data.frame(cors)
    colnames(cors) = paste("parameter_", 1:11, sep = "")
    cors$draw = 1:max_dr
    cors$N1 = length(unique(sub2$model_id))
    cors$N2 = length(unique(sub1$model_id))
    cors$lake = i
    
    Cors_res[[i]] <- cors
  }
  
  # plot
  cors_tmp = do.call(rbind, Cors_res)
  saveRDS(cors_tmp, file = "utils/data/processed_data/cors_tmp.rds")
}

cors_tmp = readRDS(file = "utils/data/processed_data/cors_tmp.rds")

cors_tmp %>% 
  group_by(lake) %>% 
  dplyr::summarize(N = mean(N1)) %>% 
  dplyr::summarize(N = mean(N))
cors_tmp_plot = cors_tmp %>%
  pivot_longer(2:3) %>% 
  filter(name %in% c("parameter_2",
                     "parameter_3")) %>% 
  mutate(name = factor(name, levels = c("parameter_3", "parameter_2"))) %>% 
  ggplot(aes(x = name, y = value, group = lake)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = .6)+
  stat_pointinterval(position = position_dodge(.9), size = .6, point_size = 1, alpha = .3) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = rev(c("Fish\nDiscovery", "Time w/o\nCatch"))) + 
  coord_flip(ylim = c(-1, 1))+ 
  ylab("Posterior Estimate (95% CrI)\nCorrelation Ind. Estimates Across Lakes")+
  ggtitle("Parameter Correlation Patch Leaving")+
  xlab("Feature")+
  #scale_color_manual(values = c("green4", "lightblue3"), breaks = c("Lake", "Id"))+
  #scale_color_npg()+
  labs(tag = "B")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, angle = 50, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.6, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans"))
cors_tmp_plot

################################################################################
# 2c. Predicting Ind. level strategies from person variables
################################################################################
# note: run these scripts in succession. don't first run tmp variables and then
# spatial variables. some variables are named the same in both parts and thus this
# might cause the code not functioning properly.
id_effects = tmp_model %>% 
  spread_draws(v_id[model_participant_id, parameter])
identifiers = readRDS(file = "utils/data/processed_data/identifiers_patch_leaving.rds")

# join participant ids
identifiers$participant_id = as.character(identifiers$participant_id)
id_effects = left_join(id_effects, identifiers %>% group_by(participant_id) %>% filter(spot_id >1) %>% slice(1) %>% dplyr::select(participant_id, model_participant_id))

# join id level predictors to identifiers
catch_data <- read_csv("utils/data/raw_data/catch_data.csv")
variables = catch_data %>% 
  group_by(participant_id) %>% 
  slice(1) %>% 
  dplyr::select(sex, year_of_birth)
variables$age = 2023-variables$year_of_birth
variables$sex = as.numeric(as.factor(variables$sex))-1
variables$participant_id = as.character(variables$participant_id)
variables = left_join(identifiers, variables)

survey_data = read_csv("utils/data/raw_data/survey_data.csv")
skill = survey_data %>% dplyr::select(participant_id, angling_skill)
skill$participant_id = as.character(skill$participant_id)
dem_dat = left_join(skill, variables %>% group_by(participant_id) %>% slice(1) %>% dplyr::select(participant_id, sex, age, model_participant_id)) %>% dplyr::select(participant_id, angling_skill, sex, age, model_participant_id)

# join demographic data to original dataset
id_effects = left_join(id_effects, dem_dat)
id_effects = id_effects[!is.na(id_effects$age),]
id_effects = id_effects[!is.na(id_effects$angling_skill),]
id_effects = id_effects[!is.na(id_effects$sex),]

if (!file.exists("utils/data/processed_data/ind_regression_tmp.rds")){
  
start_time = Sys.time()
regression_list = list()
for (draw in unique(id_effects$.draw)){
  
  b = matrix(nrow = 3, ncol = 4)
  for (par_tmp in 1:3){
    
    # subset draw
    sub = id_effects[id_effects$.draw == draw & id_effects$parameter == par_tmp,]
    
    # compute X and y
    X = cbind(1, scale(sub$age)[,1], sub$sex, scale(sub$angling_skill)[,1])
    y = scale(sub$v_id)[,1]
    
    # compute normal equation
    b[par_tmp,] <- t(solve(t(X)%*%X)%*%t(X)%*%y)
    
  }
  
  # tidy results
  colnames(b) = c("Intercept", "age", "sex (ref. = m)", "skill")
  b = as.data.frame(b)
  b$parameter = 1:3
  b$draw = draw
  b$N = nrow(sub)
  
  # save in list
  regression_list[[draw]] = b
  
}
end_time = Sys.time()
end_time - start_time
regression_list = do.call(rbind, regression_list)
parameter_key = data.frame(parameter = 1:3, name = c("Baseline", "Fish\nDiscovery", "Time w/o\nCatch"))
regression_list_plot_tmp = left_join(regression_list, parameter_key)
regression_list_plot_tmp = regression_list_plot_tmp %>% pivot_longer(2:4, names_to = "Predictor", values_to = "Estimate")
# accidentally mislabelled gender. relabel here. m = 1, f = 0.
regression_list_plot_tmp = regression_list_plot_tmp %>% mutate(Predictor = ifelse(Predictor == "sex (ref. = m)", "gender\n(ref. = f)", Predictor))

# save
saveRDS(regression_list_plot_tmp, file = "utils/data/processed_data/ind_regression_tmp.rds")
}

regression_list_plot_tmp = readRDS(file = "utils/data/processed_data/ind_regression_tmp.rds")
regression_list_plot_tmp %>%
  filter(parameter %in% c(2, 3)) %>% 
  group_by(name, Predictor) %>% 
  dplyr::summarize(mean = mean(Estimate),
                   lower = quantile(Estimate, .025),
                   upper = quantile(Estimate, .975))

# regression results
regression_plot_tmp = regression_list_plot_tmp %>% 
  filter(name %in% c("Fish\nDiscovery", "Time w/o\nCatch")) %>% 
  mutate(name = factor(name, levels = c("Fish\nDiscovery", "Time w/o\nCatch"))) %>% 
  ggplot(aes(x = as.factor(-parameter), y = Estimate, group = Predictor, color = Predictor)) + 
  stat_pointinterval(position = position_dodge(.4),size = .6, point_size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = rev(c("Fish\nDiscovery", "Time w/o\nCatch"))) + 
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  scale_color_npg(breaks = c("skill", "gender\n(ref. = f)", "age")) + 
  coord_flip(ylim = c(-1.5, 1.5)) + 
  scale_y_continuous(breaks = seq(from = -1.5, to = 1.5, by = .5))+
  geom_hline(yintercept = 0, linetype = "dotted") + 
  xlab("Parameter") + 
  ylab("Posterior Estimate (95% CrI)\nInd. Parameter ~ Predictor") + 
  labs(tag = "D", color = "Predictor:")+
  ggtitle("Individual Variation (N = 72)")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_blank(),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.position = "bottom",
        legend.key.height = unit(0, "cm"),
        legend.margin = margin(2.8, 2.8, 2.8, 2.8))
regression_plot_tmp
################################################################################
# success
# average individual level catch (scaled per day to account for nesting).
rank_data = catch_data %>% 
  group_by(day, year) %>% 
  mutate(rank = rank(catch)/n(),
         catch_s = scale(catch)[, 1]) %>% 
  group_by(participant_id) %>% 
  dplyr::summarize(rank = mean(rank),
                   catch = mean(catch_s))
rank_data$scaled_rank = scale(rank_data$rank)[,1]
rank_data = rank_data %>% filter(!is.na(participant_id))
rank_data$participant_id = as.character(rank_data$participant_id)
id_effects = left_join(id_effects, rank_data)

# compute posterior regression results predicting catch from strategies
# similar results when using avg rank.
start_time = Sys.time()
b_tmp = matrix(nrow = max(ind_effects$.draw), ncol = 4)
for (draw in unique(ind_effects$.draw)){
  
  # subset draw
  sub = id_effects %>% 
    filter(.draw == draw) %>% 
    pivot_wider(names_from = parameter, values_from = v_id)
  
  # compute X and y
  X = cbind(1, scale(sub[, 12:14])[, 1:3])
  y = scale(sub$catch)[,1]
  
  # compute normal equation
  b_tmp[draw,] <- t(solve(t(X)%*%X)%*%t(X)%*%y)
  
}

# tidy results
colnames(b_tmp) = c("Intercept", "Baseline", "Fish\nDiscovery", "Time w/o\nCatch")
b_tmp = as.data.frame(b_tmp)
b_tmp$draw = 1:draw

success_prediction_tmp = b_tmp %>% 
  pivot_longer(3:4) %>% 
  ggplot(aes(x = name, y = value)) + 
  stat_pointinterval(size = .6, point_size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = (c("Fish\nDiscovery", "Time w/o\nCatch"))) + 
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  #scale_color_npg(breaks = c("skill", "gender (ref. = f)", "age")) + 
  coord_flip() + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ylab("Posterior Estimate (95% CrI)\nForaging Success ~ Ind. Feature Weight") + 
  xlab("Predictor") + 
  labs(tag = "B")+
  ggtitle("Foraging Success (N = 71)")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, angle = 55, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.2, .5),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans"),
        legend.box.margin = margin(3, 3, 3, 3))
success_prediction_tmp

################################################################################
# 2g. Predicting Lake--level variation by lake--characteristics
################################################################################
# lake characteristics (need model data)
spot_selection_data = read.csv("utils/data/raw_data/spot_selection_data.csv")
spot_selection_data = spot_selection_data %>% 
  filter(!exclude)
key = identifiers %>% group_by(day, year) %>% slice(1) %>% dplyr::select(day, year, model_lake_id)

# lake estimates
lake_estimates = tmp_model %>% 
  spread_draws(v_lake[model_lake_id, parameter])
lake_estimates = left_join(lake_estimates, key)

# 1. ratio of successful angling spots [%] (here measured as ratio, not %, but bc of scaling later 
# does not affect result)
ratio_successful = spot_selection_data %>% group_by(day, year) %>% dplyr::summarize(ratio_successful = mean(Return))
ratio_successful = left_join(lake_estimates, ratio_successful)

# 2. average amount of fish caught [kg] (here measured in g, but bc of scaling later does 
# not affect result)
avg_catch = catch_data %>% group_by(day, year) %>% dplyr::summarize(avg_catch = mean(catch))
avg_catch = left_join(ratio_successful, avg_catch)

test_data = avg_catch
# 3. posterior estimate, v_lake from lake characteristics
if(!file.exists("utils/data/processed_data/eco_regression_tmp.rds")){
results = list()
for (j in 1:3){
  
  b = matrix(ncol = 2, nrow = max(test_data$.draw))
  
  for (i in 1:max(test_data$.draw)){
    # subset draw
    sub = test_data %>% 
      filter(.draw == i & parameter == j)
    
    for (k in 1:2){
      
      # compute X and y
      X = cbind(1, scale(sub$avg_catch)[, 1],
                scale(sub$ratio_successful)[, 1])
      
      X = X[, c(1, (k+1))]
      y = scale(sub$v_lake)[,1]
      
      # compute normal equation
      b[i,k] <- (t(solve(t(X)%*%X)%*%t(X)%*%y))[, 2]
      
    }
    
  }
  
  b = as.data.frame(b)
  colnames(b) <- c("% Successful", "Avg. Catch (kg)")
  b$parameter = j
  b$draw = 1:max(test_data$.draw)
  
  results[[j]] <- b
  
}
results = do.call(rbind, results)

# make figure
parameter_key = data.frame(parameter = 1:3, name = c("Baseline", "Fish\nDiscovery", "Time w/o\nCatch"))
regression_list_plot_eco_tmp = left_join(results, parameter_key)
regression_list_plot_eco_tmp = regression_list_plot_eco_tmp %>% pivot_longer(1:2, names_to = "Predictor", values_to = "Estimate")

# save
saveRDS(regression_list_plot_eco_tmp, file = "utils/data/processed_data/eco_regression_tmp.rds")
} else {
  regression_list_plot_eco_tmp = readRDS(file = "utils/data/processed_data/eco_regression_tmp.rds")
}
regression_list_plot_eco_tmp %>% 
  filter(parameter %in% c(2, 3)) %>% 
  group_by(name, Predictor) %>% 
  dplyr::summarize(mean = mean(Estimate),
                   lower = quantile(Estimate, .025),
                   upper = quantile(Estimate, .975))

# regression results
regression_plot_tmp_eco = regression_list_plot_eco_tmp %>% 
  filter(name %in% c("Fish\nDiscovery", "Time w/o\nCatch")) %>% 
  mutate(name = factor(name, levels = c("Fish\nDiscovery", "Time w/o\nCatch"))) %>% 
  ggplot(aes(x = as.factor(-parameter), y = Estimate, group = Predictor, color = Predictor)) + 
  stat_pointinterval(position = position_dodge(.2),size = .6, point_size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = rev(c("Fish\nDiscovery", "Time w/o\nCatch"))) + 
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  scale_color_nejm(breaks = c("Avg. Catch (kg)", "% Successful")) + 
  coord_flip(ylim = c(-1, 1)) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  xlab("Parameter") + 
  ylab("Posterior Estimate (95% CrI)\nLake Parameter ~ Predictor") + 
  labs(tag = "F", color = "Predictor:")+
  ggtitle("Ecological Variation (N = 10)")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_blank(),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(0, 0, 0, 0),
        legend.position = "bottom",
        legend.key.height = unit(.1, "cm"))
regression_plot_tmp_eco
################################################################################
# make plot
library(gridExtra)
panel_int_diff = arrangeGrob(
  grobs = list(variance_decomposition_spatial,
               regression_plot_spatial,
               regression_plot_spatial_eco,
               variance_decomposition_tmp,
               regression_plot_tmp,
               regression_plot_tmp_eco),
  widths = c(1, 1, 1),
  heights = c(2, 1.85),
  layout_matrix = rbind(c(1, 2, 3),
                        c(4, 5, 6))
)

ggsave(panel_int_diff, file = "5_inter_individual_differences/output/panel_int_diff.png",
       width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
ggsave(panel_int_diff, file = "5_inter_individual_differences/output/panel_int_diff.svg",
       width = 18, height = 10, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
ggsave(panel_int_diff, file = "5_inter_individual_differences/output/panel_int_diff.pdf",
       width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
ggsave(panel_int_diff, file = "5_inter_individual_differences/output/panel_int_diff.tiff",
       width = 18, height = 10, bg = "white", dpi = 600, units = "cm")

################################################################################
# supplement
# success prediction
panel_success_prediction = arrangeGrob(
  grobs = list(success_prediction_spatial,
               success_prediction_tmp,
               icc_catch),
  widths = c(1, .03, 1),
  heights = c(1, 1),
  layout_matrix = rbind(c(1, 2, 2), 
                        c(1, NA, 3))
)
ggsave(panel_success_prediction, file = "5_inter_individual_differences/output/panel_success_prediction.png", 
       width = 18, height = 9, bg = "white", dpi = 600, units = "cm")
ggsave(panel_success_prediction, file = "5_inter_individual_differences/output/panel_success_prediction.svg", 
       width = 18, height = 9, bg = "white", dpi = 600, units = "cm", svglite::svglite)
ggsave(panel_success_prediction, file = "5_inter_individual_differences/output/panel_success_prediction.pdf", 
       width = 18, height = 9, bg = "white", dpi = 600, units = "cm")
ggsave(panel_success_prediction, file = "5_inter_individual_differences/output/panel_success_prediction.tiff", 
       width = 18, height = 9, bg = "white", dpi = 600, units = "cm")

# correlations parameters across lakes
panel_cors = arrangeGrob(
  grobs = list(cors_spatial_plot,
               cors_tmp_plot),
  widths = c(1),
  heights = c(2, 1.35),
  layout_matrix = rbind(c(1), 
                        c(2))
)
ggsave(panel_cors, file = "5_inter_individual_differences/output/panel_cors.png", 
       width = 9, height = 9, bg = "white", dpi = 600, units = "cm")
ggsave(panel_cors, file = "5_inter_individual_differences/output/panel_cors.svg", 
       width = 9, height = 9, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
ggsave(panel_cors, file = "5_inter_individual_differences/output/panel_cors.pdf", 
       width = 9, height = 9, bg = "white", dpi = 600, units = "cm")
ggsave(panel_cors, file = "5_inter_individual_differences/output/panel_cors.tiff", 
       width = 9, height = 9, bg = "white", dpi = 600, units = "cm")

################################################################################
# avgerage giving up times and lake predictions
patch_leaving_data = fread("utils/data/raw_data/patch_leaving_data.csv")
patch_leaving_data = patch_leaving_data %>% filter(!exclude)

# compute GUT
gut_data = patch_leaving_data %>%
  filter(spot_id > 1) %>% 
  filter(!exclude) %>% 
  group_by(day, year, camera_id, sequence, participant_id) %>% 
  summarise(GUT = ifelse(sum(reward)>0, max(bin)-max(bin[reward == 1]), max(bin)),
            Success = as.numeric(sum(reward)>0))

# compute lake characteristics
avg_catch = catch_data %>% group_by(day, year) %>% dplyr::summarize(avg_catch = mean(catch)/1000)
percent_successful = gut_data %>% group_by(day, year) %>% dplyr::summarize(percent_successful = mean(Success)*100)

# join to data
gut_data = left_join(gut_data, avg_catch)
gut_data = left_join(gut_data, percent_successful)
gut_data$GUT = gut_data$GUT * 10 #gut in seconds
gut_data$GUT = ifelse(gut_data$GUT == 0, 1, gut_data$GUT)

if (!file.exists( "utils/data/processed_data/model_gut.rds")){
  
  # used Gaussian model here, but lognormal (commented out) yields same results.
  model_gut = brm(data = gut_data, 
                  formula = GUT ~ Success*avg_catch + (1 + Success*avg_catch | participant_id),
                  chains = 4, 
                  cores = 4)
  #model_gut_ln = brm(data = gut_data, 
  #                 formula = bf(GUT ~ Success*avg_catch + (1 + Success*avg_catch | participant_id)),
  #                 chains = 4, 
  #                 family = lognormal(),
  #                 cores = 4)
  
  model_gut2 = brm(data = gut_data, 
                   formula = GUT ~ Success*percent_successful + (1 + Success*percent_successful | participant_id),
                   chains = 4, 
                   cores = 4)
  
  #model_gut2_ln = brm(data = gut_data, 
  #                 formula = bf(GUT ~ Success*percent_successful + (1 + Success*percent_successful | participant_id)),
  #                 chains = 4, 
  #                 family = lognormal(),
  #                 cores = 4)
  
  saveRDS(model_gut, "utils/data/processed_data/model_gut.rds")
  saveRDS(model_gut2, "utils/data/processed_data/model_gut2.rds")
}
model_gut = readRDS("utils/data/processed_data/model_gut.rds")
model_gut2 = readRDS("utils/data/processed_data/model_gut2.rds")

coef = coef(model_gut)
ranef = ranef(model_gut)
fixef = fixef(model_gut)

# make posterior predictions for each group
get_variables(model_gut)
posterior_predictions = model_gut %>% 
  spread_draws(b_Intercept,
               b_Success,
               b_avg_catch,
               `b_Success:avg_catch`) %>%
  group_by(.draw) %>% 
  slice(rep(row_number(), 400)) %>% 
  mutate(avg_catch = rep(seq(from = 0, to = 2.100, length.out = 200), 2),
         catch = rep(0:1, each = 200)) %>% 
  mutate(gut_prediction = b_Intercept + b_Success*catch + b_avg_catch*avg_catch + `b_Success:avg_catch`*catch*avg_catch)

posterior_predictions2 = model_gut2 %>% 
  spread_draws(b_Intercept,
               b_Success,
               b_percent_successful,
               `b_Success:percent_successful`) %>%
  group_by(.draw) %>% 
  slice(rep(row_number(), 400)) %>% 
  mutate(percent_successful = rep(seq(from = 0, to = 35, length.out = 200), 2),
         catch = rep(0:1, each = 200)) %>% 
  mutate(gut_prediction = b_Intercept + b_Success*catch + b_percent_successful*percent_successful + `b_Success:percent_successful`*catch*percent_successful)

# visualize effect
gut_data = gut_data %>% 
  mutate(catch = ifelse(Success == 0, "Unsuccessful Spots     ", "Successful Spots"))

n_fun <- function(x){
  return(data.frame(y = mean(x) + rnorm(1, 0, .15),
                    label = paste0(length(x),"\n\n")))
}
gut_data = gut_data %>% 
  group_by(percent_successful, catch) %>% 
  mutate(n = n())

avg_catch_prediction = ggplot() + 
  stat_lineribbon(data = posterior_predictions %>% filter(catch == 0), aes(x = avg_catch, y = gut_prediction/60), .width = c(.95), color = "darkgreen", fill = "lightgrey", alpha = .45)+
  stat_lineribbon(data = posterior_predictions %>% filter(catch == 1), aes(x = avg_catch, y = gut_prediction/60), .width = c(.95), color = "cyan3", fill = "lightgrey", alpha = .45)+
  stat_summary(data = gut_data, aes(x = avg_catch, y = GUT/60, color = catch, size = n), alpha = .4, show.legend = F) + 
  stat_summary(data = gut_data, aes(x = avg_catch, y = GUT/60, color = catch),fun.data = mean_se, size = .1) +
  stat_summary(data = gut_data, aes(x = avg_catch, y = GUT/60, group = catch), fun.data = n_fun, geom = "text", size = 4/.pt) + 
  scale_color_manual(values = c("cyan3","darkgreen")) + 
  xlab("Average Catch (kg)") + 
  coord_cartesian(ylim = c(0, 9))+
  scale_y_continuous(breaks = seq(from = 0, to = 9, by = 1))+
  ylab("Avg. Giving-up Time (min)") + 
  labs(tag = "A", color = "")+
  scale_size(range=c(0.1, 2))+
  #ggtitle("Behavioral Flexibility [N = 10]")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans"),
        legend.box.margin = margin(1, 1, 1, 1))
avg_catch_prediction

ratio_successful_prediction = ggplot() + 
  stat_lineribbon(data = posterior_predictions2 %>% filter(catch == 0), aes(x = percent_successful, y = gut_prediction/60), .width = c(.95), color = "darkgreen", fill = "lightgrey", alpha = .45)+
  stat_lineribbon(data = posterior_predictions2 %>% filter(catch == 1), aes(x =percent_successful, y = gut_prediction/60), .width = c(.95), color = "cyan3", fill = "lightgrey", alpha = .45)+
  stat_summary(data = gut_data, aes(x = percent_successful, y = GUT/60, color = catch, size = n), alpha = .4, show.legend = F) + 
  stat_summary(data = gut_data, aes(x = percent_successful, y = GUT/60, color = catch),fun.data = mean_se, size = .1) +
  stat_summary(data = gut_data, aes(x = percent_successful, y = GUT/60, group = catch), fun.data = n_fun, geom = "text", size = 4/.pt) + 
  scale_color_manual(values = c("cyan3","darkgreen")) + 
  xlab("% Successful Spots") + 
  scale_size(range=c(0.1, 2))+
  coord_cartesian(ylim = c(0, 9))+
  scale_y_continuous(breaks = seq(from = 0, to = 9, by = 1))+
  ylab("Avg. Giving-up Time (min)") + 
  labs(tag = "B", color = "")+
  #ggtitle("Behavioral Flexibility [N = 10]")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .8),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(1, 1, 1, 1))
ratio_successful_prediction

# make panel
panel_gut_edf = arrangeGrob(
  grobs = list(avg_catch_prediction,
               ratio_successful_prediction),
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2))
)
ggsave(panel_gut_edf, file = "5_inter_individual_differences/output/panel_gut_lake.png", 
       width = 18, height = 8, bg = "white", dpi = 300, units = "cm")
ggsave(panel_gut_edf, file = "5_inter_individual_differences/output/panel_gut_lake.svg", 
       width = 18, height = 6, bg = "white", dpi = 300, units = "cm")
ggsave(panel_gut_edf, file = "5_inter_individual_differences/output/panel_gut_lake.pdf", 
       width = 18, height = 6, bg = "white", dpi = 300, units = "cm")
ggsave(panel_gut_edf, file = "5_inter_individual_differences/output/panel_gut_lake.tiff", 
       width = 18, height = 6, bg = "white", dpi = 300, units = "cm")

################################################################################
# REVISION: Travel times
################################################################################
# compute average time between two subsequent successful spots for each foraging trip
spot_selection_data = fread(file = "utils/data/raw_data/spot_selection_data.csv")
spot_selection_data = spot_selection_data %>% filter(!exclude)
catch_data = fread("utils/data/raw_data/catch_data.csv")
avg_catch = catch_data %>% group_by(day, year) %>% dplyr::summarize(avg_catch = mean(catch)/1000)
percent_successful = spot_selection_data %>% group_by(day, year) %>% dplyr::summarize(percent_successful = mean(Return)*100)
travel_time = spot_selection_data %>% 
  group_by(day, year, camera_id, participant_id) %>% 
  filter(Return == 1) %>% 
  mutate(travel_time = time_start - lag(time_end)) %>% 
  dplyr::summarize(travel_time = mean(travel_time, na.rm = T)) %>% 
  filter(!is.na(travel_time))
travel_time = left_join(left_join(travel_time, percent_successful), avg_catch)

# add spatial spatial features
fitted_features = readRDS(file = "utils/data/processed_data/fitted_features.rds")
avg_features = fitted_features %>% group_by(day, year, camera_id) %>% dplyr::summarize(avg_success = mean(success_feature, na.rm = T),
                                                                                avg_loss = mean(loss_feature, na.rm = T))

# join back
avg_features$day= as.numeric(avg_features$day)
avg_features$year = as.numeric(avg_features$year)
travel_time = left_join(travel_time, avg_features)

# predict travel time from avg catch and percentage successful spots
travel_time$travel_time2 = travel_time$travel_time/60

model_travel_time1 = brm(data = travel_time, 
                formula = bf(travel_time2 ~ avg_catch + (1 + avg_catch | participant_id),
                             sigma~1),
                family = "lognormal",
                chains = 4, 
                cores = 4,
                iter = 6000,
                warmup = 4000)
model_travel_time2 = brm(data = travel_time, 
                         formula = bf(travel_time2 ~ percent_successful + (1 + percent_successful | participant_id),
                                      sigma~1),
                         family = "lognormal",
                         chains = 4, 
                         cores = 4,
                         iter = 6000,
                         warmup = 4000)

travel_time$lake_id = as.numeric(as.factor(paste(travel_time$day, travel_time$year)))
travel_time$avg_success2 = travel_time$avg_success-mean(travel_time$avg_success)
travel_time$avg_loss2 = travel_time$avg_loss-mean(travel_time$avg_loss)
model_travel_time3 = brm(data = travel_time, 
                         formula = bf(travel_time2 ~ avg_success2 + avg_loss2 + (1 + avg_success2 + avg_loss2 | participant_id) + (1 + avg_success2 + avg_loss2 | lake_id) ,
                                      sigma~1),
                         family = "lognormal",
                         chains = 4, 
                         cores = 4,
                         iter = 6000,
                         warmup = 4000)

# visualize
n_fun <- function(x){
  return(data.frame(y = mean(x) + 2*sd(x)/sqrt(length(x)) + 2,
                    label = paste0(length(x),"\n\n")))
}
travel_time = travel_time %>% 
  group_by(avg_catch) %>% 
  mutate(n = n())



get_variables(model_travel_time1)
posterior_predictions = model_travel_time1 %>% 
  spread_draws(b_Intercept,
               b_avg_catch,
               b_sigma_Intercept) %>%
  group_by(.draw) %>% 
  slice(rep(row_number(), 200)) %>% 
  mutate(avg_catch = seq(from = 0, to = 2100, length.out = 200)) %>% 
  mutate(travel_time_prediction = exp(b_Intercept + b_avg_catch*avg_catch + exp(b_sigma_Intercept)^2/2))

posterior_predictions2 = model_travel_time2 %>% 
  spread_draws(b_Intercept,
               b_percent_successful,
               b_sigma_Intercept) %>%
  group_by(.draw) %>% 
  slice(rep(row_number(), 200)) %>% 
  mutate(percent_successful = seq(from = 0, to = 35, length.out = 200)) %>% 
  mutate(travel_time_prediction = exp(b_Intercept + b_percent_successful*percent_successful + exp(b_sigma_Intercept)^2/2))

#visualize
travel_time_prediction1 = ggplot() + 
  stat_lineribbon(data = posterior_predictions, aes(x = avg_catch/1000, y = travel_time_prediction), .width = c(.95), color = "black", fill = "lightgrey", alpha = .45)+
  stat_summary(data = travel_time, aes(x = avg_catch/1000, y = travel_time2, size = n), alpha = .2, show.legend = F) + 
  stat_summary(data = travel_time, aes(x = avg_catch/1000, y = travel_time2),fun.data = mean_se, fun.args = list(mult = 2), size = .1) +
  stat_summary(data = travel_time, aes(x = avg_catch/1000, y = travel_time2), fun.data = n_fun, geom = "text", size = 4/.pt) + 
  scale_color_manual(values = c("cyan3","darkgreen")) + 
  xlab("Avg. Catch (kg)") + 
  scale_size(range=c(0.1, 2))+
  coord_cartesian(ylim = c(0, 50))+
  scale_y_continuous(breaks = seq(from = 0, to = 50, by = 10))+
  ylab("Avg. Travel Time Between Successful Spots (min)") + 
  labs(tag = "A", color = "")+
  #ggtitle("Behavioral Flexibility [N = 10]")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .8),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(1, 1, 1, 1))
travel_time_prediction1

travel_time_prediction2 = ggplot() + 
  stat_lineribbon(data = posterior_predictions2, aes(x = percent_successful, y = travel_time_prediction), .width = c(.95), color = "black", fill = "lightgrey", alpha = .45)+
  stat_summary(data = travel_time, aes(x = percent_successful, y = travel_time2, size = n), alpha = .2, show.legend = F) + 
  stat_summary(data = travel_time, aes(x = percent_successful, y = travel_time2),fun.data = mean_se, fun.args = list(mult = 2), size = .1) +
  stat_summary(data = travel_time, aes(x = percent_successful, y = travel_time2), fun.data = n_fun, geom = "text", size = 4/.pt) + 
  scale_color_manual(values = c("cyan3","darkgreen")) + 
  xlab("% Successful Spots") + 
  scale_size(range=c(0.1, 2))+
  coord_cartesian(ylim = c(0, 50))+
  scale_y_continuous(breaks = seq(from = 0, to = 50, by = 10))+
  ylab("Avg. Travel Time Between Successful Spots (min)") + 
  labs(tag = "B", color = "")+
  #ggtitle("Behavioral Flexibility [N = 10]")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .8),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        text = element_text(family = "sans"),
        legend.box.margin = margin(1, 1, 1, 1))
travel_time_prediction2

# make panel
panel_travel_time = arrangeGrob(
  grobs = list(travel_time_prediction1,
               travel_time_prediction2),
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2))
)
ggsave(panel_travel_time, file = "5_inter_individual_differences/output/panel_travel_time.png", 
       width = 18, height = 8, bg = "white", dpi = 600, units = "cm")
ggsave(panel_travel_time, file = "5_inter_individual_differences/output/panel_travel_time.svg", 
       width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
ggsave(panel_travel_time, file = "5_inter_individual_differences/output/panel_travel_time.pdf", 
       width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
ggsave(panel_travel_time, file = "5_inter_individual_differences/output/panel_travel_time.tiff", 
       width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
}
################################################################################
# END
################################################################################
