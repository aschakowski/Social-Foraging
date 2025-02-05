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
#detach("package:weights", unload=TRUE)
################################################################################
# Part 0: Load Demographic Data and Catch Results
################################################################################
# load participant data and demographics
catch_data <- read_delim("utils/data/raw_data/catch_data.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
  
# survey data
library(readr)
survey_data <- read_delim("utils/data/raw_data/survey_data.csv")
catch_data$P_ID = catch_data$participant_id
catch_data = left_join(catch_data, survey_data)

################################################################################
# Part 1: Analyses Spatial Model
################################################################################
# load spatial model fit
spatial_model = readRDS(file = "utils/data/processed_data/model_fit/spot_selection_models/f_conditional.rds")

# create key linking unique participant id to model estimates (recreate prepare stan data)
# spatial model
identifiers = readRDS(file = "utils/data/processed_data/identifiers_spot_selection.rds")
social_distances =  readRDS(file = "utils/data/processed_data/social_distance_matrix.rds")
social_distances = social_distances[(identifiers$spot_id!= 1)]
sub_ids = identifiers[identifiers$spot_id != 1, ][1:length(social_distances),]

# lake index
sub_ids$lake_index = as.numeric(as.factor(sub_ids$lake_id))
sub_ids$model_lake = sub_ids$lake_index

# trip index
sub_ids$trip_index = as.numeric(as.factor(sub_ids$trip_id))
sub_ids$model_trip = sub_ids$trip_index

# participant index
sub_ids$id_index = as.numeric(as.factor(sub_ids$participant_id))
sub_ids$model_id = sub_ids$id_index

# add participant id to data
sub_ids$day = as.numeric(sub_ids$day)
sub_ids$year = as.numeric(sub_ids$year)
catch_data = left_join(catch_data, sub_ids %>% group_by(model_trip) %>% slice(1))

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
  labs(tag = "a", color = "Component  ")+
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
  spread_draws(v_id[model_id, parameter])
trip_offsets = spatial_model %>% 
  spread_draws(v_trip[model_trip, parameter])
lake_offsets = spatial_model %>% 
  spread_draws(v_lake[model_lake, parameter])

# key
key = catch_data %>% 
  dplyr::select(model_lake, model_trip, model_id, day, year, sex, catch, year_of_birth)

# join together
trip_offsets = left_join(trip_offsets, key)
trip_offsets = left_join(trip_offsets, lake_offsets)
trip_offsets = left_join(trip_offsets, id_offsets)

# mutate
trip_offsets$trip_estimate = trip_offsets$v_lake + trip_offsets$v_id + trip_offsets$v_trip

# loop through all lakes
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
  cors$N1 = length(unique(sub2$model_id))
  cors$N2 = length(unique(sub1$model_id))
  cors$lake = i
  
  Cors_res[[i]] <- cors
}

# plot
cors_spatial = do.call(rbind, Cors_res)
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
  ggtitle("Individual Consistency")+
  xlab("Feature")+
  #scale_color_manual(values = c("green4", "lightblue3"), breaks = c("Lake", "Id"))+
  #scale_color_npg()+
  labs(tag = "a")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, angle = 50.5, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.6, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans"))
cors_spatial_plot
################################################################################
# 1c. Predicting Ind. level strategies from person variables
################################################################################
id_effects = spatial_model %>% 
  spread_draws(v_id[model_id, parameter])
dem_dat = catch_data %>% dplyr::select(model_id, age, sex, angling_skill) %>% 
  group_by(model_id) %>% 
  slice(1)

id_effects = left_join(id_effects, dem_dat)
id_effects = id_effects[!is.na(id_effects$age),]
id_effects = id_effects[!is.na(id_effects$angling_skill),]
id_effects = id_effects[!is.na(id_effects$sex),]
id_effects$sex = ifelse(id_effects$sex == "M", 1, 0)

start_time = Sys.time()
regression_list = list()
for (draw in unique(id_effects$.draw)){#unique(id_effects$.draw)
  
  b = matrix(nrow = 4, ncol = 4)
  for (par_tmp in 1:4){
    
    # subset draw
    sub = id_effects[id_effects$.draw == draw & id_effects$parameter == par_tmp,]
    #sub = id_effects %>% filter(.draw == draw & parameter == par_tmp)
    
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
regression_list_plot = regression_list_plot %>% mutate(Predictor = ifelse(Predictor == "sex (ref. = m)", "gender\n(ref. = f)", Predictor))

# regression results
regression_plot_spatial = regression_list_plot %>% 
  filter(name %in% c("Social", "Roughness", "Successful", "Unsuccessful")) %>% 
  mutate(name = factor(name, levels = c("Social", "Roughness", "Unsuccessful", "Successful"))) %>% 
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
  labs(tag = "b", color = "Predictor")+
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
rank_data = catch_data %>% 
  group_by(day, year) %>% 
  mutate(rank = rank(catch)/n(),
         catch_s = scale(catch)[, 1]) %>% 
  group_by(model_id) %>% 
  dplyr::summarize(rank = mean(rank),
            catch = mean(catch_s))
rank_data$scaled_rank = scale(rank_data$rank)[,1]
rank_data = rank_data %>% filter(!is.na(model_id))

# compute posterior regression results predicting rank from strategies
start_time = Sys.time()
b = matrix(nrow = max(ind_effects$.draw), ncol = 5)
for (draw in unique(ind_effects$.draw)){#unique(id_effects$.draw)
  
  # subset draw
  sub = ind_effects %>% 
    filter(.draw == draw) %>% 
    pivot_wider(names_from = parameter, values_from = v_id)
  
  # compute X and y
  X = cbind(1, scale(sub[, 5:8])[, 1:4])
  y = rank_data$scaled_rank
  
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
  #scale_x_discrete(labels = c("Loss", "Success", "Roughness", "Social")) + 
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  #scale_color_npg(breaks = c("skill", "gender (ref. = f)", "age")) + 
  coord_flip() + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ylab("Posterior Estimate (95% CrI)\nForaging Success ~ Ind. Feature Weight") + 
  xlab("Feature") + 
  labs(tag = "a")+
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
# 1e. Predicting Success From trip level weights
################################################################################
# catch results
catch_res = catch_data %>% 
  group_by(model_trip) %>% 
  slice(1) %>% 
  dplyr::select(model_trip, day, year, catch, model_id)##
################################################################################
# 1f. ICC catch success
################################################################################
catch_res$lake = paste(catch_res$day, catch_res$year)
mod_catch = brm(data = catch_res, 
                formula = catch/1000 ~ (1|model_id) + (1|lake))

icc_catch = mod_catch %>% 
  spread_draws(sd_lake__Intercept,
               sd_model_id__Intercept,
               sigma) %>% 
  mutate(var_lake = sd_lake__Intercept^2 / (sd_lake__Intercept^2 + sd_model_id__Intercept^2 + sigma^2),
         var_id = sd_model_id__Intercept^2 / (sd_lake__Intercept^2 + sd_model_id__Intercept^2 + sigma^2)) %>% 
  pivot_longer(7:8) %>% 
  ggplot(aes(x = name, y = value)) + 
  stat_pointinterval(size = .6, point_size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = c("Individual", "Lake")) + 
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  #scale_color_npg(breaks = c("skill", "gender (ref. = f)", "age")) + 
  coord_flip() + 
  labs(tag = "c")+
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

key = catch_data %>% group_by(model_lake) %>% slice(1) %>% dplyr::select(model_lake, day, year)
key = key[!is.na(key$model_lake),]

# lake estimates
lake_estimates = spatial_model %>% 
  spread_draws(v_lake[model_lake, parameter])
lake_estimates = left_join(lake_estimates, key)

# 1. ratio of successful angling spots [%]
ratio_successful = spot_selection_data %>% group_by(day, year) %>% summarize(ratio_successful = mean(Return))
ratio_successful = left_join(lake_estimates, ratio_successful)

# 2. average amount of fish caught [kg]
avg_catch = catch_data %>% group_by(day, year) %>% summarize(avg_catch = mean(catch))
avg_catch = left_join(ratio_successful, avg_catch)

# 3. variation in depth profile
fitted_features=readRDS(file = "utils/data/processed_data/fitted_features.rds")
var_roughness = fitted_features %>% 
  group_by(model_lake = lake) %>% 
  summarize(var = sd(roughness_feature))
test_data = left_join(avg_catch, var_roughness)  


# 3. posterior estimate, v_lake from lake characteristics
results = list()
for (j in 1:4){
  
  b = matrix(ncol = 3, nrow = max(test_data$.draw))
  
  for (i in 1:max(test_data$.draw)){
    # subset draw
    sub = test_data %>% 
      filter(.draw == i & parameter == j)
    
    for (k in 1:3){
      
      # compute X and y
      X = cbind(1, scale(sub$avg_catch)[, 1],
                scale(sub$ratio_successful)[, 1],
                scale(sub$var)[,1])
      
      X = X[, c(1, (k+1))]
      y = scale(sub$v_lake)[,1]
      
      # compute normal equation
      b[i,k] <- (t(solve(t(X)%*%X)%*%t(X)%*%y))[, 2]
      
    }
    
  }
  
  b = as.data.frame(b)
  colnames(b) <- c("% Successful", "Avg. Catch (kg)     ", "Var. Roughness")
  b$parameter = j
  b$draw = 1:max(test_data$.draw)
  
  results[[j]] <- b
  
}
results = do.call(rbind, results)

# make figure
parameter_key = data.frame(parameter = 1:8, name = c("Social", "Roughness", "Successful", "Unsuccessful"))
regression_list_plot_eco = left_join(results, parameter_key)
regression_list_plot_eco = regression_list_plot_eco %>% pivot_longer(1:2, names_to = "Predictor", values_to = "Estimate")

# regression results
regression_plot_spatial_eco = regression_list_plot_eco %>% 
  filter(name %in% c("Social", "Roughness", "Successful", "Unsuccessful")) %>% 
  mutate(name = factor(name, levels = c("Social", "Roughness", "Unsuccessful", "Successful"))) %>% 
  ggplot(aes(x = name, y = Estimate, group = Predictor, color = Predictor)) + 
  stat_pointinterval(position = position_dodge(.4),size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  scale_color_nejm(breaks = c("% Successful", "Avg. Catch (kg)     ")) + 
  coord_flip(ylim = c(-1, 1)) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  xlab("Feature") + 
  ylab("Posterior Estimate (95% CrI)\nLake Feature Weight ~ Predictor") + 
  labs(tag = "c", color = "Predictor")+
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

################################################################################
# Part 2: Analyses Patch Leaving Model
################################################################################
# load patch leaving model fit
tmp_model = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/e_patch_leaving_model_spatial_features_fit.rds")

# create key linking unique participant id to model estimates (recreate prepare stan data)
# spatial model
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
  scale_x_discrete(labels = rev(c("Patch\nDiscovery", "Time w/o\nCatch"))) + 
  coord_flip(ylim = c(0, 1))+ 
  ylab("Posterior Estimate (95% CrI)\nProp. Expl. Variance ~ Component")+
  labs(tag = "d", color = "Component:")+
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

################################################################################
# 2b. Correlations parameters (loo analysis)
################################################################################
# load id, trip, lake offsets
id_offsets = tmp_model %>% 
  spread_draws(v_id[model_id, parameter])
trip_offsets = tmp_model %>% 
  spread_draws(v_trip[model_trip, parameter])
lake_offsets = tmp_model %>% 
  spread_draws(v_lake[model_lake, parameter])

# key
key = catch_data %>% 
  dplyr::select(model_lake, model_trip, model_id, day, year, sex, catch, year_of_birth)

# join together
trip_offsets = left_join(trip_offsets, key)
trip_offsets = left_join(trip_offsets, lake_offsets)
trip_offsets = left_join(trip_offsets, id_offsets)

# mutate
trip_offsets$trip_estimate = trip_offsets$v_lake + trip_offsets$v_id + trip_offsets$v_trip

# loop through all lakes
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

cors_tmp %>% 
  group_by(lake) %>% 
  dplyr::summarize(N = mean(N1)) %>% 
  dplyr::summarize(N = mean(N))
cors_tmp_plot = cors_tmp %>%
  pivot_longer(2:3) %>% 
  filter(name %in% c("parameter_2",
                     "parameter_3")) %>% 
  ggplot(aes(x = name, y = value, group = lake)) + 
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey", size = .6)+
  stat_pointinterval(position = position_dodge(.9), size = .6, point_size = 1, alpha = .3) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = (c("Patch\nDiscovery", "Time w/o\nCatch"))) + 
  coord_flip(ylim = c(-1, 1))+ 
  ylab("Posterior Estimate (95% CrI)\nCorrelation Ind. Estimates Across Lakes")+
  ggtitle("Parameter Correlation Patch Leaving")+
  xlab("Feature")+
  #scale_color_manual(values = c("green4", "lightblue3"), breaks = c("Lake", "Id"))+
  #scale_color_npg()+
  labs(tag = "b")+
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
id_effects = tmp_model %>% 
  spread_draws(v_id[model_id, parameter])
dem_dat = catch_data %>% dplyr::select(model_id, age, sex, angling_skill) %>% 
  group_by(model_id) %>% 
  slice(1)

id_effects = left_join(id_effects, dem_dat)
id_effects = id_effects[!is.na(id_effects$age),]
id_effects = id_effects[!is.na(id_effects$angling_skill),]
id_effects = id_effects[!is.na(id_effects$sex),]
id_effects$sex = ifelse(id_effects$sex == "M", 1, 0)

start_time = Sys.time()
regression_list = list()
for (draw in unique(id_effects$.draw)){#unique(id_effects$.draw)
  
  b = matrix(nrow = 11, ncol = 4)
  for (par_tmp in 1:11){
    
    # subset draw
    sub = id_effects[id_effects$.draw == draw & id_effects$parameter == par_tmp,]
    #sub = id_effects %>% filter(.draw == draw & parameter == par_tmp)
    
    # compute X and y
    X = cbind(1, scale(sub$age)[,1], sub$sex, scale(sub$angling_skill)[,1])
    y = scale(sub$v_id)[,1]
    
    # compute normal equation
    b[par_tmp,] <- t(solve(t(X)%*%X)%*%t(X)%*%y)
    
  }
  
  # tidy results
  colnames(b) = c("Intercept", "age", "sex (ref. = m)", "skill")
  b = as.data.frame(b)
  b$parameter = 1:11
  b$draw = draw
  b$N = nrow(sub)
  
  # save in list
  regression_list[[draw]] = b
  
}
end_time = Sys.time()
end_time - start_time
regression_list = do.call(rbind, regression_list)
parameter_key = data.frame(parameter = 1:11, name = c("Baseline", "Patch\nDiscovery", "Time w/o\nCatch", "Waiting  Time\nLocal", "Waiting Time\nGlobal", "Success", "Loss", "Roughness", "Social", "Updating\nGlobal", "Updating\nLocal"))
regression_list_plot_tmp = left_join(regression_list, parameter_key)
regression_list_plot_tmp = regression_list_plot_tmp %>% pivot_longer(2:4, names_to = "Predictor", values_to = "Estimate")
regression_list_plot_tmp = regression_list_plot_tmp %>% mutate(Predictor = ifelse(Predictor == "sex (ref. = m)", "gender\n(ref. = f)", Predictor))

# regression results
regression_plot_tmp = regression_list_plot_tmp %>% 
  filter(name %in% c("Patch\nDiscovery", "Time w/o\nCatch")) %>% 
  mutate(name = factor(name, levels = c("Patch\nDiscovery", "Time w/o\nCatch"))) %>% 
  ggplot(aes(x = as.factor(-parameter), y = Estimate, group = Predictor, color = Predictor)) + 
  stat_pointinterval(position = position_dodge(.4),size = .6, point_size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = rev(c("Patch\nDiscovery", "Time w/o\nCatch"))) + 
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  scale_color_npg(breaks = c("skill", "gender\n(ref. = f)", "age")) + 
  coord_flip(ylim = c(-1.5, 1.5)) + 
  scale_y_continuous(breaks = seq(from = -1.5, to = 1.5, by = .5))+
  geom_hline(yintercept = 0, linetype = "dotted") + 
  xlab("Parameter") + 
  ylab("Posterior Estimate (95% CrI)\nInd. Parameter ~ Predictor") + 
  labs(tag = "e", color = "Predictor:")+
  ggtitle("Individual Variation (N = 71)")+
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
# individual level parameters
ind_effects = tmp_model %>% 
  spread_draws(v_id[model_id, parameter])

# average individual level rank
rank_data = catch_data %>% 
  group_by(day, year) %>% 
  mutate(rank = rank(catch)/n(),
         catch_s = scale(catch)[, 1]) %>% 
  group_by(model_id) %>% 
  dplyr::summarize(rank = mean(rank),
            catch = mean(catch_s))
rank_data$scaled_rank = scale(rank_data$rank)[,1]
rank_data = rank_data %>% filter(!is.na(model_id))

# compute posterior regression results predicting rank from strategies
start_time = Sys.time()
b_tmp = matrix(nrow = max(ind_effects$.draw), ncol = 4)
for (draw in unique(ind_effects$.draw)){#unique(id_effects$.draw)
  
  # subset draw
  sub = ind_effects %>% 
    filter(.draw == draw) %>% 
    pivot_wider(names_from = parameter, values_from = v_id)
  
  # compute X and y
  X = cbind(1, scale(sub[, 5:7])[, 1:3])
  y = rank_data$scaled_rank
  
  # compute normal equation
  b_tmp[draw,] <- t(solve(t(X)%*%X)%*%t(X)%*%y)
  
}

# tidy results
colnames(b_tmp) = c("Intercept", "Baseline", "Patch\nDiscovery", "Time w/o\nCatch")
b_tmp = as.data.frame(b_tmp)
b_tmp$draw = 1:draw

success_prediction_tmp = b_tmp %>% 
  pivot_longer(3:4) %>% 
  ggplot(aes(x = name, y = value)) + 
  stat_pointinterval(size = .6, point_size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = c("Patch\nDiscovery", "Time w/o\nCatch")) + 
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  #scale_color_npg(breaks = c("skill", "gender (ref. = f)", "age")) + 
  coord_flip() + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  ylab("Posterior Estimate (95% CrI)\nForaging Success ~ Ind. Cue Weight") + 
  xlab("Predictor") + 
  labs(tag = "b")+
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

key = catch_data %>% group_by(model_lake) %>% slice(1) %>% dplyr::select(model_lake, day, year)
key = key[!is.na(key$model_lake),]

# lake estimates
lake_estimates = tmp_model %>% 
  spread_draws(v_lake[model_lake, parameter])
lake_estimates = left_join(lake_estimates, key)

# 1. ratio of successful angling spots [%]
ratio_successful = spot_selection_data %>% group_by(day, year) %>% summarize(ratio_successful = mean(Return))
ratio_successful = left_join(lake_estimates, ratio_successful)

# 2. average amount of fish caught [kg]
avg_catch = catch_data %>% group_by(day, year) %>% summarize(avg_catch = mean(catch))
avg_catch = left_join(ratio_successful, avg_catch)

# 3. variation in depth profile
fitted_features=readRDS(file = "utils/data/processed_data/fitted_features.rds")
var_roughness = fitted_features %>% 
  group_by(model_lake = lake) %>% 
  summarize(var = sd(roughness_feature))
test_data = left_join(avg_catch, var_roughness)  

# 3. posterior estimate, v_lake from lake characteristics
results = list()
for (j in 1:11){
  
  b = matrix(ncol = 3, nrow = max(test_data$.draw))
  
  for (i in 1:max(test_data$.draw)){
    # subset draw
    sub = test_data %>% 
      filter(.draw == i & parameter == j)
    
    for (k in 1:3){
      
      # compute X and y
      X = cbind(1, scale(sub$avg_catch)[, 1],
                scale(sub$ratio_successful)[, 1],
                scale(sub$var)[,1])
      
      X = X[, c(1, (k+1))]
      y = scale(sub$v_lake)[,1]
      
      # compute normal equation
      b[i,k] <- (t(solve(t(X)%*%X)%*%t(X)%*%y))[, 2]
      
    }
    
  }
  
  b = as.data.frame(b)
  colnames(b) <- c("% Successful", "Avg. Catch (kg)", "Var. Roughness")
  b$parameter = j
  b$draw = 1:max(test_data$.draw)
  
  results[[j]] <- b
  
}
results = do.call(rbind, results)
#results %>% 
#  ggplot(aes(x = parameter, y = `Avg. Catch [kg]`)) + 
#  stat_pointinterval()
#results %>% 
#  ggplot(aes(x = parameter, y = `Ratio Successful [%]`)) + 
#  stat_pointinterval()
#results %>% 
#  ggplot(aes(x = parameter, y = `Var. Roughness`)) + 
#  stat_pointinterval()

# make figure
parameter_key = data.frame(parameter = 1:11, name = c("Baseline", "Patch\nDiscovery", "Time w/o\nCatch", "Waiting Time\nLocal", "Waiting Time\nGlobal", "Success", "Loss", "Roughness", "Social", "Updating\nLocal", "Updating\nGlobal"))
regression_list_plot_eco_tmp = left_join(results, parameter_key)
regression_list_plot_eco_tmp = regression_list_plot_eco_tmp %>% pivot_longer(1:2, names_to = "Predictor", values_to = "Estimate")

# regression results
regression_plot_tmp_eco = regression_list_plot_eco_tmp %>% 
  filter(name %in% c("Patch\nDiscovery", "Time w/o\nCatch")) %>% 
  mutate(name = factor(name, levels = c("Patch\nDiscovery", "Time w/o\nCatch"))) %>% 
  ggplot(aes(x = as.factor(-parameter), y = Estimate, group = Predictor, color = Predictor)) + 
  stat_pointinterval(position = position_dodge(.2),size = .6, point_size = 1) + 
  scale_size_continuous(range = c(3, 10))+
  scale_x_discrete(labels = rev(c("Patch\nDiscovery", "Time w/o\nCatch"))) + 
  #scale_color_manual(values = c("green4", "lightblue3", "cyan4"), breaks = c("skill", "gender (ref. = f)", "age"))+
  scale_color_nejm(breaks = c("Avg. Catch (kg)", "% Successful")) + 
  coord_flip(ylim = c(-1, 1)) + 
  geom_hline(yintercept = 0, linetype = "dotted") + 
  xlab("Parameter") + 
  ylab("Posterior Estimate (95% CrI)\nLake Parameter ~ Predictor") + 
  labs(tag = "f", color = "Predictor:")+
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
       width = 18, height = 10, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_int_diff, file = "5_inter_individual_differences/output/panel_int_diff.svg",
       width = 18, height = 10, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_int_diff, file = "5_inter_individual_differences/output/panel_int_diff.pdf",
       width = 18, height = 10, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_int_diff, file = "5_inter_individual_differences/output/panel_int_diff.tiff",
       width = 18, height = 10, bg = "white", dpi = 1200, units = "cm")

################################################################################
# supplement
# success prediction
panel_success_prediction = arrangeGrob(
  grobs = list(success_prediction_spatial,
               success_prediction_tmp,
               icc_catch),
  widths = c(1, 1),
  heights = c(1, 1),
  layout_matrix = rbind(c(1, 2), 
                        c(1, 3))
)
ggsave(panel_success_prediction, file = "5_inter_individual_differences/output/panel_success_prediction.png", 
       width = 18, height = 9, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_success_prediction, file = "5_inter_individual_differences/output/panel_success_prediction.svg", 
       width = 18, height = 9, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_success_prediction, file = "5_inter_individual_differences/output/panel_success_prediction.pdf", 
       width = 18, height = 9, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_success_prediction, file = "5_inter_individual_differences/output/panel_success_prediction.tiff", 
       width = 18, height = 9, bg = "white", dpi = 1200, units = "cm")

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
       width = 9, height = 9, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_cors, file = "5_inter_individual_differences/output/panel_cors.svg", 
       width = 9, height = 9, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_cors, file = "5_inter_individual_differences/output/panel_cors.pdf", 
       width = 9, height = 9, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_cors, file = "5_inter_individual_differences/output/panel_cors.tiff", 
       width = 9, height = 9, bg = "white", dpi = 1200, units = "cm")

################################################################################
# avgerage giving up times and lake predictions
gut = readRDS(file = "utils/data/processed_data/patch_leaving_simulations/a_sim_baseline.rds")
gut = do.call(rbind, gut)

gut_data = gut %>%
  filter(iter == 1) %>% 
  mutate(catch = as.numeric(catch > 0))
lake_data = left_join(avg_catch %>% group_by(lake = model_lake) %>% slice(1) %>% dplyr::select(lake, ratio_successful, avg_catch), gut_data)

model_gut = brm(data = lake_data, 
                formula = gut_data ~ catch*avg_catch + (1 + catch*avg_catch | id),
                chains = 4, 
                cores = 4)
model_gut2 = brm(data = lake_data, 
                 formula = gut_data ~ catch*ratio_successful + (1 + catch*ratio_successful | id),
                 chains = 4, 
                 cores = 4)

coef = coef(model_gut)
ranef = ranef(model_gut)
fixef = fixef(model_gut)

# make posterior predictions for each group
get_variables(model_gut)
posterior_predictions = model_gut %>% 
  spread_draws(b_Intercept,
               b_catch,
               b_avg_catch,
               `b_catch:avg_catch`) %>%
  group_by(.draw) %>% 
  slice(rep(row_number(), 400)) %>% 
  mutate(avg_catch = rep(seq(from = 0, to = 2100, length.out = 200), 2),
         catch = rep(0:1, each = 200)) %>% 
  mutate(gut_prediction = b_Intercept + b_catch*catch + b_avg_catch*avg_catch + `b_catch:avg_catch`*catch*avg_catch)

posterior_predictions2 = model_gut2 %>% 
  spread_draws(b_Intercept,
               b_catch,
               b_ratio_successful,
               `b_catch:ratio_successful`) %>%
  group_by(.draw) %>% 
  slice(rep(row_number(), 400)) %>% 
  mutate(ratio_successful = rep(seq(from = 0, to = .35, length.out = 200), 2),
         catch = rep(0:1, each = 200)) %>% 
  mutate(gut_prediction = b_Intercept + b_catch*catch + b_ratio_successful*ratio_successful + `b_catch:ratio_successful`*catch*ratio_successful)

# visualize effect
lake_data = lake_data %>% 
  mutate(catch = ifelse(catch == 0, "Unsuccessful Spots     ", "Successful Spots"))

n_fun <- function(x){
  return(data.frame(y = mean(x) + rnorm(1, 0, .15),
                    label = paste0(length(x),"\n\n")))
}
lake_data = lake_data %>% 
  group_by(ratio_successful, catch) %>% 
  mutate(n = n())

avg_catch_prediction = ggplot() + 
  stat_lineribbon(data = posterior_predictions %>% filter(catch == 0), aes(x = avg_catch/1000, y = gut_prediction*10/60), .width = c(.95), color = "darkgreen", fill = "lightgrey", alpha = .45)+
  stat_lineribbon(data = posterior_predictions %>% filter(catch == 1), aes(x = avg_catch/1000, y = gut_prediction*10/60), .width = c(.95), color = "cyan3", fill = "lightgrey", alpha = .45)+
  stat_summary(data = lake_data, aes(x = avg_catch/1000, y = gut_data*10/60, color = catch, size = n), alpha = .4, show.legend = F) + 
  stat_summary(data = lake_data, aes(x = avg_catch/1000, y = gut_data*10/60, color = catch),fun.args = list(mult = 2), size = .1) +
  stat_summary(data = lake_data, aes(x = avg_catch/1000, y = gut_data*10/60, group = catch), fun.data = n_fun, geom = "text", size = 4/.pt) + 
  scale_color_manual(values = c("cyan3","darkgreen")) + 
  xlab("Average Catch (kg)") + 
  coord_cartesian(ylim = c(0, 9))+
  scale_y_continuous(breaks = seq(from = 0, to = 9, by = 1))+
  ylab("Avg. Giving up Time (min)") + 
  labs(tag = "a", color = "")+
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
  stat_lineribbon(data = posterior_predictions2 %>% filter(catch == 0), aes(x = ratio_successful*100, y = gut_prediction*10/60), .width = c(.95), color = "darkgreen", fill = "lightgrey", alpha = .45)+
  stat_lineribbon(data = posterior_predictions2 %>% filter(catch == 1), aes(x = ratio_successful*100, y = gut_prediction*10/60), .width = c(.95), color = "cyan3", fill = "lightgrey", alpha = .45)+
  stat_summary(data = lake_data, aes(x = ratio_successful*100, y = gut_data*10/60, color = catch, size = n), alpha = .4, show.legend = F) + 
  stat_summary(data = lake_data, aes(x = ratio_successful*100, y = gut_data*10/60, color = catch),fun.args = list(mult = 2), size = .1) +
  stat_summary(data = lake_data, aes(x = ratio_successful*100, y = gut_data*10/60, group = catch), fun.data = n_fun, geom = "text", size = 4/.pt) + 
  scale_color_manual(values = c("cyan3","darkgreen")) + 
  xlab("% Successful Spots") + 
  scale_size(range=c(0.1, 2))+
  coord_cartesian(ylim = c(0, 9))+
  scale_y_continuous(breaks = seq(from = 0, to = 9, by = 1))+
  ylab("Avg. Giving up Time (min)") + 
  labs(tag = "b", color = "")+
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
       width = 18, height = 8, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_gut_edf, file = "5_inter_individual_differences/output/panel_gut_lake.svg", 
       width = 18, height = 6, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_gut_edf, file = "5_inter_individual_differences/output/panel_gut_lake.pdf", 
       width = 18, height = 6, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_gut_edf, file = "5_inter_individual_differences/output/panel_gut_lake.tiff", 
       width = 18, height = 6, bg = "white", dpi = 1200, units = "cm")

}
################################################################################
# END
################################################################################