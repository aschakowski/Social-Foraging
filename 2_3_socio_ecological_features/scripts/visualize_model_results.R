################################################################################
#
# Title: Visualize model results
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 02/04/2024 (DD/MM/YYYY)
#
################################################################################
visualize_spot_selection_model <- function(){
  
  ################################################################################
  # load data
  ################################################################################
  # load spatial features
  fitted_features = readRDS(file = "utils/data/processed_data/fitted_features.rds")
  
  # load spatial model fit
  spot_selection_model_fit = readRDS(file = "utils/data/processed_data/model_fit/spot_selection_models/f_conditional.rds")
  print(spot_selection_model_fit, max_rows = 20)
  spot_selection_model_fit$print(NULL,~quantile(.x, probs = c(0.025, 0.975)),max_rows = 20)
  
  # spot data containing all angling spots
  spot_selection_data = fread("utils/data/raw_data/spot_selection_data.csv")
  spot_selection_data = spot_selection_data %>% 
    filter(!exclude)
  
  ################################################################################
  # Visualize feature weights
  ################################################################################
  lake_effects = spot_selection_model_fit %>% 
    spread_draws(weights[feature], v_lake[lake, feature]) %>% 
    mutate(lake_effect = weights + v_lake)
  lake_effects$feature = lake_effects$feature %>% recode("Social", "Roughness", "Successful","Unsuccessful", "Distance","Distance x Time", "Social x Unsuccessful","Social x Successful")
  
  # compute standardized predictors
  # by weighing unstand. estimates with avg sd of features encountered on each trip on each lake
  scale_features = fitted_features %>% 
    unnest() %>% 
    group_by(trip, lake) %>% 
    mutate(time_since_start = as.numeric(time_since_start) / max(as.numeric(time_since_start))) %>% 
    summarize(sd_social = sd(social_feature),
              sd_roughness = sd(roughness_feature),
              sd_success = sd(success_feature),
              sd_loss = sd(loss_feature),
              sd_locality = sd(locality_feature),
              sd_locality_time = sd(locality_feature) * sd(time_since_start),#not best way to standardize time, but not reported in paper
              sd_social_loss = sd(social_feature) * sd(loss_feature),
              sd_social_success = sd(social_feature) * sd(success_feature)) %>% 
    ungroup() %>% 
    group_by(lake) %>% 
    summarize(sd_social = mean(sd_social),
              sd_roughness = mean(sd_roughness),
              sd_success = mean(sd_success),
              sd_loss = mean(sd_loss),
              sd_locality = mean(sd_locality),
              sd_locality_time = mean(sd_locality_time),
              sd_social_loss = mean(sd_social_loss),
              sd_social_success = mean(sd_social_success)) %>% 
    unnest()
  
  # compute standardized predictors
  lake_effects = left_join(lake_effects, scale_features)
  
  lake_effects = lake_effects %>% 
    mutate(standardized_weight = ifelse(feature == "Social", lake_effect*sd_social, 
                                        ifelse(feature == "Roughness", lake_effect*sd_roughness, 
                                               ifelse(feature == "Successful", lake_effect*sd_success, 
                                                      ifelse(feature == "Unsuccessful", lake_effect*sd_loss,
                                                             ifelse(feature == "Social x Unsuccessful", lake_effect*sd_social_loss,
                                                                    ifelse(feature == "Social x Successful", lake_effect*sd_social_success,
                                                                           ifelse(feature == "Distance x Time", lake_effect*sd_locality_time, 
                                                                                  ifelse(feature == "Distance", lake_effect*sd_locality, 0)))))))),
           standardized_main_effect = ifelse(feature == "Social", (weights)*mean(scale_features$sd_social), 
                                             ifelse(feature == "Roughness", (weights)*mean(scale_features$sd_roughness), 
                                                    ifelse(feature == "Successful", (weights)*mean(scale_features$sd_success), 
                                                           ifelse(feature == "Unsuccessful", weights*mean(scale_features$sd_loss),
                                                                  ifelse(feature == "Social x Unsuccessful", weights * mean(scale_features$sd_social_loss),
                                                                         ifelse(feature == "Social x Successful", weights * mean(scale_features$sd_social_success),
                                                                                ifelse(feature == "Distance x Time", weights*mean(scale_features$sd_locality_time), 
                                                                                       ifelse(feature == "Distance", weights*mean(scale_features$sd_locality), 0)))))))))
  # visualize estimates
  plot_dat1 = lake_effects %>% 
    filter(!feature %in% c("Distance","Distance x Time", "Social x Unsuccessful", "Social x Successful")) %>% 
    mutate(name = factor(feature, levels=c("Roughness", "Social", "Unsuccessful", "Successful")))
  model_results = 
    ggplot(plot_dat1,aes(x = standardized_weight, y = name, group = as.factor(lake))) + 
    geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
    # .95
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .975), ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "#BAE4B3")+
    #.66
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .83), ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .5)-.005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Successful"], .5)+.005,ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "black")+
    
    # .5
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .975), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .83),ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .5)-.005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social"], .5)+.005,ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Roughness"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Roughness"], .975), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Roughness"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Roughness"], .83), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Roughness"], .5)-.005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Roughness"], .5)+.005,ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .975), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .83), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#31A354")+ 
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .5)-.005, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Unsuccessful"], .5)+.005, ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "black")+
    ylab("Feature")+
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0) +
    
    #stat_pointinterval(color = "black", alpha = .15) +
    labs(tag = "E", color = "value")+
    ggtitle("Feature Weights")+
    xlab("Posterior Est. (95%CrI)\nSpot Selection ~ Feature")+
    scale_x_continuous(breaks = seq(from = -3, to = 6, by = 1))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5, angle = 55, hjust = .5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  model_results
  
  ################################################################################
  # Visualize interaction effects
  ################################################################################
  model_results_interaction = lake_effects %>% 
    filter(feature %in% c("Social x Unsuccessful", "Social x Successful")) %>% 
    mutate(name = factor(feature, levels=c("Social x Unsuccessful", "Social x Successful"))) %>%
    ggplot(aes(x = standardized_weight, y = name, group = as.factor(lake))) + 
    geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
    # .95
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Unsuccessful"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Unsuccessful"], .975), ymin = .65, ymax = 1.35,
             alpha = 1,fill = "#BAE4B3")+
    #.66
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Unsuccessful"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Unsuccessful"], .83), ymin = .65, ymax = 1.35,
             alpha = 1,fill = "#31A354")+
    # .5
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Unsuccessful"], .5)-.001, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Unsuccessful"], .5)+.001, ymin = .65, ymax = 1.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Successful"], .025), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Successful"], .975), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Successful"], .17), xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Successful"], .83),  ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Successful"], .5)-.001, xmax =  quantile(lake_effects$standardized_main_effect[lake_effects$feature == "Social x Successful"], .5)+.001,  ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "black")+
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0)  +
    
    
    ylab("Feature")+
    labs(tag = "G", color = "value")+
    ggtitle("Conditional Weights")+
    xlab("Posterior Est. (95%CrI)\nSpot Selection ~ Feature")+
    scale_x_continuous(breaks = seq(from = -2, to = 2, by = .2))+
    scale_y_discrete(labels = c("Unsuccessful x\nSocial", "Successful x\nSocial"))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5, angle = 90, hjust = .5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  model_results_interaction
  
  ################################################################################
  # Visualize bandwidth estimates
  ################################################################################
  # extract kernel bandwidths for each lake
  bandwidths = spot_selection_model_fit %>% 
    spread_draws(lambdas[feature]) %>% 
    group_by(feature) %>% 
    summarize(h = mean(lambdas))
  bandwidths_offset = spot_selection_model_fit %>% 
    spread_draws(v_lake[lake, feature]) %>% 
    filter(feature > 8) %>% 
    group_by(lake, feature) %>% 
    summarize(offset = mean(v_lake)) %>% 
    mutate(feature = feature - 8)
  bandwidths = left_join(bandwidths, bandwidths_offset)
  bandwidths = bandwidths %>% 
    mutate(bandwidth = exp(h + offset))
  
  # full posterior for plot
  bandwidths_full = spot_selection_model_fit %>% 
    spread_draws(lambdas[feature])
  bandwidths_offset_full = spot_selection_model_fit %>% 
    spread_draws(v_lake[lake, feature]) %>%  
    mutate(feature = feature - 8)
  bandwidths_offset_full = left_join(bandwidths_full, bandwidths_offset_full)
  bandwidths_offset_full = bandwidths_offset_full %>% 
    mutate(bandwidth = lambdas + v_lake)
  bandwidths_offset_full$feature = bandwidths_offset_full$feature %>% recode("Social\nInformation","Personal\nInformation")
  bandwidths_full$feature = bandwidths_full$feature %>% recode("Social\nInformation","Personal\nInformation")
  
  bandwidth_results = bandwidths_offset_full %>% 
    filter(feature %in% c("Personal\nInformation", "Social\nInformation")) %>% 
    #filter(exp(bandwidth)<200) %>% 
    mutate(name = factor(feature, levels=c("Social\nInformation","Personal\nInformation"))) %>%
    ggplot(aes(x = exp(bandwidth), y = name, group = as.factor(lake))) + 
    #geom_vline(xintercept = 0, linetype = "dashed") + 
    # .95
    annotate("rect", xmin = quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Personal\nInformation"]), .025), xmax =  quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Personal\nInformation"]), .975), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#BAE4B3")+
    #.66
    annotate("rect", xmin = quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Personal\nInformation"]), .17), xmax =  quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Personal\nInformation"]), .83), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#31A354")+
    # .5
    annotate("rect", xmin = quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Personal\nInformation"]), .5)-.05, xmax =  quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Personal\nInformation"]), .5)+.05, ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Social\nInformation"]), .025), xmax =  quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Social\nInformation"]), .975), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#BAE4B3")+
    annotate("rect", xmin = quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Social\nInformation"]), .17), xmax =  quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Social\nInformation"]), .83), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Social\nInformation"]), .5)-.05, xmax =  quantile(exp(bandwidths_full$lambdas[bandwidths_full$feature == "Social\nInformation"]), .5)+.05,  ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "black")+
    
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0) +
    
    
    
    ylab("Feature")+
    labs(tag = "H", color = "value")+
    ggtitle("Generalization")+
    coord_cartesian(xlim = c(0, 40))+
    xlab("Posterior Est. (95%CrI)\nBandwidth")+
    scale_x_continuous(breaks = seq(from = 0, to = 50, by = 10))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5, angle = 90, hjust = .5),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  bandwidth_results
  
  ################################################################################
  # Visualize adaptive cues
  ################################################################################
  # evaluate success as function of feature values at chosen locations
  subset = fitted_features %>% filter(chosen == 1)
  
  # scale features on each lake (ests mean: given a change in sd in one feature
  # value on one lake, how much more likely is a participant to catch fish there)
  # scaling does not affect results.
  subset = subset %>% 
    group_by(day, year) %>% 
    mutate(success_scaled = scale(success_feature)[,1],
           social_scaled = scale(social_feature)[,1],
           roughness_scaled = scale(roughness_feature)[,1],
           loss_scaled = scale(loss_feature)[,1])
  
  # add success information
  subset$day = as.numeric(subset$day)
  subset$year = as.numeric(subset$year)
  subset$participant_id = as.numeric(subset$participant_id)
  subset$track_id = as.numeric(subset$track_id)
  subset$unique_spot_id = as.numeric(subset$unique_spot_id)
  subset$lake_id = as.numeric(subset$lake_id)
  subset = left_join(subset, spot_selection_data)
  
  # visualize
  #subset %>% 
  #  ggplot(aes(x = social_feature, y = Return)) + 
  #  stat_summary_bin(bins = 10) + 
  #  geom_smooth(method = "lm") + 
  #  facet_wrap(~lake, scales = "free")
  
  #subset %>% 
  #  ggplot(aes(x = social_feature, y = total_catch)) + 
  #  stat_summary_bin(bins = 10) + 
  #  geom_smooth(method = "lm") + 
  #  facet_wrap(~lake, scales = "free")
  
  if (!file.exists("utils/data/processed_data/cue_validity_model.rds")){
    model = brms::brm(data = subset, 
                      Return ~ success_scaled + social_scaled + roughness_scaled + loss_scaled + 
                        (success_scaled + social_scaled + roughness_scaled + loss_scaled|lake) + 
                        (success_scaled + social_scaled + roughness_scaled + loss_scaled|trip) + 
                        (success_scaled + social_scaled + roughness_scaled + loss_scaled|id),
                      family = "bernoulli",
                      chains = 4, 
                      cores = 4)
    #model_unscaled = brms::brm(data = subset, 
    #                  Return ~ success_feature + social_feature + roughness_feature + loss_feature + 
    #                    (success_feature + social_feature + roughness_feature + loss_feature|lake) + 
    #                    (success_feature + social_feature + roughness_feature + loss_feature|trip) + 
    #                    (success_feature + social_feature + roughness_feature + loss_feature|id),
    #                  family = "bernoulli",
    #                  chains = 4, 
    #                  cores = 4)
    saveRDS(model, "utils/data/processed_data/cue_validity_model.rds")
  } else {
    model = readRDS("utils/data/processed_data/cue_validity_model.rds")
  }
  
  # extract estimates for table
  summary(model)
  
  # visualize
  roughness = model %>%
    spread_draws(b_success_scaled,
                 b_social_scaled,
                 b_roughness_scaled,
                 b_loss_scaled,
                 r_lake[lake, parameter]) %>% 
    filter(parameter %in% c("roughness_scaled")) %>% 
    mutate(or = r_lake + b_roughness_scaled) #name or, but on logit scale
  social = model %>%
    spread_draws(b_success_scaled,
                 b_social_scaled,
                 b_roughness_scaled,
                 b_loss_scaled,
                 r_lake[lake, parameter]) %>% 
    filter(parameter %in% c("social_scaled")) %>% 
    mutate(or = r_lake + b_social_scaled)
  success = model %>%
    spread_draws(b_success_scaled,
                 b_social_scaled,
                 b_roughness_scaled,
                 b_loss_scaled,
                 r_lake[lake, parameter]) %>% 
    filter(parameter %in% c("success_scaled")) %>% 
    mutate(or = r_lake + b_success_scaled)
  loss = model %>%
    spread_draws(b_success_scaled,
                 b_social_scaled,
                 b_roughness_scaled,
                 b_loss_scaled,
                 r_lake[lake, parameter]) %>% 
    filter(parameter %in% c("loss_scaled")) %>% 
    mutate(or = r_lake + b_loss_scaled)
  
  parameters_all = rbind(roughness, social)
  parameters_all = rbind(parameters_all, loss)
  parameters_all = rbind(parameters_all, success)
  
  parameters_all$parameter = recode(parameters_all$parameter, "social_scaled" = "Social", "roughness_scaled" = "Roughness", "success_scaled" = "Successful", "loss_scaled" = "Unsuccessful")
  
  main_effect = model %>% 
    spread_draws(b_success_scaled,
                 b_social_scaled,
                 b_roughness_scaled,
                 b_loss_scaled)
  
  library(RColorBrewer)
  cols = brewer.pal(5, "Greens")
  #"#EDF8E9" "#BAE4B3" "#74C476" "#31A354" "#006D2C"
  adaptive_cues =
    parameters_all %>% 
    mutate(parameter = factor(parameter, levels=c("Roughness", "Social", "Unsuccessful", "Successful"))) %>%
    ggplot(aes(x = or, y =parameter, group = as.factor(lake))) + 
    geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
    # .95
    annotate("rect", xmin = quantile(main_effect$b_success_scaled, .025), xmax =  quantile(main_effect$b_success_scaled, .975), ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "#BAE4B3")+
    #.66
    annotate("rect", xmin = quantile(main_effect$b_success_scaled, .17), xmax =  quantile(main_effect$b_success_scaled, .83), ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "#31A354")+
    # .5
    annotate("rect", xmin = quantile(main_effect$b_success_scaled, .5)-.0025, xmax =  quantile(main_effect$b_success_scaled, .5)+.0025, ymin = 3.65, ymax = 4.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(main_effect$b_social_scaled, .025), xmax =  quantile(main_effect$b_social_scaled, .975), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill ="#BAE4B3")+
    annotate("rect", xmin = quantile(main_effect$b_social_scaled, .17), xmax =  quantile(main_effect$b_social_scaled, .83), ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "#31A354")+
    annotate("rect", xmin = quantile(main_effect$b_social_scaled, .5)-.0025, xmax =  quantile(main_effect$b_social_scaled, .5)+.0025, ymin = 1.65, ymax = 2.35,
             alpha = 1,fill = "black")+
    
    annotate("rect", xmin = quantile(main_effect$b_roughness_scaled, .025), xmax =  quantile(main_effect$b_roughness_scaled, .975), ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "#BAE4B3") +
    annotate("rect", xmin = quantile(main_effect$b_roughness_scaled, .17), xmax =  quantile(main_effect$b_roughness_scaled, .83),ymin = 0.65, ymax = 1.35,
             alpha = 1,fill ="#31A354") +
    annotate("rect", xmin = quantile(main_effect$b_roughness_scaled, .5)-.0025, xmax =  quantile(main_effect$b_roughness_scaled, .5)+.0025, ymin = 0.65, ymax = 1.35,
             alpha = 1,fill = "black") +
    
    annotate("rect", xmin = quantile(main_effect$b_loss_scaled, .025), xmax =  quantile(main_effect$b_loss_scaled, .975), ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "#BAE4B3") + 
    annotate("rect", xmin = quantile(main_effect$b_loss_scaled, .17), xmax =  quantile(main_effect$b_loss_scaled, .83),  ymin = 2.65, ymax = 3.35,
             alpha = 1,fill ="#31A354") + 
    annotate("rect", xmin = quantile(main_effect$b_loss_scaled, .5)-.0025, xmax =  quantile(main_effect$b_loss_scaled, .5)+.0025,  ymin = 2.65, ymax = 3.35,
             alpha = 1,fill = "black") + 
    stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0) +
    
    ylab("Feature")+
    labs(tag = "F", color = "value")+
    ggtitle("Feature Validity")+
    xlab("Posterior Est. (95%CrI)\nP(Successful) ~ Feature")+
    scale_x_continuous(breaks = seq(from = -2, to = 2, by = .5))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 6),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_blank(),#, hjust = .5
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  adaptive_cues
  
  ##############################################################################
  # REVISION: Visualize alternative roughness features 
  ##############################################################################
    # load fast dist function
    sourceCpp("utils/library/fast_dist.cpp")
    
    # simulate evaluation points across whole subset
    day5 = spot_selection_data[spot_selection_data$day == 5 & spot_selection_data$year == 2022,]
    box = c(min(min(day5$x), min(day5$x)), max(max(day5$x), max(day5$x)),
            min(min(day5$y), min(day5$y)), max(max(day5$y), max(day5$y)))
    
    box[2]-box[1]
    box[4]-box[3]
    
    bbox = rbind(c(box[1], box[3]),
                 c(box[2], box[3]),
                 c(box[1], box[4]),
                 c(box[2], box[4]))
    plot(bbox)
    sim_x = seq(from = box[1], to = box[2], by = 10)
    sim_y = seq(from = box[3], to = box[4], by = 10)
    sim_coords = expand_grid(sim_x, sim_y)
    sim_coords$x = sim_coords$sim_x
    sim_coords$y = sim_coords$sim_y
    
    # compute roughness feature
    lake_data = readRDS("utils/data/raw_data/lake_data/depth_profiles.rds")
    # get in same order
    depth_profiles = list()
    depth_profiles[[1]] <- lake_data[[7]]
    depth_profiles[[2]] <- lake_data[[10]]
    depth_profiles[[3]] <- lake_data[[6]]
    depth_profiles[[4]] <- lake_data[[9]]
    depth_profiles[[5]] <- lake_data[[8]]
    depth_profiles[[6]] <- lake_data[[1]]
    depth_profiles[[7]] <- lake_data[[2]]
    depth_profiles[[8]] <- lake_data[[3]]
    depth_profiles[[9]] <- lake_data[[5]]
    depth_profiles[[10]] <- lake_data[[4]]
    
    profile = depth_profiles[[9]]
    profile = profile %>% filter(depth > 0)
    
    # create depth model 
    pts = list()
    for (i in 1:nrow(profile)){
      pts[[i]] = st_line_sample(profile[i,], density = 1/10)
      #pts[[i]] = profile[i,] %>% st_cast("MULTIPOINT")
      pts[[i]] = pts[[i]] %>% st_cast("POINT")
      pts[[i]] = st_coordinates(pts[[i]])
      pts[[i]] = as.data.frame(pts[[i]])
      pts[[i]]$depth = pull(profile[i, "depth"] %>% st_drop_geometry())
    }
    pts = do.call(rbind, pts)
    colnames(pts) <- c("x", "y", "depth")
    
    # depth model linear interpolation: 
    gs <- gstat(formula = depth ~ 1, locations= ~x+y, data = pts, nmax = 50, set = list(idp = 0))
    
    sim_coords <- predict(gs, sim_coords)
    
    ##############################################################################
    if (!file.exists("utils/data/processed_data/roughness_features_plot.rds")){
    # simulate N = 100 points from disk of radius r=100
    N_pts = 100
    # simulate points using polar coordinate transformation
    roughness_matrix = matrix(nrow = nrow(sim_coords), ncol = 4)
    start = Sys.time()
    for (pt in 1:nrow(sim_coords)){
      
      # simulate random angle and radius within disk
      angles = runif(N_pts, -pi, pi)
      r = runif(N_pts, 0, 50)
      
      # simulate N points from location
      eval_points_x = sim_coords$x[pt] + r * cos(angles)
      eval_points_y = sim_coords$y[pt] + r * sin(angles)
      eval_df = data.frame(x = eval_points_x, y = eval_points_y)  
      
      # depth
      pred = predict(gs, eval_df)
      
      # turn into proper sf dataframe
      sf_df = st_as_sf(data.frame(x = eval_points_x, y = eval_points_y), coords = c("x", "y"))
      st_crs(sf_df) <- crs(profile)
      
      # distances to isoclines
      dists = st_distance(sf_df, profile)
      dist_to_isocline = apply(dists, 1, min)
      
      # compute variance in depth
      roughness_matrix[pt,1] = sd(pred$var1.pred)
      roughness_matrix[pt,2] = max(pred$var1.pred)-min(pred$var1.pred)
      roughness_matrix[pt,3] = mean(dist_to_isocline)
      roughness_matrix[pt,4] = mean(pred$var1.pred)
      
    }
    end = Sys.time()
    end-start
    saveRDS(roughness_matrix, file = "utils/data/processed_data/roughness_features_plot.rds")
  } else {
    
  }
  
  roughness_feature_matrix = readRDS("utils/data/processed_data/roughness_features_plot.rds")
  sim_coords$depth_variance = roughness_feature_matrix[,1]
  sim_coords$max_min_depth = roughness_feature_matrix[,2]
  sim_coords$mean_dist = roughness_feature_matrix[,3]
  sim_coords$mean_depth = roughness_feature_matrix[,4]
  
  roughness_feature_depth_variance = ggplot() + 
    geom_tile(data = sim_coords, aes(x = x, y = y, fill = depth_variance)) +
    geom_sf(data = profile, color ="lightgrey", alpha = .2, size = .6) +
    coord_sf(xlim = c(box[1], box[2]), ylim = c(box[3],box[4]))+
    #coord_fixed(ratio = 1)+
    scale_fill_viridis_c() + 
    ylab("Northing")+
    ggtitle("SD(Depth)")+
    labs(tag = "B", fill = expression(italic(f)))+
    xlab("Easting")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  roughness_feature_depth_variance
  
  roughness_feature_min_max = ggplot() + 
    geom_tile(data = sim_coords, aes(x = x, y = y, fill = max_min_depth)) +
    geom_sf(data = profile, color ="lightgrey", alpha = .2, size = .6) +
    coord_sf(xlim = c(box[1], box[2]), ylim = c(box[3],box[4]))+
    #coord_fixed(ratio = 1)+
    scale_fill_viridis_c() + 
    ylab("Northing")+
    ggtitle("Max(Depth)-Min(Depth)")+
    labs(tag = "C", fill = expression(italic(f)))+
    xlab("Easting")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  roughness_feature_min_max
  
  roughness_feature_depth = ggplot() + 
    geom_tile(data = sim_coords, aes(x = x, y = y, fill = mean_depth)) +
    geom_sf(data = profile, color ="lightgrey", alpha = .2, size = .6) +
    #geom_point(data = points, aes(x = X, y = Y), color = "darkgrey", alpha = .5)+
    coord_sf(xlim = c(box[1], box[2]), ylim = c(box[3],box[4]))+
    #coord_fixed(ratio = 1)+
    scale_fill_viridis_c() + 
    ylab("Northing")+
    ggtitle("Avg. Depth")+
    labs(tag = "A", fill = expression(italic(f)))+
    xlab("Easting")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  roughness_feature_depth
  
  # panel
  panel_ecological_features_suppl = gridExtra::arrangeGrob(
    grobs = list(roughness_feature_depth,
                 roughness_feature_depth_variance,
                 roughness_feature_min_max),
    layout_matrix = rbind(c(1, 2, 3))
  )
  
  ggsave(panel_ecological_features_suppl, file = "2_3_socio_ecological_features/output/panel_ecological_features_suppl.png", 
         width = 13.5, height = 5, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_ecological_features_suppl, file = "2_3_socio_ecological_features/output/panel_ecological_features_suppl.svg", 
         width = 13.5, height = 5, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_ecological_features_suppl, file = "2_3_socio_ecological_features/output/panel_ecological_features_suppl.pdf", 
         width = 13.5, height = 5, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_ecological_features_suppl, file = "2_3_socio_ecological_features/output/panel_ecological_features_suppl.tiff", 
         width = 13.5, height = 5, bg = "white", dpi = 600, units = "cm")
  
  ##############################################################################
  # visualize model comparison alternative ecological features
  comparison = readRDS(file = "utils/data/processed_data/model_fit/spot_selection_models/comparison/comparison_ecological_features.rds")
  comparison = as.data.frame(comparison)
  comparison$model_id = c(4, 1, 2, 3, 5)
  comparison_ecological_features = comparison %>%
    ggplot(aes(x = elpd_diff, y = fct_reorder(name, -elpd_diff))) + 
    geom_col(color = "black", width = .6) +
    geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
    #ggtitle("Model Comparison Spot-Selection") + 
    xlab(expression(Delta~"ELPD (rel. to best)")) + 
    ylab("Model") + 
    #labs(tag = "a")+
    scale_x_continuous(labels = scales::comma_format(big.mark = ','))+
    geom_vline(xintercept = 0)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7, angle = 65, hjust = .5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_blank(),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  comparison_ecological_features
  
  ggsave(comparison_ecological_features, file = "2_3_socio_ecological_features/output/comparison_ecological_features.png", 
         width = 8, height = 8, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_ecological_features, file = "2_3_socio_ecological_features/output/comparison_ecological_features.svg", 
         width = 8, height = 8, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_ecological_features, file = "2_3_socio_ecological_features/output/comparison_ecological_features.pdf", 
         width = 8, height = 8, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_ecological_features, file = "2_3_socio_ecological_features/output/comparison_ecological_features.tiff", 
         width = 8, height = 8, bg = "white", dpi = 600, units = "cm")
  
  ################################################################################
  # Visualize model fit procedure + spatial features
  ################################################################################
  # take example trajectory
  traj = readRDS("utils/data/raw_data/data_figure2/sub.rds")
  
  # load fast dist function
  sourceCpp("utils/library/fast_dist.cpp")
  
  # simulate evaluation points across whole lake
  sim_spots = readRDS(file = "utils/data/raw_data/sim_spots.rds")
  
  # simulate evaluation points across whole subset
  box = c(min(min(traj$x), min(sim_spots$sim_spots_x)) - 37.5, max(max(traj$x), max(sim_spots$sim_spots_x)) + 37.5,
          min(min(traj$y), min(sim_spots$sim_spots_y)) - 4, max(max(traj$y), max(sim_spots$sim_spots_y)) + 4)
  box = c(min(min(traj$x), min(sim_spots$sim_spots_x)) - 37.5, max(max(traj$x), max(sim_spots$sim_spots_x)) + 37.5,
          min(min(traj$y), min(sim_spots$sim_spots_y)) - 4, max(max(traj$y), max(sim_spots$sim_spots_y)) + 4)
  
  box[2]-box[1]
  box[4]-box[3]
  
  bbox = rbind(c(box[1], box[3]),
               c(box[2], box[3]),
               c(box[1], box[4]),
               c(box[2], box[4]))
  plot(bbox)
  points(sim_spots)
  sim_x = seq(from = box[1], to = box[2], by = .5)
  sim_y = seq(from = box[3], to = box[4], by = .5)
  sim_coords = expand_grid(sim_x, sim_y)
  
  # compute success feature
  last_successful_spot = traj[c(3, 4), c("x", "y")]
  dist = fastPdist2(as.matrix(sim_coords), as.matrix(last_successful_spot))
  success_feature = kern(d = dist,h = exp(2.22), option = "success", ncol = 2)
  sim_coords$success_feature = success_feature
  
  # compute loss feature
  last_loss = traj[c(1, 2, 5), c("x", "y")]
  dist = fastPdist2(as.matrix(sim_coords), as.matrix(last_loss))
  loss_feature = kern(d = dist,h = exp(2.22), option = "success", ncol = 3)
  
  sim_coords$loss_feature = loss_feature
  
  # compute social feature
  coords = readRDS("utils/data/raw_data/data_figure2/coords.rds")
  dist = fastPdist2(as.matrix(sim_coords), as.matrix(coords[, 2:3]))
  social_feature = kern(d = dist,h = exp(3.19), option = "success", ncol = 50)
  sim_coords$social_feature = social_feature
  
  # compute roughness feature
  lake_data = readRDS("utils/data/raw_data/lake_data/depth_profiles.rds")
  # get in same order
  depth_profiles = list()
  depth_profiles[[1]] <- lake_data[[7]]
  depth_profiles[[2]] <- lake_data[[10]]
  depth_profiles[[3]] <- lake_data[[6]]
  depth_profiles[[4]] <- lake_data[[9]]
  depth_profiles[[5]] <- lake_data[[8]]
  depth_profiles[[6]] <- lake_data[[1]]
  depth_profiles[[7]] <- lake_data[[2]]
  depth_profiles[[8]] <- lake_data[[3]]
  depth_profiles[[9]] <- lake_data[[5]]
  depth_profiles[[10]] <- lake_data[[4]]
  
  profile = depth_profiles[[9]]
  profile = profile %>% filter(depth > 0)
  
  depth_profile_points = list()
  for (i in 1:10){
    tmp = st_line_sample(depth_profiles[[i]] %>% filter(depth>0), density = .1, type = "regular")
    depth_profile_points[[i]] = tmp[!st_is_empty(tmp)] %>% st_cast("POINT")
  }
  
  point_profile = depth_profile_points[[9]]
  
  # compute distances
  points = st_coordinates(point_profile)
  d = fastPdist2(as.matrix(sim_coords)[, c("sim_x", "sim_y")], points)
  
  # if very small sometimes na
  d[which(is.na(d))] <- .00001
  
  # keep only 50 closest features for storage
  d = t(apply(d, 1, sort))[, 1:50]
  
  roughness_feature = as.matrix(d[, 1])
  sim_coords$roughness_feature = roughness_feature/1000
  ##############################################################################
  # REVISION:
  # softmax visualization
  sim_coords$social_feature = ifelse(is.na(sim_coords$social_feature), 0 , sim_coords$social_feature)
  sum(is.na(sim_coords$loss_feature))
  sim_coords$loss_feature = ifelse(is.na(sim_coords$loss_feature), 0.001, sim_coords$loss_feature)
  spot_selection_model_fit
  s = softmax(2.35 * sim_coords$success_feature -1.03 * sim_coords$loss_feature + .89 * sim_coords$social_feature - 4.48 * sim_coords$roughness_feature)
  sim_coords$softmax = s
  
  softmax_feature = ggplot() + 
    geom_tile(data = sim_coords, aes(x = sim_x, y = sim_y, fill = -1/softmax), alpha = 1) +
    #geom_point(data = sim_coords, aes(x = sim_x, y = sim_y, color = success_feature), alpha = .2) +
    geom_path(data = traj[1:5,], aes(x = x, y = y), color = "chartreuse", size = .3)+
    geom_path(data = traj[5:6,], aes(x = x, y = y), color = "chartreuse", linetype = "dotted", size = .3)+
    geom_point(data = sim_spots, aes(x = sim_spots_x, y = sim_spots_y), size = .6, stroke = .6,color = "darkgrey", alpha = .75)+
    #geom_point(data = traj[c(1, 2, 5),], aes(x = x, y = y), size = .6, stroke = .6, color = "chartreuse")+
    geom_point(data = traj[c(1, 2,3, 4, 5),], aes(x = x, y = y), size = .6, stroke = .6, color = "chartreuse")+
    geom_point(data = traj[1:5,], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 21)+
    #geom_point(data = traj[c(3, 4),], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 8)+ #, shape = 8
    #geom_point(data = traj[c(1, 2, 5),], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 21)+ #, shape = 21
    geom_point(data = traj[6,], aes(x = x, y = y), stroke = .6,color = "red", shape = 4, size = .6)+
    geom_point(data = coords %>% filter(between(X, box[1], box[2]) & between(Y, box[3], box[4])), aes(x = X, y = Y), size = .6, shape = 17, color = "black", alpha = .8)+
    geom_sf(data = profile, size = 1, color = "lightgrey") + 
    scale_fill_viridis_c()+
    #scale_color_distiller()+
    ylab("Northing")+
    ggtitle("Model Prediction")+
    #labs(tag = "A")+
    xlab("Easting")+
    #theme_classic()+
    coord_sf(xlim = c(box[1], box[2]), ylim = c(box[3],box[4]))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white")) + 
    annotation_scale(text_col = "white")
  softmax_feature
  ggsave(softmax_feature, file = "2_3_socio_ecological_features/output/softmax_feature.png", 
         width = 6, height = 6, bg = "white", dpi = 600, units = "cm")
  ggsave(softmax_feature, file = "2_3_socio_ecological_features/output/softmax_feature.svg", 
         width = 6, height = 6, bg = "white", dpi = 600, units = "cm",device = svglite::svglite)
  ggsave(softmax_feature, file = "2_3_socio_ecological_features/output/softmax_feature.pdf", 
         width = 6, height = 6, bg = "white", dpi = 600, units = "cm")
  ggsave(softmax_feature, file = "2_3_socio_ecological_features/output/softmax_feature.tiff", 
         width = 6, height = 6, bg = "white", dpi = 600, units = "cm")
  
  # visualize
  gc()
  gc()
  success_feature = ggplot() + 
    geom_tile(data = sim_coords, aes(x = sim_x, y = sim_y, fill = success_feature), alpha = 1) +
    #geom_point(data = sim_coords, aes(x = sim_x, y = sim_y, color = success_feature), alpha = .2) +
    geom_path(data = traj[1:5,], aes(x = x, y = y), color = "chartreuse", size = .3)+
    geom_path(data = traj[5:6,], aes(x = x, y = y), color = "chartreuse", linetype = "dotted", size = .3)+
    geom_point(data = sim_spots, aes(x = sim_spots_x, y = sim_spots_y), size = .6, stroke = .6,color = "darkgrey", alpha = .75)+
    #geom_point(data = traj[c(1, 2, 5),], aes(x = x, y = y), size = .6, stroke = .6, color = "chartreuse")+
    geom_point(data = traj[c(1, 2,3, 4, 5),], aes(x = x, y = y), size = .6, stroke = .6, color = "chartreuse")+
    geom_point(data = traj[1:5,], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 21)+
    #geom_point(data = traj[c(3, 4),], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 8)+ #, shape = 8
    #geom_point(data = traj[c(1, 2, 5),], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 21)+ #, shape = 21
    geom_point(data = traj[6,], aes(x = x, y = y), stroke = .6,color = "red", shape = 4, size = .6)+
    scale_fill_viridis_c()+
    #scale_color_distiller()+
    ylab("Northing")+
    ggtitle("Successful Spots")+
    labs(tag = "A")+
    xlab("Easting")+
    #theme_classic()+
    coord_sf(xlim = c(box[1], box[2]), ylim = c(box[3],box[4]))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white")) + 
    annotation_scale(text_col = "white")
  success_feature
  
  loss_feature = ggplot() + 
    geom_tile(data = sim_coords, aes(x = sim_x, y = sim_y, fill = loss_feature), alpha = 1) +
    #geom_point(data = sim_coords, aes(x = sim_x, y = sim_y, color = success_feature), alpha = .2) +
    geom_path(data = traj[1:5,], aes(x = x, y = y), color = "chartreuse", size = .3)+
    geom_path(data = traj[5:6,], aes(x = x, y = y), color = "chartreuse", linetype = "dotted", size = .3)+
    geom_point(data = sim_spots, aes(x = sim_spots_x, y = sim_spots_y), size = .6, stroke = .6,color = "darkgrey", alpha = .75)+
    geom_point(data = traj[c(1, 2, 3, 4, 5),], aes(x = x, y = y), size = .6, stroke = .6, color = "chartreuse")+
    geom_point(data = traj[1:5,], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 21)+
    #geom_point(data = traj[c(3, 4),], aes(x = x, y = y), size = .6, stroke = .6, color = "chartreuse")+
    #geom_point(data = traj[c(1, 2, 5),], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 13)+ #, shape = 8
    #geom_point(data = traj[c(3, 4),], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 21)+ #, shape = 21
    geom_point(data = traj[6,], aes(x = x, y = y), stroke = .6,color = "red", shape = 4, size = .6)+
    scale_fill_viridis_c()+
    #scale_color_distiller()+
    ylab("Northing")+
    ggtitle("Unsuccessful Spots")+
    labs(tag = "B", fill = expression(italic(f)))+
    xlab("Easting")+
    #theme_classic()+
    coord_sf(xlim = c(box[1], box[2]), ylim = c(box[3],box[4]))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  loss_feature
  
  social_feature = ggplot() + 
    geom_tile(data = sim_coords, aes(x = sim_x, y = sim_y, fill = social_feature)) +
    geom_path(data = traj[1:5,], aes(x = x, y = y), color = "chartreuse", size = .3)+
    geom_path(data = traj[5:6,], aes(x = x, y = y), color = "chartreuse", linetype = "dotted", size = .3)+
    geom_point(data = sim_spots, aes(x = sim_spots_x, y = sim_spots_y), size = .6, stroke = .6,color = "darkgrey", alpha = .75)+
    geom_point(data = traj[1:5,], aes(x = x, y = y), size = .6, stroke = .6, color = "chartreuse")+
    geom_point(data = traj[1:5,], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 21)+
    geom_point(data = traj[6,], aes(x = x, y = y), size = .6, stroke = .6,color = "red", shape = 4)+
    geom_point(data = coords %>% filter(between(X, box[1], box[2]) & between(Y, box[3], box[4])), aes(x = X, y = Y), size = .6, shape = 17, color = "black", alpha = .8)+
    scale_fill_viridis_c() + 
    ylab("Northing")+
    coord_sf(xlim = c(box[1], box[2]), ylim = c(box[3],box[4]))+
    ggtitle("Social Density")+
    labs(tag = "C", fill = expression(italic(f)))+
    xlab("Easting")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  social_feature
  
  roughness_feature = ggplot() + 
    geom_tile(data = sim_coords, aes(x = sim_x, y = sim_y, fill = roughness_feature)) +
    geom_path(data = traj[1:5,], aes(x = x, y = y), color = "chartreuse", size = .3)+
    geom_path(data = traj[5:6,], aes(x = x, y = y), color = "chartreuse", linetype = "dotted", size = .3)+
    geom_point(data = sim_spots, aes(x = sim_spots_x, y = sim_spots_y), size = .6, stroke = .6,color = "darkgrey", alpha = .75)+
    geom_point(data = traj[1:5,], aes(x = x, y = y), size = .6, stroke =.6, color = "chartreuse")+
    geom_point(data = traj[1:5,], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 21)+
    #geom_point(data = traj[1:5,], aes(x = x, y = y), size = .6, stroke = .6, color = "black", shape = 21)+
    geom_point(data = traj[6,], aes(x = x, y = y), size = .6, stroke = .6,color = "red", shape = 4)+
    geom_sf(data = profile, size = 1, color = "lightgrey") + 
    #geom_point(data = points, aes(x = X, y = Y), color = "darkgrey", alpha = .5)+
    coord_sf(xlim = c(box[1], box[2]), ylim = c(box[3],box[4]))+
    #coord_fixed(ratio = 1)+
    scale_fill_viridis_c() + 
    ylab("Northing")+
    ggtitle("Lake-bed Roughness")+
    labs(tag = "D", fill = expression(italic(f)))+
    xlab("Easting")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          axis.text.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  roughness_feature
  
  
  
  # panel
  panel_model_results = gridExtra::arrangeGrob(
    grobs = list(success_feature,
                 loss_feature,
                 social_feature,
                 roughness_feature,
                 model_results,
                 model_results_interaction,
                 bandwidth_results,
                 adaptive_cues),
    widths = rep(1, 20),
    heigths = c(1, .8),
    layout_matrix = rbind(rep(1:4, each = 5),
                          rep(c(5, 8, 6, 7), each = 5))
  )
  
  ggsave(panel_model_results, file = "2_3_socio_ecological_features/output/panel_model_results.png", 
         width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_model_results, file = "2_3_socio_ecological_features/output/panel_model_results.svg", 
         width = 18, height = 10, bg = "white", dpi = 600, units = "cm",device = svglite::svglite)
  ggsave(panel_model_results, file = "2_3_socio_ecological_features/output/panel_model_results.pdf", 
         width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_model_results, file = "2_3_socio_ecological_features/output/panel_model_results.tiff", 
         width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
  
  ################################################################################
  # Supplementary plot: visualize bandwidth generalization
  generalization_plot = spot_selection_model_fit %>% 
    spread_draws(lambdas[parameter]) %>% 
    slice(rep(1:n(), each = 500)) %>% 
    group_by(parameter, .chain, .iteration, .draw) %>% 
    mutate(d = seq(from = 0, to = 100, length.out = 500)) %>% 
    mutate(g = exp(-d^2 / (2 * exp(lambdas)^2))) %>% 
    mutate(parameter = ifelse(parameter == 1, "Social", "Personal")) %>% 
    ggplot(aes(x = d, y = g, color = parameter)) + 
    stat_lineribbon(fill = "lightgrey") +
    scale_color_manual(values = c("green4","orange2"))+
    labs(color = "Information", fill = "CrI")+
    xlab("Distance to Feature (m)")+
    ylab("Feature Value")+
    #ggtitle("Spatial Generalization") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          axis.text.y = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .5),
          legend.text = element_text(size = 7),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = .2))
  ggsave(generalization_plot, file = "2_3_socio_ecological_features/output/generalization_vis.png", 
         width = 9, height = 9, bg = "white", dpi = 600, units = "cm")
  ggsave(generalization_plot, file = "2_3_socio_ecological_features/output/generalization_vis.svg", 
         width = 9, height = 9, bg = "white", dpi = 600, units = "cm")
  ggsave(generalization_plot, file = "2_3_socio_ecological_features/output/generalization_vis.pdf", 
         width = 9, height = 9, bg = "white", dpi = 600, units = "cm")
  ggsave(generalization_plot, file = "2_3_socio_ecological_features/output/generalization_vis.tiff", 
         width = 9, height = 9, bg = "white", dpi = 600, units = "cm")
  
  ################################################################################
  # EDF: Model comparison
  comparison = readRDS(file = "utils/data/processed_data/model_fit/spot_selection_models/comparison/comparison.rds")
  comparison = as.data.frame(comparison)
  comparison$name = c("+Conditional\nSocial Info",
                      "+Social",
                      "+Personal",
                      "+Ecological",
                      "+Locality",
                      "Random")
  comparison$model_id = 1:6
  comparison_spot_selection = comparison %>%
    ggplot(aes(x = elpd_diff, y = fct_reorder(name, model_id))) + 
    geom_col(color = "black", width = .6) +
    geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
    #ggtitle("Model Comparison Spot-Selection") + 
    xlab(expression(Delta~"ELPD (rel. to best)")) + 
    ylab("Model") + 
    #labs(tag = "a")+
    scale_x_continuous(labels = scales::comma_format(big.mark = ','))+
    geom_vline(xintercept = 0)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7, angle = 65, hjust = .5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_blank(),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  comparison_spot_selection
  
  ggsave(comparison_spot_selection, file = "2_3_socio_ecological_features/output/comparison_spot_selection.png", 
         width = 8, height = 8, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_spot_selection, file = "2_3_socio_ecological_features/output/comparison_spot_selection.svg", 
         width = 8, height = 8, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_spot_selection, file = "2_3_socio_ecological_features/output/comparison_spot_selection.pdf", 
         width = 8, height = 8, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_spot_selection, file = "2_3_socio_ecological_features/output/comparison_spot_selection.tiff", 
         width = 8, height = 8, bg = "white", dpi = 600, units = "cm")
  
  ################################################################################
  # REVISION: Model comparison across lakes
  comparison_lake = readRDS(file = "utils/data/processed_data/model_fit/spot_selection_models/comparison/comparison_spot_selection_lakes.rds")
  comparison_spot_selection_lake = comparison_lake %>%
    ggplot(aes(y = fct_reorder(model_name, -elpd_diff), x = elpd_diff)) +
    geom_col(color = "black", width = .6) +
    geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
    xlab(expression(Delta~"ELPD (rel. to best)")) + 
    ylab("Model") + 
    facet_wrap(~lake_id,nrow = 2)+
    scale_x_continuous(labels = scales::comma_format(big.mark = ','))+
    geom_vline(xintercept = 0)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 5, angle = 65, hjust = .5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_blank(),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  comparison_spot_selection_lake
  
  ggsave(comparison_spot_selection_lake, file = "2_3_socio_ecological_features/output/comparison_spot_selection_lake.png", 
         width = 18, height = 11, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_spot_selection, file = "2_3_socio_ecological_features/output/comparison_spot_selection_lake.svg", 
         width = 18, height = 11, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_spot_selection, file = "2_3_socio_ecological_features/output/comparison_spot_selection_lake.pdf", 
         width = 18, height = 11, bg = "white", dpi = 600, units = "cm")
  ggsave(comparison_spot_selection, file = "2_3_socio_ecological_features/output/comparison_spot_selection_lake.tiff", 
         width = 18, height = 11, bg = "white", dpi = 600, units = "cm")
  
  ##############################################################################
  # supplement: visualize distributions for model
  steps = spot_selection_data$step[spot_selection_data$step > -99 & spot_selection_data$angle > -99]
  steps = ifelse(steps == 0, .1, steps)
  angles = spot_selection_data$angle[spot_selection_data$step > -99 & spot_selection_data$angle > -99]
  sum(steps>300)
  # plot for supplement
  vis = data.frame(steps = steps,
                   angles = angles)
  
  a=ggplot(vis, aes(x = steps)) +
    geom_density(color = "black", fill = "darkgrey") + 
    coord_cartesian(xlim = c(0, 300)) + 
    ylab("Density")+
    labs(tag = "A", color = "value")+
    xlab("Relocation Distance (m)")+
    scale_x_continuous(breaks = seq(from = 0, to = 300, by = 100))+
    scale_y_continuous(breaks = seq(from = 0, to = .03, by = .01))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  b=ggplot(vis, aes(x = angles)) +
    geom_density(color = "black", fill = "darkgrey")+
    #geom_histogram(bins = 50, color = "black", fill = "darkgrey") + 
    #coord_cartesian(xlim = c(-3.5, 3.5)) + 
    ylab("Density")+
    labs(tag = "B", color = "value")+
    xlab("Turning Angle (rad.)")+
    scale_x_continuous(breaks = seq(from = -4, to = 4, by = 2))+
    scale_y_continuous(breaks = seq(from = 0, to = .4, by = .1))+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  
  # panel
  panel_distributions = arrangeGrob(
    grobs = list(a, b),
    widths = c(1, 1),
    layout_matrix = matrix(1:2, nrow = 1, ncol = 2)
  )
  
  ggsave(panel_distributions, file = "2_3_socio_ecological_features/output/panel_distributions.png", 
         width = 9, height = 5, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_distributions, file = "2_3_socio_ecological_features/output/panel_distributions.svg", 
         width = 9, height = 5, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_distributions, file = "2_3_socio_ecological_features/output/panel_distributions.pdf", 
         width = 9, height = 5, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_distributions, file = "2_3_socio_ecological_features/output/panel_distributions.tiff", 
         width = 9, height = 5, bg = "white", dpi = 600, units = "cm")
  
}
################################################################################
# END
################################################################################
