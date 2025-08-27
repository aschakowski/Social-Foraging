################################################################################
#
# Title: Visualize patch leaving results
#
# Description: 
#
# Authors: Schakowski, A.
#
# Last updated: 02/04/2024 (DD/MM/YYYY)
#
################################################################################
visualize_patch_leaving <- function(){

# load model estimates
patch_leaving_model_fit = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/e_patch_leaving_model_bandwidth_fit.rds")
stan_data = readRDS(file = "utils/data/processed_data/stan_data_patch_leaving.rds")

# summarize as table
summary = patch_leaving_model_fit$summary(variables = c("beta", "alpha","lambda", "sigma_lake", "sigma_trip", "sigma_id"),~mean(.x), ~quantile(.x, probs = c(0.025, 0.975)))
print(summary, n = 50)
# parameter estimates and convergence diagnostics
# parameter
# beta1 = leave baseline
# beta2 = patch discovery
# beta3 = time without fish
# beta4 = global reward
# beta5 = local reward
# beta6 = success feature
# beta7 = loss feature
# beta8 = roughness feature
# beta9 = social feature
# alpha1 = global updating rate
# alpha2 = local updating rate

#patch_leaving_model_fit = test_fit
# visualize parameter estimates
if(!file.exists("utils/data/processed_data/main_effects_patch_leaving.rds")){
  main_effects = patch_leaving_model_fit %>% 
    spread_draws(beta[parameter]) %>% 
    pivot_wider(names_from = parameter, values_from = beta)
  saveRDS(main_effects, file = "utils/data/processed_data/main_effects_patch_leaving.rds")
  lake_effects = patch_leaving_model_fit %>% 
    spread_draws(v_lake[lake, parameter]) %>% 
    filter(parameter < 10)
  lake_effects = left_join(lake_effects, main_effects)
  saveRDS(lake_effects, file = "utils/data/processed_data/lake_effects_patch_leaving.rds")
  
} else {
  main_effects = readRDS("utils/data/processed_data/main_effects_patch_leaving.rds")
  lake_effects = readRDS("utils/data/processed_data/lake_effects_patch_leaving.rds")
}
gc()
gc()

lake_effects = lake_effects %>% 
  mutate(lake_effect = ifelse(parameter == 1, v_lake + `1`, 
                              ifelse(parameter == 2, v_lake + `2`,
                                     ifelse(parameter == 3, v_lake + `3`,
                                            ifelse(parameter == 4, v_lake + `4`,
                                                   ifelse(parameter == 5, v_lake + `5`,
                                                          ifelse(parameter == 6, v_lake + `6`,
                                                                 ifelse(parameter == 7, v_lake + `7`,
                                                                        ifelse(parameter == 8, v_lake + `8`, v_lake + `9`)))))))))

# scale effect sizes
if (!file.exists("utils/data/processed_data/fitted_features_patch_leaving.rds")){
  source("4_incremental_strategies/scripts/simulate_fitted_features_patch_leaving.r")
  simulate_fitted_features_patch_leaving()
}
fitted_features = readRDS(file = "utils/data/processed_data/fitted_features_patch_leaving.rds")

scale_features = fitted_features %>% 
  unnest() %>% 
  group_by(unique_trip, lake) %>% 
  summarize(sd_social = sd(social_feature),
            sd_roughness = sd(roughness_feature),
            sd_success = sd(success_feature),
            sd_loss = sd(loss_feature)) %>% 
  ungroup() %>% 
  group_by(lake) %>% 
  summarize(sd_social = mean(sd_social),
            sd_roughness = mean(sd_roughness),
            sd_success = mean(sd_success),
            sd_loss = mean(sd_loss)) %>% 
  unnest()

lake_effects = left_join(lake_effects, scale_features)

# compute standardized predictors
lake_effects = lake_effects %>% 
  ungroup() %>% 
  mutate(standardized_weight = ifelse(parameter == 6, lake_effect*sd_success, 
                                      ifelse(parameter == 7, lake_effect*sd_loss, 
                                             ifelse(parameter == 8, lake_effect*sd_roughness, 
                                                    ifelse(parameter == 9, lake_effect*sd_social, 0)))),
         standardized_main_effect = ifelse(parameter == 6, `6` * mean(scale_features$sd_success), 
                                           ifelse(parameter == 7, `7` * mean(scale_features$sd_loss),
                                                  ifelse(parameter == 8, `8` * mean(scale_features$sd_roughness),
                                                         ifelse(parameter == 9, `9` * mean(scale_features$sd_social),
                                                                0)))))

lake_effects$parameter = lake_effects$parameter %>% recode("Baseline", "Fish\nDiscovery", "Time w/o\nCatch (min)", "Waiting Time\nGlobal (min)", "Waiting Time\nLocal (min)", "Successful", "Unsuccessful", "Roughness", "Social")
catch_bias = lake_effects %>% 
  mutate(name = factor(parameter, levels = rev(c("Time w/o\nCatch (min)", "Fish\nDiscovery")))) %>%
  filter(name %in% (c("Time w/o\nCatch (min)", "Fish\nDiscovery"))) %>% 
  mutate(lake_effect = ifelse(parameter == "Time w/o\nCatch (min)", lake_effect * 6, lake_effect)) %>% 
  ggplot(aes(x = lake_effect, y = name, group = as.factor(lake))) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
  # .95
  annotate("rect", xmin = (quantile(main_effects$`3`*6, .025)), xmax =  (quantile(main_effects$`3`*6, .975)), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = (quantile(main_effects$`3`*6, .17)), xmax =  (quantile(main_effects$`3`*6, .83)), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = (quantile(main_effects$`3`*6, .5))-.005, xmax =  (quantile(main_effects$`3`*6, .5))+.005, ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "black")+
  
  annotate("rect", xmin = (quantile(main_effects$`2`, .025)), xmax =  (quantile(main_effects$`2`, .975)), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#BAE4B3")+
  annotate("rect", xmin = (quantile(main_effects$`2`, .17)), xmax =  (quantile(main_effects$`2`, .83)), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#31A354")+
  annotate("rect", xmin = (quantile(main_effects$`2`, .5))-.005, xmax =  (quantile(main_effects$`2`, .5))+.005, ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "black")+
  
  stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0) +
  
  #
  #ylab("")+
  labs(tag = "A")+
  ggtitle("Success-dependent GUT")+
  xlab("Posterior Est. (95%CrI)\nlogit P(Leave) ~ Feature")+
  ylab("Predictor")+
  #scale_x_continuous(breaks = seq(from = 0, to = 0.08, by = .02))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
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
catch_bias

cue_weights = lake_effects %>% 
  filter(parameter %in% c("Waiting Time\nGlobal (min)", "Waiting Time\nLocal (min)")) %>% 
  mutate(name = factor(parameter, levels = c("Waiting Time\nGlobal (min)", "Waiting Time\nLocal (min)"))) %>%
  ggplot(aes(x = (lake_effect)*6, y = name, group = as.factor(lake))) + 
  #geom_vline(xintercept = 0, linetype = "dashed") + 
  #stat_pointinterval(position = position_dodge(width = .6), color = "black", alpha = .25, size = 5) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
  # .95
  annotate("rect", xmin = (quantile(main_effects$`4`*6, .025)), xmax =  (quantile(main_effects$`4`*6, .975)), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = (quantile(main_effects$`4`*6, .17)), xmax =  (quantile(main_effects$`4`*6, .83)), ymin = .65, ymax = 1.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = (quantile(main_effects$`4`*6, .5))-.0025, xmax =  (quantile(main_effects$`4`*6, .5))+.0025, ymin = .65, ymax = 1.35,
           alpha = 1,fill = "black")+
  
  annotate("rect", xmin = (quantile(main_effects$`5`*6, .025)), xmax =  (quantile(main_effects$`5`*6, .975)), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = (quantile(main_effects$`5`*6, .17)), xmax =  (quantile(main_effects$`5`*6, .83)), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = (quantile(main_effects$`5`*6, .5))-.0025, xmax =  (quantile(main_effects$`5`*6, .5))+.0025, ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "black")+
  stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0) +
  
  ylab("Predictor")+
  labs(tag = "C")+
  ggtitle("Updating")+
  xlab("Posterior Est. (95%CrI)\nlogit P(Leave) ~ Feature")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
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
cue_weights

spatial_weights = lake_effects %>% 
  filter(parameter %in% c("Successful","Unsuccessful", "Social", "Roughness")) %>% 
  mutate(name = factor(parameter, levels = c("Roughness","Social", "Unsuccessful", "Successful"))) %>%
  ggplot(aes(x = (standardized_weight), y = name, group = as.factor(lake))) + 
  #stat_pointinterval(position = position_dodge(width = .6), color = "black", alpha = .25, size = 5) + 
  geom_vline(xintercept = 0, linetype = "dashed", size = .4) + 
  # .95
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Successful"], .025)), xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Successful"], .975)), ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Successful"], .17)), xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Successful"], .83)), ymin = 3.65, ymax = 4.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Successful"], .5))-.001, xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Successful"], .50))+.001, ymin =3.65, ymax = 4.35,
           alpha = 1,fill = "black")+
  
  
  # .95
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Unsuccessful"], .025)), xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Unsuccessful"], .975)), ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Unsuccessful"], .17)), xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Unsuccessful"], .83)), ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Unsuccessful"], .5))-.001, xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Unsuccessful"], .5))+.001, ymin = 2.65, ymax = 3.35,
           alpha = 1,fill = "black")+
  
  # .95
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Roughness"], .025)), xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Roughness"], .975)), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Roughness"], .17)), xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Roughness"], .83)), ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Roughness"], .5))-.001, xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Roughness"], .50))+.001, ymin = 0.65, ymax = 1.35,
           alpha = 1,fill = "black")+
  
  # .95
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Social"], .025)), xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Social"], .975)), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#BAE4B3")+
  #.66
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Social"], .17)), xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Social"], .83)), ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "#31A354")+
  # .5
  annotate("rect", xmin = (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Social"], .5))-.001, xmax =  (quantile(lake_effects$standardized_main_effect[lake_effects$parameter == "Social"], .50))+.001, ymin = 1.65, ymax = 2.35,
           alpha = 1,fill = "black")+
  stat_pointinterval(position = position_dodge(width = .8), color = "black", alpha = .10, size = 0) +
  ylab("Feature")+
  labs(tag = "D")+
  ggtitle("Spatial Features")+
  xlab("Posterior Est. (95%CrI)\nlogit P(Leave) ~ Feature")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6, angle = 55, hjust = .5),#, hjust = .5
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.9, .2),
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
spatial_weights

gc()
################################################################################
# REVISION: visualize model comparisons
comparison_basic = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/comparison/comparison_basic_strategies.rds")
comparison_basic = as.data.frame(comparison_basic)
comparison_basic$name = c("GUT","Fixed-Time","Fixed-n","Incremental", "Baseline")
comparison_basic$model_id = c(2, 5, 4, 3, 1)
model_comparison = comparison_basic %>% 
  filter(model_id <= 5) %>% 
  ggplot(aes(x = elpd_diff, y = fct_reorder(name, -elpd_diff))) + 
  geom_col(color = "black", width = .6) +
  geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
  #ggtitle("Model Comparison Patch-Leaving") + 
  xlab(expression(Delta~"ELPD (rel. to best)")) + 
  ylab("Model") + 
  labs(tag = "A")+
  ggtitle("Basic Strategies")+
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
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
model_comparison

#
comparison_extended = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/comparison/comparison_extended_strategies.rds")
comparison_extended = as.data.frame(comparison_extended)
comparison_extended$name = c("+Spatial\nFeatures", "+Updating", "+Fish\nDiscovery", "GUT")
comparison_extended$model_id = 1:4
model_comparison_extended = comparison_extended %>% 
  filter(model_id <= 4) %>% 
  ggplot(aes(x = elpd_diff, y = fct_reorder(name, -elpd_diff))) + 
  geom_col(color = "black", width = .6) +
  geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
  #ggtitle("Model Comparison Patch-Leaving") + 
  xlab(expression(Delta~"ELPD (rel. to best)")) + 
  ylab("Model") + 
  labs(tag = "B")+
  ggtitle("Extended GUT Rule")+
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
        title = element_text(size = 7),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "white", fill = "white"),
        text = element_text(family = "sans"))
model_comparison_extended

# make panel
panel_tmp_comparison = arrangeGrob(
  grobs = list(model_comparison,
               model_comparison_extended),
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 2))
)



# save figure
ggsave(panel_tmp_comparison, file = "4_incremental_strategies/output/model_comparison_patch_leaving.png", 
       width = 16, height = 8, bg = "white", dpi = 600, units = "cm")
ggsave(panel_tmp_comparison, file = "4_incremental_strategies/output/model_comparison_patch_leaving.svg", 
       width = 16, height = 8, bg = "white", dpi = 600, units = "cm")
ggsave(panel_tmp_comparison, file = "4_incremental_strategies/output/model_comparison_patch_leaving.pdf", 
       width = 16, height = 8, bg = "white", dpi = 600, units = "cm")
ggsave(panel_tmp_comparison, file = "4_incremental_strategies/output/model_comparison_patch_leaving.tiff", 
       width = 16, height = 8, bg = "white", dpi = 600, units = "cm")

################################################################################
# Revision: Comparisons by lake:
comparison_lake_basic = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/comparison/comparison_basic_strategies_lakes.rds")
comparison_patch_leaving_lake = comparison_lake_basic %>%
  ggplot(aes(y = fct_reorder(model_name, -elpd_diff), x = elpd_diff)) +
  geom_col(color = "black", width = .6) +
  geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
  #ggtitle("Model Comparison Spot-Selection") + 
  xlab(expression(Delta~"ELPD (rel. to best)")) + 
  ylab("Model") + 
  #labs(tag = "a")+
  facet_wrap(~lake,nrow = 2)+
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
comparison_patch_leaving_lake

ggsave(comparison_patch_leaving_lake, file = "4_incremental_strategies/output/comparison_patch_leaving_lake.png", 
       width = 18, height = 11, bg = "white", dpi = 600, units = "cm")
ggsave(comparison_patch_leaving_lake, file = "4_incremental_strategies/output/comparison_patch_leaving_lake.svg", 
       width = 18, height = 11, bg = "white", dpi = 600, units = "cm")
ggsave(comparison_patch_leaving_lake, file = "4_incremental_strategies/output/comparison_patch_leaving_lake.pdf", 
       width = 18, height = 11, bg = "white", dpi = 600, units = "cm")
ggsave(comparison_patch_leaving_lake, file = "4_incremental_strategies/output/comparison_patch_leaving_lake.tiff", 
       width = 18, height = 11, bg = "white", dpi = 600, units = "cm")

comparison_lake_extended = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/comparison/comparison_extended_strategies_lakes.rds")
comparison_patch_leaving_lake_extended = comparison_lake_extended %>%
  ggplot(aes(y = fct_reorder(model_name, -elpd_diff), x = elpd_diff)) +
  geom_col(color = "black", width = .6) +
  geom_errorbarh(aes(xmin = elpd_diff - 4*se_diff, xmax = elpd_diff + 4*se_diff), height = .3, size = .8) + 
  #ggtitle("Model Comparison Spot-Selection") + 
  xlab(expression(Delta~"ELPD (rel. to best)")) + 
  ylab("Model") + 
  #labs(tag = "a")+
  facet_wrap(~lake,nrow = 2)+
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
comparison_patch_leaving_lake_extended

ggsave(comparison_patch_leaving_lake_extended, file = "4_incremental_strategies/output/comparison_patch_leaving_lake_extended.png", 
       width = 18, height = 11, bg = "white", dpi = 600, units = "cm")
ggsave(comparison_patch_leaving_lake_extended, file = "4_incremental_strategies/output/comparison_patch_leaving_lake_extended.svg", 
       width = 18, height = 11, bg = "white", dpi = 600, units = "cm")
ggsave(comparison_patch_leaving_lake_extended, file = "4_incremental_strategies/output/comparison_patch_leaving_lake_extended.pdf", 
       width = 18, height = 11, bg = "white", dpi = 600, units = "cm")
ggsave(comparison_patch_leaving_lake_extended, file = "4_incremental_strategies/output/comparison_patch_leaving_lake_extended.tiff", 
       width = 18, height = 11, bg = "white", dpi = 600, units = "cm")

################################################################################
# visualize model predictions
################################################################################
# simulate rewards
rewards = rep(0, 100)
t = c(10, 20, 30, 50, 70, 85)
rewards[t]= 1
time = 1:100
cumulative_catch=cumsum(rewards)
time_since_event = vector()
time_since_event[1] = 0
for (i in 2:100){
  time_since_event[i] = ifelse(rewards[i] == 0, time_since_event[i-1]+1, 0) 
}

# load heuristic model
binary_gut_fit = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/c_patch_leaving_model_patch_discovery_fit.rds")
#print(binary_gut_fit, max_rows = 7)
draws = binary_gut_fit$draws(format = "df")
before_catch = (draws$`beta[1]`)
tss = (draws$`beta[3]`)
after_catch = (draws$`beta[1]` + draws$`beta[2]`)
tsc = (draws$`beta[3]`)

df_binary_gut = data.frame(rewards = rep(rewards, each = length(tsc)),
                           time = rep(time, each = length(tsc)),
                           time_since_event = rep(time_since_event, each = length(tsc)),
                           cumulative_catch = rep(cumulative_catch, each = length(tsc)),
                           before_catch = rep(before_catch, length(time)),
                           tss = rep(tss, length(time)),
                           after_catch = rep(after_catch, length(time)),
                           tsc = rep(tsc, length(time)))

# compute p
df_binary_gut$p = ifelse(df_binary_gut$cumulative_catch == 0,
                         inv_logit(df_binary_gut$before_catch + df_binary_gut$tss * df_binary_gut$time_since_event),
                         inv_logit(df_binary_gut$after_catch + df_binary_gut$tsc * df_binary_gut$time_since_event))


binary_gut_prediction = df_binary_gut %>% 
  ggplot(aes(x = time*10/60, y = p)) + 
  stat_lineribbon(.width = c(.66, .95), size = .6) + 
  labs(tag = "B")+
  ylab("Posterior Pred. (95%CrI)\nP(Leave)")+
  xlab("Angling Time (min)")+
  ggtitle("Model Predictions")+
  geom_vline(xintercept = (which(rewards==1)-.5)*10/60, linetype = "dashed", color = "darkgrey")+
  coord_cartesian(ylim = c(0, .15)) + 
  scale_x_continuous(breaks = seq(from = 0, to = 20, by = 2))+
  #scale_fill_brewer(palette = 2) +
  scale_fill_manual(values = c("#BAE4B3","#31A354"))+#brewer.pal(5, "Greens")[])+
  annotate("segment", x = 0, y = inv_logit(mean(before_catch)), xend = which(rewards == 1)[1]*10/60, yend = inv_logit(mean(before_catch)), size = 1.2, col = "chartreuse") +
  annotate("segment", x = which(rewards == 1)[1]*10/60, y = inv_logit(mean(after_catch)), xend = 100*10/60, yend = inv_logit(mean(after_catch)), size = 1.2, col = "chartreuse") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        plot.tag = element_text(size = 7, face = "bold"),
        title = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"),
        text = element_text(family = "sans")) + 
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = (which(rewards==1)-.5)[1]*10/60 - 2*10/60,xmax = which(rewards==1)[1]*10/60 + 2*10/60,
              ymin = .12, fill = "cadetblue",alpha = 0.9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = (which(rewards==1)-.5)[2]*10/60 - 2*10/60,xmax = which(rewards==1)[2]*10/60 + 2*10/60,
              ymin = .12, fill = "cadetblue",alpha = 0.9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = (which(rewards==1)-.5)[3]*10/60 - 2*10/60,xmax = which(rewards==1)[3]*10/60 + 2*10/60,
              ymin = .12, fill = "cadetblue",alpha = 0.9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = (which(rewards==1)-.5)[4]*10/60 - 2*10/60,xmax = which(rewards==1)[4]*10/60 + 2*10/60,
              ymin = .12, fill = "cadetblue",alpha = 0.9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = (which(rewards==1)-.5)[5]*10/60 - 2*10/60,xmax = which(rewards==1)[5]*10/60 + 2*10/60,
              ymin = .12, fill = "cadetblue",alpha = 0.9)+
  add_fishape(family = "Percidae",option = "Perca_flavescens",
              xmin = (which(rewards==1)-.5)[6]*10/60 - 2*10/60,xmax = which(rewards==1)[6]*10/60 + 2*10/60,
              ymin = .12, fill = "cadetblue",alpha = 0.9)
binary_gut_prediction

# make panel
panel_tmp_results = arrangeGrob(
  grobs = list(catch_bias, 
               binary_gut_prediction,
               cue_weights,
               spatial_weights),
  widths = c(1, 1, 1),
  heights = c(1.2, 1),
  layout_matrix = rbind(c(1, 3, 4),
                        c(2, 2, 2))
)


ggsave(panel_tmp_results, file = "4_incremental_strategies/output/panel_tmp_results.png", 
       width = 18,height = 11, bg = "white", dpi = 600, units = "cm")
ggsave(panel_tmp_results, file = "4_incremental_strategies/output/panel_tmp_results.svg", 
       width = 18,height = 11, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
ggsave(panel_tmp_results, file = "4_incremental_strategies/output/panel_tmp_results.pdf", 
       width = 18,height = 11, bg = "white", dpi = 600, units = "cm")
ggsave(panel_tmp_results, file = "4_incremental_strategies/output/panel_tmp_results.tiff", 
       width = 18,height = 11, bg = "white", dpi = 600, units = "cm")

################################################################################
# REVISION: Make figure illustrating generalization for both tasks
lambda_patch_leaving = as.vector(patch_leaving_model_fit$draws("lambda"))
spot_selection_fit = readRDS("utils/data/processed_data/model_fit/spot_selection_models/f_conditional.rds")
lambda_spot_selection = as.vector(spot_selection_fit$draws("lambdas[1]"))
lambda_df = data.frame(lambda_patch_leaving = lambda_patch_leaving,
                       lambda_spot_selection = lambda_spot_selection)
lambda_df = lambda_df %>% 
  mutate(draw = 1:n()) %>% 
  slice(rep(1:n(), each = 500)) %>% 
  pivot_longer(1:2, names_to = "condition", values_to = "lambda") %>%
  group_by(draw, condition) %>% 
  mutate(d = seq(from = 0, to = 100, length.out = 500)) %>% 
  mutate(g = exp(-d^2 / (2 * exp(lambda)^2))) %>% 
  mutate(condition = ifelse(condition == "lambda_patch_leaving", "Spot-Leaving", "Spot-Selection"))

# make plot
generalization_plot_comparison = lambda_df %>% 
  ggplot(aes(x = d, y = g, color = condition)) + 
  stat_lineribbon(fill = "lightgrey") +
  scale_color_manual(values = c("turquoise4","cadetblue2"))+
  labs(color = "Decision", fill = "CrI")+
  xlab("Distance to Competitor (m)")+
  ylab("Social Density")+
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
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2))
ggsave(generalization_plot_comparison, file = "4_incremental_strategies/output/generalization_vis_comparison.png", 
       width = 9, height = 9, bg = "white", dpi = 600, units = "cm")
ggsave(generalization_plot_comparison, file = "4_incremental_strategies/output/generalization_vis_comparison.svg", 
       width = 9, height = 9, bg = "white", dpi = 600, units = "cm")
ggsave(generalization_plot_comparison, file = "4_incremental_strategies/output/generalization_vis_comparison.pdf", 
       width = 9, height = 9, bg = "white", dpi = 600, units = "cm")
ggsave(generalization_plot_comparison, file = "4_incremental_strategies/output/generalization_vis_comparison.tiff", 
       width = 9, height = 9, bg = "white", dpi = 600, units = "cm")

}


################################################################################
# END
################################################################################