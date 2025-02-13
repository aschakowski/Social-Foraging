################################################################################
#
# Title: Visualize giving-up time simulations
#
# Description: Visualize simulations
#
# Authors: Schakowski, A.
#
# Last updated: 02/04/2024 (DD/MM/YYYY)
#
################################################################################
visualize_gut_simulations <- function(){
  
# Load simulation results
sim_baseline = readRDS(file = "utils/data/raw_data/patch_leaving_simulations/a_sim_baseline.rds")
sim_baseline = do.call(rbind, sim_baseline)
sim_time_wo_catch = readRDS(file = "utils/data/raw_data/patch_leaving_simulations/b_sim_time_wo_catch.rds")
sim_time_wo_catch = do.call(rbind, sim_time_wo_catch)
sim_patch_discovery = readRDS(file = "utils/data/raw_data/patch_leaving_simulations/c_sim_patch_discovery.rds")
sim_patch_discovery = do.call(rbind, sim_patch_discovery)
sim_global_local = readRDS(file = "utils/data/raw_data/patch_leaving_simulations/d_sim_global_local.rds")
sim_global_local = do.call(rbind, sim_global_local)
sim_spatial_features = readRDS(file = "utils/data/raw_data/patch_leaving_simulations/e_sim_spatial_features.rds")
sim_spatial_features = do.call(rbind, sim_spatial_features)

# visualize results
fun_median = function(z) {
  data.frame(
    y = median(z),
    ymin = median(z)-2*(1.2533 * sd(z)/length(z)^2),
    ymax = median(z)+2*(1.2533 * sd(z)/length(z)^2)
  )
}

################################################################################
library(weights)
plots = list()
sim_list = list(sim_baseline,
                sim_time_wo_catch,
                sim_patch_discovery,
                sim_global_local,
                sim_spatial_features)
title = c("Baseline Model",
          "+Time w/o Catch",
          "+Patch Discovery",
          "+Global/Local Updating",
          "+Spatial Features")
for (model in 1:5){
  
  # compute errorbars
  plot_data1 = sim_list[[model]] %>% 
    mutate(unique_spot = paste(trip, spot)) %>% 
    filter(iter == 1) %>% 
    #group_by(trip) %>% 
    group_by(unique_spot) %>% 
    dplyr::summarize(gut_data_mean = mean(gut_data),
              gut_data_sd = sd(gut_data),
              n = n()) %>% 
    mutate(lowerx = gut_data_mean - 2*gut_data_sd / sqrt(n),
           upperx = gut_data_mean + 2*gut_data_sd / sqrt(n)) %>% 
    mutate(lowerx = ifelse(lowerx < 0, 0, lowerx))
  
  plot_data2 = sim_list[[model]] %>% 
    mutate(unique_spot = paste(trip, spot)) %>% 
    #group_by(trip, iter) %>% 
    group_by(unique_spot, iter) %>% 
    dplyr::summarize(gut_sim_mean = mean(gut_sim)) %>% 
    #group_by(trip) %>% 
    group_by(unique_spot) %>% 
    dplyr::summarize(lowery = quantile(gut_sim_mean, .025),
              uppery = quantile(gut_sim_mean, .975),
              gut_sim_mean = mean(gut_sim_mean))
  plot_data = left_join(plot_data1, plot_data2)
  
  # compute correlation coefficient across iterations
  correlation = sim_list[[model]] %>%
    mutate(unique_spot = paste(trip, spot)) %>% 
    group_by(unique_spot, iter) %>% 
    #group_by(trip, iter) %>% 
    dplyr::summarize(gut_data = mean(gut_data),
              gut_sim = mean(gut_sim)) %>% 
    group_by(iter) %>% 
    dplyr::summarize(cor = cor(gut_data, gut_sim)) %>% 
    dplyr::summarize(lower = quantile(cor, .025),
           upper = quantile(cor, .975),
           mean = mean(cor))
  
  plots[[model]] = ggplot(data = plot_data, aes(x = gut_data_mean, y = gut_sim_mean)) + 
    geom_point(alpha = .2, size = .6)+
    #geom_errorbar(aes(xmin = lowerx, xmax = upperx), alpha = .2)+
    geom_errorbar(aes(ymin = lowery, ymax = uppery), alpha = .1) + 
    annotate("text", x = 95, y = 137, label = paste("r = ", rd(correlation$mean, digits = 2), " [", rd(correlation$lower, digits = 2), " ,", rd(correlation$upper, digits = 2), "]", sep = ""), 
             size = 6/.pt) + 
    ggtitle(title[model])+
    labs(tag = letters[model])+
    coord_fixed(ratio = 1,
                ylim = c(0, 150),
                xlim = c(0, 150))+
    geom_abline(slope = 1, linetype = "dotted")+
    scale_x_continuous(breaks = seq(from = 0, to = 150, by = 25))+
    scale_y_continuous(breaks = seq(from = 0, to = 150, by = 25))+
    xlab("Empirical Giving Up Time\nSpot (10s)") + 
    ylab("Simulated Giving Up Time\nSpot (10s)") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  
}

################################################################################
plots_trip = list()
sim_list = list(sim_baseline,
                sim_time_wo_catch,
                sim_patch_discovery,
                sim_global_local,
                sim_spatial_features)
title = c("Baseline Model",
          "+Time w/o Catch",
          "+Patch Discovery",
          "+Global/Local Updating",
          "+Spatial Features")
for (model in 1:5){
  
  # compute errorbars
  plot_data1 = sim_list[[model]] %>% 
    mutate(unique_spot = paste(trip, spot)) %>% 
    filter(iter == 1) %>% 
    group_by(trip) %>% 
    #group_by(unique_spot) %>% 
    dplyr::summarize(gut_data_mean = mean(gut_data),
                     gut_data_sd = sd(gut_data),
                     n = n()) %>% 
    mutate(lowerx = gut_data_mean - 2*gut_data_sd / sqrt(n),
           upperx = gut_data_mean + 2*gut_data_sd / sqrt(n)) %>% 
    mutate(lowerx = ifelse(lowerx < 0, 0, lowerx))
  
  plot_data2 = sim_list[[model]] %>% 
    mutate(unique_spot = paste(trip, spot)) %>% 
    group_by(trip, iter) %>% 
    #group_by(unique_spot, iter) %>% 
    dplyr::summarize(gut_sim_mean = mean(gut_sim)) %>% 
    group_by(trip) %>% 
    #group_by(unique_spot) %>% 
    dplyr::summarize(lowery = quantile(gut_sim_mean, .025),
                     uppery = quantile(gut_sim_mean, .975),
                     gut_sim_mean = mean(gut_sim_mean))
  plot_data = left_join(plot_data1, plot_data2)
  
  # compute correlation coefficient across iterations
  correlation = sim_list[[model]] %>%
    mutate(unique_spot = paste(trip, spot)) %>% 
    #group_by(unique_spot, iter) %>% 
    group_by(trip, iter) %>% 
    dplyr::summarize(gut_data = mean(gut_data),
                     gut_sim = mean(gut_sim)) %>% 
    group_by(iter) %>% 
    dplyr::summarize(cor = cor(gut_data, gut_sim)) %>% 
    dplyr::summarize(lower = quantile(cor, .025),
                     upper = quantile(cor, .975),
                     mean = mean(cor))
  
  plots_trip[[model]] = ggplot(data = plot_data, aes(x = gut_data_mean, y = gut_sim_mean)) + 
    geom_point(alpha = .4)+
    geom_errorbar(aes(xmin = lowerx, xmax = upperx), alpha = .2)+
    geom_errorbar(aes(ymin = lowery, ymax = uppery), alpha = .2) + 
    annotate("text", x = 95, y = 137, label = paste("r = ", rd(correlation$mean, digits = 2), " [", rd(correlation$lower, digits = 2), " ,", rd(correlation$upper, digits = 2), "]", sep = ""), 
             size = 6/.pt) + 
    ggtitle(title[model])+
    coord_fixed(ratio = 1,
                ylim = c(0, 150),
                xlim = c(0, 150))+
    labs(tag = letters[model])+
    geom_abline(slope = 1, linetype = "dotted")+
    scale_x_continuous(breaks = seq(from = 0, to = 150, by = 25))+
    scale_y_continuous(breaks = seq(from = 0, to = 150, by = 25))+
    xlab("Empirical Giving Up Time\nForaging Trip (10s)") + 
    ylab("Simulated Giving Up Time\nForaging Trip (10s)") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  
}

################################################################################plots = list()
plots_id = list()
sim_list = list(sim_baseline,
                sim_time_wo_catch,
                sim_patch_discovery,
                sim_global_local,
                sim_spatial_features)
title = c("Baseline Model",
          "+Time w/o Catch",
          "+Patch Discovery",
          "+Global/Local Updating",
          "+Spatial Features")
for (model in 1:5){
  
  # compute errorbars
  plot_data1 = sim_list[[model]] %>% 
    mutate(unique_spot = paste(trip, spot)) %>% 
    filter(iter == 1) %>% 
    group_by(id) %>% 
    #group_by(unique_spot) %>% 
    dplyr::summarize(gut_data_mean = mean(gut_data),
                     gut_data_sd = sd(gut_data),
                     n = n()) %>% 
    mutate(lowerx = gut_data_mean - 2*gut_data_sd / sqrt(n),
           upperx = gut_data_mean + 2*gut_data_sd / sqrt(n)) %>% 
    mutate(lowerx = ifelse(lowerx < 0, 0, lowerx))
  
  plot_data2 = sim_list[[model]] %>% 
    mutate(unique_spot = paste(trip, spot)) %>% 
    group_by(id, iter) %>% 
    #group_by(unique_spot, iter) %>% 
    dplyr::summarize(gut_sim_mean = mean(gut_sim)) %>% 
    group_by(id) %>% 
    #group_by(unique_spot) %>% 
    dplyr::summarize(lowery = quantile(gut_sim_mean, .025),
                     uppery = quantile(gut_sim_mean, .975),
                     gut_sim_mean = mean(gut_sim_mean))
  plot_data = left_join(plot_data1, plot_data2)
  
  # compute correlation coefficient across iterations
  correlation = sim_list[[model]] %>%
    mutate(unique_spot = paste(trip, spot)) %>% 
    #group_by(unique_spot, iter) %>% 
    group_by(id, iter) %>% 
    dplyr::summarize(gut_data = mean(gut_data),
                     gut_sim = mean(gut_sim)) %>% 
    group_by(iter) %>% 
    dplyr::summarize(cor = cor(gut_data, gut_sim)) %>% 
    dplyr::summarize(lower = quantile(cor, .025),
                     upper = quantile(cor, .975),
                     mean = mean(cor))
  
  plots_id[[model]] = ggplot(data = plot_data, aes(x = gut_data_mean, y = gut_sim_mean)) + 
    geom_point(alpha = .4)+
    geom_errorbar(aes(xmin = lowerx, xmax = upperx), alpha = .2)+
    geom_errorbar(aes(ymin = lowery, ymax = uppery), alpha = .2) + 
    annotate("text", x = 95, y = 137, label = paste("r = ", rd(correlation$mean, digits = 2), " [", rd(correlation$lower, digits = 2), " ,", rd(correlation$upper, digits = 2), "]", sep = ""), 
             size = 6/.pt) + 
    ggtitle(title[model])+
    labs(tag = letters[model])+
    coord_fixed(ratio = 1,
                ylim = c(0, 150),
                xlim = c(0, 150))+
    geom_abline(slope = 1, linetype = "dotted")+
    scale_x_continuous(breaks = seq(from = 0, to = 150, by = 25))+
    scale_y_continuous(breaks = seq(from = 0, to = 150, by = 25))+
    xlab("Empirical Giving Up Time\nID (10s)") + 
    ylab("Simulated Giving Up Time\nID (10s)") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  
}

################################################################################
plots_lake = list()
sim_list = list(sim_baseline,
                sim_time_wo_catch,
                sim_patch_discovery,
                sim_global_local,
                sim_spatial_features)
title = c("Baseline Model",
          "+Time w/o Catch",
          "+Patch Discovery",
          "+Global/Local Updating",
          "+Spatial Features")
for (model in 1:5){
  
  # compute errorbars
  plot_data1 = sim_list[[model]] %>% 
    mutate(unique_spot = paste(trip, spot)) %>% 
    filter(iter == 1) %>% 
    group_by(lake) %>% 
    #group_by(unique_spot) %>% 
    dplyr::summarize(gut_data_mean = mean(gut_data),
                     gut_data_sd = sd(gut_data),
                     n = n()) %>% 
    mutate(lowerx = gut_data_mean - 2*gut_data_sd / sqrt(n),
           upperx = gut_data_mean + 2*gut_data_sd / sqrt(n)) %>% 
    mutate(lowerx = ifelse(lowerx < 0, 0, lowerx))
  
  plot_data2 = sim_list[[model]] %>% 
    mutate(unique_spot = paste(trip, spot)) %>% 
    group_by(lake, iter) %>% 
    #group_by(unique_spot, iter) %>% 
    dplyr::summarize(gut_sim_mean = mean(gut_sim)) %>% 
    group_by(lake) %>% 
    #group_by(unique_spot) %>% 
    dplyr::summarize(lowery = quantile(gut_sim_mean, .025),
                     uppery = quantile(gut_sim_mean, .975),
                     gut_sim_mean = mean(gut_sim_mean))
  plot_data = left_join(plot_data1, plot_data2)
  
  # compute correlation coefficient across iterations
  correlation = sim_list[[model]] %>%
    mutate(unique_spot = paste(trip, spot)) %>% 
    #group_by(unique_spot, iter) %>% 
    group_by(lake, iter) %>% 
    dplyr::summarize(gut_data = mean(gut_data),
                     gut_sim = mean(gut_sim)) %>% 
    group_by(iter) %>% 
    dplyr::summarize(cor = cor(gut_data, gut_sim)) %>% 
    dplyr::summarize(lower = quantile(cor, .025),
                     upper = quantile(cor, .975),
                     mean = mean(cor))
  
  plots_lake[[model]] = ggplot(data = plot_data, aes(x = gut_data_mean, y = gut_sim_mean)) + 
    geom_point(alpha = .4)+
    geom_errorbar(aes(xmin = lowerx, xmax = upperx), alpha = .2)+
    geom_errorbar(aes(ymin = lowery, ymax = uppery), alpha = .2) + 
    annotate("text", x = 23, y = 38, label = paste("r = ", rd(correlation$mean, digits = 2), " [", rd(correlation$lower, digits = 2), " ,", rd(correlation$upper, digits = 2), "]", sep = ""), 
             size = 6/.pt) + 
    ggtitle(title[model])+
    labs(tag = letters[model])+
    coord_fixed(ratio = 1,
                ylim = c(0, 40),
                xlim = c(0, 40))+
    geom_abline(slope = 1, linetype = "dotted")+
    scale_x_continuous(breaks = seq(from = 0, to = 150, by = 10))+
    scale_y_continuous(breaks = seq(from = 0, to = 150, by = 10))+
    xlab("Empirical Giving Up Time\nLake (10s)") + 
    ylab("Simulated Giving Up Time\nLake (10s)") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.9, .2),
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  
}

################################################################################
# make panel
# make panel with model predictions
panel_spot_predictions = arrangeGrob(
  grobs = list(plots[[1]],
               plots[[2]],
               plots[[3]],
               plots[[4]],
               plots[[5]]),
  widths = c(1, 1, 1, 1, 1, 1),
  heights = c(1, 1),
  layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                        c(NA, 4, 4, 5, 5, NA))
)
panel_trip_predictions = arrangeGrob(
  grobs = list(plots_trip[[1]],
               plots_trip[[2]],
               plots_trip[[3]],
               plots_trip[[4]],
               plots_trip[[5]]),
  widths = c(1, 1, 1, 1, 1, 1),
  heights = c(1, 1),
  layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                        c(NA, 4, 4, 5, 5, NA))
)
panel_id_predictions = arrangeGrob(
  grobs = list(plots_id[[1]],
               plots_id[[2]],
               plots_id[[3]],
               plots_id[[4]],
               plots_id[[5]]),
  widths = c(1, 1, 1, 1, 1, 1),
  heights = c(1, 1),
  layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                        c(NA, 4, 4, 5, 5, NA))
)
panel_lake_predictions = arrangeGrob(
  grobs = list(plots_lake[[1]],
               plots_lake[[2]],
               plots_lake[[3]],
               plots_lake[[4]],
               plots_lake[[5]]),
  widths = c(1, 1, 1, 1, 1, 1),
  heights = c(1, 1),
  layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                        c(NA, 4, 4, 5, 5, NA))
)

ggsave(panel_spot_predictions, file = "4_incremental_strategies/output/panel_patch_leaving_spot_predictions.png", 
       width = 18, height = 13, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_trip_predictions, file = "4_incremental_strategies/output/panel_patch_leaving_trip_predictions.png", 
       width = 18, height = 13, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_id_predictions, file = "4_incremental_strategies/output/panel_patch_leaving_id_predictions.png", 
       width = 18, height = 13, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_lake_predictions, file = "4_incremental_strategies/output/panel_patch_leaving_lake_predictions.png", 
       width = 18, height = 13, bg = "white", dpi = 1200, units = "cm")

################################################################################
# visualize gut by success/no success
plots_spot_success = list()
sim_list = list(sim_baseline,
                sim_time_wo_catch,
                sim_patch_discovery,
                sim_global_local,
                sim_spatial_features)
title = c("Baseline Model",
          "+Time w/o Catch",
          "+Patch Discovery",
          "+Global/Local Updating",
          "+Spatial Features")
for (model in 1:5){
  
  # compute errorbars
  plot_data1 = sim_list[[model]] %>% 
    mutate(unique_spot = paste(trip, spot)) %>% 
    filter(iter == 1) %>% 
    #group_by(trip) %>% 
    group_by(unique_spot) %>% 
    dplyr::summarize(gut_data_mean = mean(gut_data),
                     gut_data_sd = sd(gut_data),
                     n = n()) %>% 
    mutate(lowerx = gut_data_mean - 2*gut_data_sd / sqrt(n),
           upperx = gut_data_mean + 2*gut_data_sd / sqrt(n)) %>% 
    mutate(lowerx = ifelse(lowerx < 0, 0, lowerx))
  
  plot_data2 = sim_list[[model]] %>% 
    mutate(unique_spot = paste(trip, spot)) %>% 
    #group_by(trip, iter) %>% 
    group_by(unique_spot, iter) %>% 
    dplyr::summarize(gut_sim_mean = mean(gut_sim),
                     catch = as.numeric(sum(catch)>0),
                     catch_sim = as.numeric(sum(catch_sim)>0)) %>% 
    #group_by(trip) %>% 
    group_by(unique_spot) %>% 
    dplyr::summarize(lowery = quantile(gut_sim_mean, .025),
                     uppery = quantile(gut_sim_mean, .975),
                     gut_sim_mean = mean(gut_sim_mean),
                     catch = mean(catch),
                     catch_sim = mean(catch_sim))
  plot_data = left_join(plot_data1, plot_data2)
  
  # compute correlation coefficient across iterations split by success
  correlation = sim_list[[model]] %>%
    mutate(unique_spot = paste(trip, spot)) %>%
    group_by(unique_spot, iter) %>% 
    mutate(catch = as.numeric(catch > 0),
           catch_sim = as.numeric(catch_sim > 0)) %>% 
    filter(catch == catch_sim) %>% 
    #group_by(trip, iter) %>% 
    dplyr::summarize(gut_data = mean(gut_data),
                     gut_sim = mean(gut_sim),
                     catch = mean(catch)) %>% 
    group_by(iter, catch) %>% 
    dplyr::summarize(cor = cor(gut_data, gut_sim)) %>% 
    group_by(catch) %>% 
    dplyr::summarize(lower = quantile(cor, .025),
                     upper = quantile(cor, .975),
                     mean = mean(cor))
  
  plot_data_sum = plot_data %>% 
    group_by(catch) %>% 
    dplyr::summarize(gut_data_mean_mean = mean(gut_data_mean),
                     lower_data = mean(gut_data_mean) - 2 * sd(gut_data_mean)/sqrt(n()),
                     upper_data = mean(gut_data_mean) + 2 * sd(gut_data_mean)/sqrt(n()))
  
  plot_sim_sum = sim_list[[model]] %>%
    mutate(unique_spot = paste(trip, spot)) %>%
    group_by(unique_spot, iter) %>% 
    mutate(catch = as.numeric(catch > 0),
           catch_sim = as.numeric(catch_sim > 0)) %>% 
    filter(catch == catch_sim) %>% 
    group_by(catch,iter) %>% 
    dplyr::summarize(gut_sim_mean_mean = mean(gut_sim)) %>% 
    group_by(catch) %>% 
    dplyr::summarize(gut_sim_mean = mean(gut_sim_mean_mean),
                     lower_sim = quantile(gut_sim_mean_mean, .025),
                     upper_sim = quantile(gut_sim_mean_mean, .975))
  plot_sum = left_join(plot_data_sum, plot_sim_sum)
  
  colors <- c("Successful Spots\n(N Catch > 0)" = "cyan3", "Unsuccessful Spots\n(N Catch = 0)" = "darkgreen")
  if (model == 5){
  plots_spot_success[[model]] = ggplot() + 
    
    geom_errorbar(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, ymin = lower_sim, ymax = upper_sim), alpha = 1, color = "black", width = .2) + 
    geom_errorbar(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, ymin = lower_sim, ymax = upper_sim), alpha = 1, color = "black", width = .2) + 
    geom_errorbar(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, xmin = lower_data, xmax = upper_data), alpha = 1, color = "black", width = .2) + 
    geom_errorbar(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, xmin = lower_data, xmax = upper_data), alpha = 1, color = "black", width = .2) + 
    geom_point(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, color = "Unsuccessful Spots\n(N Catch = 0)"), alpha = 1)+
    geom_point(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, color = "Successful Spots\n(N Catch > 0)"), alpha = 1)+
      
    annotate("text", color = "darkgreen", x = plot_sum$gut_data_mean_mean[plot_sum$catch == 0], y = plot_sum$gut_sim_mean[plot_sum$catch == 0] + 2, label = paste("r = ", rd(correlation$mean[correlation$catch == 0], digits = 2), " [", rd(correlation$lower[correlation$catch == 0], digits = 2), " ,", rd(correlation$upper[correlation$catch == 0], digits = 2), "]", sep = ""), 
             size = 6/.pt) + 
    annotate("text", color = "cyan3", x = plot_sum$gut_data_mean_mean[plot_sum$catch == 1], y = plot_sum$gut_sim_mean[plot_sum$catch == 1] + 2, label = paste("r = ", rd(correlation$mean[correlation$catch == 1], digits = 2), " [", rd(correlation$lower[correlation$catch == 1], digits = 2), " ,", rd(correlation$upper[correlation$catch == 1], digits = 2), "]", sep = ""), 
             size = 6/.pt) + 
    scale_color_manual(values = colors)+
    labs(color = "", tag = letters[[model]])+
    ggtitle(title[model])+
    coord_fixed(ratio = 1,
                ylim = c(5, 25),
                xlim = c(5, 25))+
    geom_abline(slope = 1, linetype = "dotted")+
    scale_x_continuous(breaks = seq(from = 0, to = 25, by = 5))+
    scale_y_continuous(breaks = seq(from = 0, to = 25, by = 5))+
    xlab("Empirical Giving Up Time\nSpot (10s)") + 
    ylab("Simulated Giving Up Time\nSpot (10s)") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_text(size = 6),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "right",
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 5),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  } else {
    plots_spot_success[[model]] = ggplot() + 
      
      geom_errorbar(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, ymin = lower_sim, ymax = upper_sim), alpha = 1, color = "black", width = .2) + 
      geom_errorbar(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, ymin = lower_sim, ymax = upper_sim), alpha = 1, color = "black", width = .2) + 
      geom_errorbar(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, xmin = lower_data, xmax = upper_data), alpha = 1, color = "black", width = .2) + 
      geom_errorbar(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, xmin = lower_data, xmax = upper_data), alpha = 1, color = "black", width = .2) + 
      geom_point(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, color = "Unsuccessful Spots\n(N Catch = 0)"), alpha = 1)+
      geom_point(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, color = "Successful Spots\n(N Catch > 0)"), alpha = 1)+
      
      annotate("text", color = "darkgreen", x = plot_sum$gut_data_mean_mean[plot_sum$catch == 0], y = plot_sum$gut_sim_mean[plot_sum$catch == 0] + 2, label = paste("r = ", rd(correlation$mean[correlation$catch == 0], digits = 2), " [", rd(correlation$lower[correlation$catch == 0], digits = 2), " ,", rd(correlation$upper[correlation$catch == 0], digits = 2), "]", sep = ""), 
               size = 6/.pt) + 
      annotate("text", color = "cyan3", x = plot_sum$gut_data_mean_mean[plot_sum$catch == 1], y = plot_sum$gut_sim_mean[plot_sum$catch == 1] + 2, label = paste("r = ", rd(correlation$mean[correlation$catch == 1], digits = 2), " [", rd(correlation$lower[correlation$catch == 1], digits = 2), " ,", rd(correlation$upper[correlation$catch == 1], digits = 2), "]", sep = ""), 
               size = 6/.pt) + 
      scale_color_manual(values = colors)+
      labs(color = "", tag = letters[[model]])+
      ggtitle(title[model])+
      coord_fixed(ratio = 1,
                  ylim = c(5, 25),
                  xlim = c(5, 25))+
      geom_abline(slope = 1, linetype = "dotted")+
      scale_x_continuous(breaks = seq(from = 0, to = 25, by = 5))+
      scale_y_continuous(breaks = seq(from = 0, to = 25, by = 5))+
      xlab("Empirical Giving Up Time\nSpot (10s)") + 
      ylab("Simulated Giving Up Time\nSpot (10s)") + 
      theme(panel.background = element_rect(fill = NA),
            axis.title.y = element_text(size = 7),
            axis.title.x = element_text(size = 7),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 6),
            panel.border = element_rect(colour = "black", fill=NA),
            legend.position = "none",
            plot.tag = element_text(size = 7, face = "bold"),
            title = element_text(size = 7),
            legend.text = element_text(size = 5),
            legend.title = element_text(size = 5),
            legend.background = element_rect(colour = "white", fill = "white"),
            text = element_text(family = "sans"))
  }
  
}

################################################################################
# model comparing avg gut at successful and non successful angling spots. 
library(dplyr)
library(brms)
library(tidybayes)
data = sim_baseline %>% 
  group_by(trip, spot) %>% 
  slice(1)
data$patch_discovery = ifelse(data$catch > 0, 1, 0)
#model_fit = brm(data = data, 
#            formula = gut_data ~ patch_discovery + (patch_discovery|id) + (patch_discovery|trip) + (patch_discovery|lake),
#            chains = 4, 
#            cores = 4)
#get_variables(model_fit)
library(ggplot2)
################################################################################

################################################################################
# visualize gut by success/no success
#plots_within_trip = list()
#sim_list = list(sim_baseline,
#                sim_time_wo_catch,
#                sim_patch_discovery,
#                sim_global_local,
#                sim_spatial_features)
#title = c("Baseline Model",
#          "+Time w/o Catch",
#          "+Patch Discovery",
#          "+Global/Local Updating",
#          "+Spatial Features")
#for (model in 1:5){
#  
#  # compute errorbars
#  plot_data1 = sim_list[[model]] %>% 
#    mutate(unique_spot = paste(trip, spot)) %>% 
#    filter(iter == 1) %>% 
#    mutate(catch = as.numeric(catch > 0 ),
#           catch_sim = as.numeric(catch_sim > 0)) %>%
#    filter(catch == catch_sim) %>% 
#    group_by(trip, catch) %>% 
#    dplyr::summarize(gut_data_mean = mean(gut_data),
#                     gut_data_sd = sd(gut_data),
#                     n = n()) %>% 
#    mutate(lowerx = gut_data_mean - 2*gut_data_sd / sqrt(n),
#           upperx = gut_data_mean + 2*gut_data_sd / sqrt(n)) %>% 
#    mutate(lowerx = ifelse(lowerx < 0, 0, lowerx))
#  
#  plot_data2 = sim_list[[model]] %>% 
#    mutate(unique_spot = paste(trip, spot)) %>% 
#    group_by(trip, iter) %>% 
#    mutate(catch = as.numeric(catch > 0 ),
#           catch_sim = as.numeric(catch_sim > 0)) %>%
#    filter(catch == catch_sim) %>% 
#    group_by(trip, iter, catch_sim) %>% 
#    dplyr::summarize(gut_sim_mean = mean(gut_sim),
#                     catch = as.numeric(sum(catch)>0),
#                     catch_sim = as.numeric(sum(catch_sim)>0)) %>% 
#    #group_by(trip) %>% 
#    group_by(trip, catch_sim) %>% 
#    dplyr::summarize(lowery = quantile(gut_sim_mean, .025),
#                     uppery = quantile(gut_sim_mean, .975),
#                     gut_sim_mean = mean(gut_sim_mean),
#                     catch = mean(catch),
#                     catch_sim = mean(catch_sim))
#  plot_data = left_join(plot_data1, plot_data2)
#  
#  # compute correlation coefficient across iterations split by success
#  correlation = sim_list[[model]] %>%
#    mutate(unique_spot = paste(trip, spot)) %>%
#    group_by(unique_spot, iter) %>% 
#    mutate(catch = as.numeric(catch > 0),
#           catch_sim = as.numeric(catch_sim > 0)) %>% 
#    filter(catch == catch_sim) %>% 
#    #group_by(trip, iter) %>% 
#    dplyr::summarize(gut_data = mean(gut_data),
#                     gut_sim = mean(gut_sim),
#                     catch = mean(catch)) %>% 
#    group_by(iter, catch) %>% 
#    dplyr::summarize(cor = cor(gut_data, gut_sim)) %>% 
#    group_by(catch) %>% 
#    dplyr::summarize(lower = quantile(cor, .025),
#                     upper = quantile(cor, .975),
#                     mean = mean(cor))
#  
#  plot_data_sum = plot_data %>% 
#    group_by(catch) %>% 
#    dplyr::summarize(gut_data_mean_mean = mean(gut_data_mean),
#                     lower_data = mean(gut_data_mean) - 2 * sd(gut_data_mean)/sqrt(n()),
#                     upper_data = mean(gut_data_mean) + 2 * sd(gut_data_mean)/sqrt(n()))
#  
#  plot_sim_sum = sim_list[[model]] %>%
#    mutate(unique_spot = paste(trip, spot)) %>%
#    group_by(unique_spot, iter) %>% 
#    mutate(catch = as.numeric(catch > 0),
#           catch_sim = as.numeric(catch_sim > 0)) %>% 
#    filter(catch == catch_sim) %>% 
#    group_by(catch,iter) %>% 
#    dplyr::summarize(gut_sim_mean_mean = mean(gut_sim)) %>% 
#    group_by(catch) %>% 
#    dplyr::summarize(gut_sim_mean = mean(gut_sim_mean_mean),
#                     lower_sim = quantile(gut_sim_mean_mean, .025),
#                     upper_sim = quantile(gut_sim_mean_mean, .975))
#  plot_sum = left_join(plot_data_sum, plot_sim_sum)
#  
#  colors <- c("Successful Spots\n(N Catch > 0)" = "cyan3", "Unsuccessful Spots\n(N Catch = 0)" = "darkgreen")
#  if (model == 5){
#    plots_spot_success[[model]] = ggplot() + 
#      
#     geom_errorbar(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, ymin = lower_sim, ymax = upper_sim), alpha = 1, color = "black", width = .2, size = .8) + 
#      geom_errorbar(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, ymin = lower_sim, ymax = upper_sim), alpha = 1, color = "black", width = .2, size = .8) + 
#      geom_errorbar(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, xmin = lower_data, xmax = upper_data), alpha = 1, color = "black", width = .2, size = .8) + 
#      geom_errorbar(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, xmin = lower_data, xmax = upper_data), alpha = 1, color = "black", width = .2, size = .8) + 
#      geom_point(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, color = "Unsuccessful Spots\n(N Catch = 0)"), alpha = 1, size = 3)+
#      geom_point(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, color = "Successful Spots\n(N Catch > 0)"), alpha = 1, size = 3)+
#      
#      annotate("text", color = "darkgreen", x = plot_sum$gut_data_mean_mean[plot_sum$catch == 0], y = plot_sum$gut_sim_mean[plot_sum$catch == 0] + 2, label = paste("r = ", rd(correlation$mean[correlation$catch == 0], digits = 2), " [", rd(correlation$lower[correlation$catch == 0], digits = 2), " ,", rd(correlation$upper[correlation$catch == 0], digits = 2), "]", sep = ""), 
#               size = 4) + 
#      annotate("text", color = "cyan3", x = plot_sum$gut_data_mean_mean[plot_sum$catch == 1], y = plot_sum$gut_sim_mean[plot_sum$catch == 1] + 2, label = paste("r = ", rd(correlation$mean[correlation$catch == 1], digits = 2), " [", rd(correlation$lower[correlation$catch == 1], digits = 2), " ,", rd(correlation$upper[correlation$catch == 1], digits = 2), "]", sep = ""), 
#               size = 4) + 
#      scale_color_manual(values = colors)+
#      labs(color = "")+
#      ggtitle(title[model])+
#      coord_fixed(ratio = 1,
#                  ylim = c(0, 25),
#                  xlim = c(0, 25))+
#      geom_abline(slope = 1, linetype = "dotted")+
#      scale_x_continuous(breaks = seq(from = 0, to = 25, by = 5))+
#      scale_y_continuous(breaks = seq(from = 0, to = 25, by = 5))+
#      xlab("Empirical Giving Up Time\nSpot [10s]") + 
#      ylab("Simulated Giving Up Time\nSpot [10s]") + 
#      theme(panel.background = element_rect(fill = NA),
#            axis.title.y = element_text(size = 18),
#            axis.title.x = element_text(size = 18),
#            axis.text.x = element_text(size = 18),
#            axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
#            panel.border = element_rect(colour = "black", fill=NA),
#            legend.position = c(.65, .2),
#            plot.tag = element_text(size = 24, face = "bold"),
#            title = element_text(size = 18),
#            legend.text = element_text(size = 14),
#            legend.title = element_text(size = 14),
#            legend.background = element_rect(colour = "white", fill = "white"),
#            text = element_text(family = "sans"))
#  } else {
#    plots_spot_success[[model]] = ggplot() + 
#      
#      geom_errorbar(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, ymin = lower_sim, ymax = upper_sim), alpha = 1, color = "black", width = .2, size = .8) + 
#      geom_errorbar(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, ymin = lower_sim, ymax = upper_sim), alpha = 1, color = "black", width = .2, size = .8) + 
#      geom_errorbar(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, xmin = lower_data, xmax = upper_data), alpha = 1, color = "black", width = .2, size = .8) + 
#      geom_errorbar(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, xmin = lower_data, xmax = upper_data), alpha = 1, color = "black", width = .2, size = .8) + 
#      geom_point(data = plot_sum %>% filter(catch == 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, color = "Unsuccessful Spots\n(N Catch = 0)"), alpha = 1, size = 3)+
#      geom_point(data = plot_sum %>% filter(catch > 0), aes(x = gut_data_mean_mean, y = gut_sim_mean, color = "Successful Spots\n(N Catch > 0)"), alpha = 1, size = 3)+
#      
#      annotate("text", color = "darkgreen", x = plot_sum$gut_data_mean_mean[plot_sum$catch == 0], y = plot_sum$gut_sim_mean[plot_sum$catch == 0] + 2, label = paste("r = ", rd(correlation$mean[correlation$catch == 0], digits = 2), " [", rd(correlation$lower[correlation$catch == 0], digits = 2), " ,", rd(correlation$upper[correlation$catch == 0], digits = 2), "]", sep = ""), 
#               size = 4) + 
#      annotate("text", color = "cyan3", x = plot_sum$gut_data_mean_mean[plot_sum$catch == 1], y = plot_sum$gut_sim_mean[plot_sum$catch == 1] + 2, label = paste("r = ", rd(correlation$mean[correlation$catch == 1], digits = 2), " [", rd(correlation$lower[correlation$catch == 1], digits = 2), " ,", rd(correlation$upper[correlation$catch == 1], digits = 2), "]", sep = ""), 
#               size = 4) + 
#      scale_color_manual(values = colors)+
#      labs(color = "")+
#      ggtitle(title[model])+
#      coord_fixed(ratio = 1,
#                  ylim = c(0, 25),
#                  xlim = c(0, 25))+
#      geom_abline(slope = 1, linetype = "dotted")+
#      scale_x_continuous(breaks = seq(from = 0, to = 25, by = 5))+
#      scale_y_continuous(breaks = seq(from = 0, to = 25, by = 5))+
#      xlab("Empirical Giving Up Time\nSpot [10s]") + 
#      ylab("Simulated Giving Up Time\nSpot [10s]") + 
#      theme(panel.background = element_rect(fill = NA),
#            axis.title.y = element_text(size = 18),
#            axis.title.x = element_text(size = 18),
#            axis.text.x = element_text(size = 18),
#            axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
#            panel.border = element_rect(colour = "black", fill=NA),
#            legend.position = "none",
#            plot.tag = element_text(size = 24, face = "bold"),
#            title = element_text(size = 18),
#            legend.text = element_text(size = 18),
#            legend.title = element_text(size = 18),
#            legend.background = element_rect(colour = "white", fill = "white"),
#            text = element_text(family = "sans"))
#  }
#  
#}

################################################################################
# make panels
panel_delta_gut_predictions = arrangeGrob(
  grobs = list(plots_spot_success[[1]],
               plots_spot_success[[2]],
               plots_spot_success[[3]],
               plots_spot_success[[4]],
               plots_spot_success[[5]]),
  widths = c(1, 1, 1, 1, 1, 1, .055),
  heights = c(1, 1),
  layout_matrix = rbind(c(1, 1, 2, 2, 3, 3, NA),
                        c(NA, 4, 4, 5, 5, 5, 5))
)

ggsave(panel_delta_gut_predictions, file = "4_incremental_strategies/output/panel_patch_leaving_delta_gut_predictions.png", 
       width = 18, height = 13, bg = "white", dpi = 1200, units = "cm")

}
################################################################################
# END
################################################################################