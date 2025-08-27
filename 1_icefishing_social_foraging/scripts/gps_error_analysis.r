################################################################################
# GPS error analysis
################################################################################
visualize_gps_error <- function(){
  
  # load gps data
  data = readRDS("utils/data/raw_data/reliability_data/gps_accuracy_data.rds")
  
  # now compute pairwise distances
  watch_1_glonass = data %>% arrange(date) %>% filter(setting == "glonass" & watch_id == 1) 
  watch_2_glonass = data %>% arrange(date) %>% filter(setting == "glonass" & watch_id == 2) 
  watch_3_glonass = data %>% arrange(date) %>% filter(setting == "glonass" & watch_id == 3) 
  watch_4_glonass = data %>% arrange(date) %>% filter(setting == "glonass" & watch_id == 4) 
  
  # Step 1: Plot all trajectories at once
  #data %>% 
  #  ggplot(aes(x = position_long, y = position_lat)) + 
  #  geom_point()+
  #  xlab("Longitude")+ylab("Latitude")+
  #  ggtitle("Accuracy Test")
  
  
  # pairwise distances
  pw_dist_fun <- function(d1, d2, n){
    
    res = vector()
    
    for (i in 1:n){
      res[i] = distm(x = c(d1$position_long[i], d1$position_lat[i]),
                     y = c(d2$position_long[i], d2$position_lat[i]),
                     fun = distHaversine)
    }
    
    return(res)
    
  }
  
  n = min(nrow(watch_1_glonass),
          nrow(watch_2_glonass),
          nrow(watch_3_glonass),
          nrow(watch_4_glonass))
  
  # comp 1
  distances12_glonass = pw_dist_fun(watch_1_glonass, watch_2_glonass, n = n)
  
  # comparison 2
  distances13_glonass = pw_dist_fun(watch_1_glonass, watch_3_glonass, n = n)
  
  # comparison 3
  distances14_glonass = pw_dist_fun(watch_1_glonass, watch_4_glonass, n = n)
  
  # comparison 4
  distances23_glonass = pw_dist_fun(watch_2_glonass, watch_3_glonass, n = n)
  
  # comparison 5
  distances24_glonass = pw_dist_fun(watch_2_glonass, watch_4_glonass, n = n)
  
  # comparison 6
  distances34_glonass = pw_dist_fun(watch_3_glonass, watch_4_glonass, n = n)
  
  # plot this: first create one big dataframe from this
  acc_data = rbind(cbind(distances12_glonass,
                         distances13_glonass,
                         distances14_glonass,
                         distances23_glonass,
                         distances24_glonass,
                         distances34_glonass,
                         "glonass"))
  colnames(acc_data) = c("1 vs 2", "1 vs 3", "1 vs 4",
                         "2 vs 3", "2 vs 4", "3 vs 4", "setting")
  acc_data = reshape2::melt(as.data.frame(acc_data), id.vars = "setting", variable.name = "comparison", value.name = "distance")
  acc_data$distance = as.numeric(acc_data$distance)
  
  ################################################################################
  # visualize
  ################################################################################
  histogram_errors = ggplot(data = acc_data %>% filter(setting == "glonass"), aes(x = distance)) + 
    geom_histogram(alpha = .8, position = "dodge", color = "black")+
    scale_fill_brewer()+
    xlab("Pairwise distance (m)") + 
    labs(tag = "A")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .8),
          strip.text.y = element_blank(), 
          strip.background = element_blank(),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "black", fill = "white"))
  histogram_errors
  
  boxplot_errors = ggplot(data = acc_data %>% filter(setting == "glonass"), aes(x = comparison, y = distance)) + 
    geom_violin(alpha = .2)+
    geom_boxplot(width=.1, size = .5, outlier.size = .2, outlier.alpha = .2)+
    geom_hline(yintercept = 0, linetype = "dotted") +
    xlab("comparison") + 
    labs(tag = "B")+
    ylab("Pairwise distance (m)") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          strip.text.y = element_blank(), 
          strip.background = element_blank(),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "black", fill = "white"))
  boxplot_errors
  
  # summary for text
  acc_data %>% group_by(setting, comparison) %>% dplyr::summarize(median_distance = median(distance, na.rm = T),
                                                           min_distance = min(distance, na.rm = T),
                                                           max_distance = max(distance, na.rm = T))
  
  
  # compute quantiles for text
  quantile(acc_data %>% filter(setting != "galileo") %>% dplyr::select(distance), na.rm = T, p = seq(from = 0, to = 1, by = .05))
  
  ################################################################################
  # compute stationary gps error
  ################################################################################
  # load data
  acc_data_angling = readRDS(file = "utils/data/raw_data/reliability_data/gps_accuracy_data_stationary.rds")
  
  ################################################################################
  # visualize
  ################################################################################
  gps_error_angling_histogram = 
    ggplot(data = acc_data_angling, aes(x = error)) + 
    geom_histogram(alpha = .8, position = "dodge", color = "black")+
    
    scale_fill_brewer()+
    labs(tag = "C")+
    scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    xlab("Distance from center (m)") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .8),
          strip.text.y = element_blank(), 
          strip.background = element_blank(),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "black", fill = "white"))
  gps_error_angling_histogram
  
  gps_error_angling_mean_histogram = 
    ggplot(data = acc_data_angling %>% group_by(day, year, watch_id, sequence) %>% slice(1), aes(x = mean_error)) + 
    geom_histogram(alpha = .8, position = "dodge", color = "black")+
    
    scale_fill_brewer()+
    labs(tag = "E")+
    #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    xlab("Avg. distance from center (m)") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .8),
          strip.text.y = element_blank(), 
          strip.background = element_blank(),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "black", fill = "white"))
  gps_error_angling_mean_histogram
  
  gps_error_angling_boxplot = ggplot(data = acc_data_angling, aes(x = interaction(day, year), y = error)) + 
    geom_violin()+
    geom_boxplot(width=.1, size = .5, outlier.size = .2, outlier.alpha = .2)+
    
    #scale_fill_brewer()+
    #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    xlab("Competition (d.yyyy)") + 
    labs(tag = "D")+
    ylab("Distance from center (m)") +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          strip.text.y = element_blank(), 
          strip.background = element_blank(),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "black", fill = "white"))
  gps_error_angling_boxplot
  
  gps_error_angling_mean_boxplot = ggplot(data = acc_data_angling %>% group_by(day, year, watch_id, sequence) %>% slice(1), aes(x = interaction(day, year), y = mean_error)) + 
    geom_violin()+
    geom_boxplot(width=.1, size = .5, outlier.size = .2, outlier.alpha = .2)+
    
    #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
    xlab("Competition (d.yyyy)") + 
    labs(tag = "F")+
    ylab("Avg. distance from center (m)") +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 7),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          strip.text.y = element_blank(), 
          strip.background = element_blank(),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "black", fill = "white"))
  gps_error_angling_mean_boxplot
  
  # statistics for text
  acc_data_angling %>% group_by(day, year) %>% dplyr::summarize(median_error = median(error, na.rm = T),
                                                         min_error = min(error, na.rm = T),
                                                         max_error = max(error, na.rm = T))
  
  acc_data_angling %>% group_by(day, year, watch_id, sequence) %>% slice(1) %>% ungroup() %>% dplyr::summarize(median_error = median(mean_error, na.rm = T),
                                                                                                        min_error = min(mean_error, na.rm = T),
                                                                                                        max_error = max(mean_error, na.rm = T))
  
  # compute quantiles for text
  quantile(acc_data_angling %>% ungroup() %>% dplyr::select(error), na.rm = T, p = seq(from = 0, to = 1, by = .05))
  quantile(acc_data_angling %>% group_by(day, year, watch_id, sequence) %>% slice(1) %>% ungroup() %>% dplyr::select(mean_error), na.rm = T, p = seq(from = 0, to = 1, by = .05))
  
  
  ################################################################################
  # make figure
  ################################################################################
  # make panel
  library(gridExtra)
  gps_error_panel = arrangeGrob(
    grobs = list(histogram_errors, boxplot_errors,
                 gps_error_angling_histogram, gps_error_angling_boxplot,
                 gps_error_angling_mean_histogram, gps_error_angling_mean_boxplot),
    widths = c(1, 1),
    layout_matrix = rbind(c(1, 2),
                          c(3, 4),
                          c(5, 6))
  )
  ggsave(gps_error_panel, file = "1_icefishing_social_foraging/output/panel_gps_erorr.png", 
         width = 18, height = 18,bg = "white", dpi = 600, units = "cm")
  ggsave(gps_error_panel, file = "1_icefishing_social_foraging/output/panel_gps_erorr.pdf", 
         width = 18, height = 18,bg = "white", dpi = 600, units = "cm")
  ggsave(gps_error_panel, file = "1_icefishing_social_foraging/output/panel_gps_erorr.svg", 
         width = 18, height = 18,bg = "white", dpi = 600, units = "cm")
  ggsave(gps_error_panel, file = "1_icefishing_social_foraging/output/panel_gps_erorr.tiff", 
         width = 18, height = 18,bg = "white", dpi = 600, units = "cm")
  
}