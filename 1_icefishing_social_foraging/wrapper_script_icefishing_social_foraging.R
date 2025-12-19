################################################################################
#
# Title: Descriptive Analyses and Figures
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 19/12/2025 (DD/MM/YYYY)
#
################################################################################
run_1_icefishing_social_foraging <- function(){
  
  ################################################################################
  # load data
  ################################################################################
  
  # dataset contains basic demographics, catch success for each competition and camera/watch ids
  catch_data <- read_csv("utils/data/raw_data/catch_data.csv")
  
  # data contains survey responses for all participants (2 missing) 
  survey_data = read_csv("utils/data/raw_data/survey_data.csv")
  
  # data contains information for all foraging locations (ID equivalent to camera_id in other datasets)
  spot_selection_data = fread(file = "utils/data/raw_data/spot_selection_data.csv")
  
  # data contains information for each angling sequence
  patch_leaving_data = fread(file = "utils/data/raw_data/patch_leaving_data.csv")
  
  # data contains video labels and gps for those with working recording
  video_data = fread(file = "utils/data/raw_data/video_data.csv")
  
  # data only contains gps coordinates for all participants (also for those with missing recordings)
  gps_data = fread(file = "utils/data/raw_data/gps_data.csv")
  
  ################################################################################
  # descriptive analyses
  ################################################################################
  # figure with participation
  participation = catch_data %>% 
    group_by(participant_id) %>% 
    mutate(n = n()) %>% 
    ggplot(aes(x = interaction(day, year), y = fct_reorder(as.factor(participant_id), -n))) + 
    geom_tile(color = "black", fill ="cadetblue", alpha = .4 )+ 
    xlab("Competition (d.yyyy)") + 
    ylab("Participant ID") +
    #theme_classic()+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.5, .15),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          axis.ticks.y = element_blank(),
          legend.background = element_rect(colour = "white", fill = "white"))
  participation
  # save figure in different file formats
  ggsave(participation, file = "1_icefishing_social_foraging/output/participation_plot.png", 
         width = 9, height = 5,bg = "white", dpi = 600, units = "cm")
  ggsave(participation, file = "1_icefishing_social_foraging/output/participation_plot.svg", 
         width = 9, height = 5,bg = "white", dpi = 600, units = "cm",device = svglite::svglite)
  ggsave(participation, file = "1_icefishing_social_foraging/output/participation_plot.pdf", 
         width = 9, height = 5,bg = "white", dpi = 600, units = "cm")
  ggsave(participation, file = "1_icefishing_social_foraging/output/participation_plot.tiff", 
         width = 9, height = 5,bg = "white", dpi = 600, units = "cm")
  
  # how many participants per competition?
  catch_data %>% 
    group_by(day, year) %>% 
    summarize(n = n())
  
  # avg age
  mean(survey_data$age, na.rm= T)
  sd(survey_data$age, na.rm= T)
  sum(is.na(survey_data$age))
  
  # gender distr.
  mean(as.numeric(as.factor(catch_data$sex))-1)
  
  # make table with detailed results for each day
  catch_data %>% 
    group_by(year, day) %>% 
    summarize(m_catch = mean(catch, na.rm = T),
              med_catch = median(catch, na.rm = T),
              sd_catch = sd(catch, na.rm = T))
  p_data = left_join(catch_data, survey_data)
  
  p_data %>% 
    group_by(day, year, participant_id) %>%
    slice(1) %>%
    group_by(year, day) %>% 
    mutate(fished_before = ifelse(fished_before == "y", 1, 0)) %>% 
    summarize(mean_sex = 1-mean(as.numeric(as.factor(sex))-1, na.rm = T),
              sex_missing = sum(is.na(sex)),
              mean_age = mean(age, na.rm = T),
              median_age = median(age, na.rm =T ),
              sd_age = sd(age, na.rm = T),
              age_missing = sum(is.na(age)),
              mean_skill = mean(angling_skill, na.rm = T),
              median_skill = median(angling_skill, na.rm = T),
              sd_skill = sd(angling_skill, na.rm = T),
              skill_missing = sum(is.na(angling_skill)),
              mean_experience = mean(years_icefishing, na.rm = T),
              median_experience = median(years_icefishing, na.rm = T),
              sd_experience = sd(years_icefishing, na.rm = T),
              experience_missing = sum(is.na(years_icefishing)),
              fished_before = mean(fished_before, na.rm = T),
              missing_fished_before = sum(is.na(fished_before)))
  
  # how much time did individuals spend angling vs relocation and other states
  rel_angling_time = spot_selection_data %>% 
    filter(!exclude) %>% 
    group_by(day, year, camera_id) %>% 
    summarize(rel_angling_time = sum(angling_time) / max(time_end))
  mean(rel_angling_time$rel_angling_time)
  sd(rel_angling_time$rel_angling_time)
  
  # how many spots visited during a competition?
  n_spots_visited = spot_selection_data %>% 
    filter(!exclude) %>% 
    group_by(day, year, camera_id) %>% 
    summarize(n_spots = max(spot_id))
  mean(n_spots_visited$n_spots)
  sd(n_spots_visited$n_spots)
  
  # total number of foraging locations analysed
  nrow(spot_selection_data %>% filter(!exclude))
  
  # network analysis with dbscan
  radii = seq(from = 0, to = 500, by = 10)
  n_iter = 50
  source("1_icefishing_social_foraging/scripts/create_network_data.r")
  if (file.exists(file = "utils/data/processed_data/network_data.rds")){
    
    iterated_results = readRDS(file = "utils/data/processed_data/network_data.rds")
    
  } else {
    create_network_data(gps_data, radii, n_iter)
    iterated_results = readRDS(file = "utils/data/processed_data/network_data.rds")
    
  }
  
  # plot results
  changepoint = iterated_results %>% 
    group_by(r) %>% 
    summarize(random = mean(random),
              data = mean(data)) %>% 
    mutate(delta = abs(random-data)) %>% 
    filter(r > 100) %>% 
    filter(delta == min(delta))
  changepoint = changepoint$r
  colors <- c("Random" = "cadetblue", "Data" = "black")
  networks = ggplot() +
    stat_summary(data = iterated_results %>% dplyr::select(data, iteration, r) %>% filter(iteration == 1), aes(x = r, y = data, color = "Data"), alpha = .3, size = .4) + 
    stat_lineribbon(data = iterated_results %>% 
                      dplyr::select(random, iteration, r, timestamp) %>% 
                      group_by(r, iteration) %>% summarize(random = mean(random)), 
                    aes(x = r, y = random, color = "Random"),.width = .95, alpha = .4, fill = "cadetblue2")+
    scale_color_manual(values = colors)+
    xlab(expression(paste("Distance threshold ", epsilon~(m), sep = "")))+
    ylab("Relative community size")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          plot.title = element_text(size = 7),
          plot.tag = element_text(size = 7, face = "bold"),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.2, .8),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = .2)) + 
    geom_vline(xintercept = changepoint, linetype = "dotted")+
    annotate("segment", x = 80, y = .25, xend = 100, yend = .15,
             arrow = arrow(type = "closed", length = unit(.15, "cm"))) + 
    annotate("segment", x = 430, y = .5, xend = 410, yend = .65,
             arrow = arrow(type = "closed", length = unit(.15, "cm"))) + 
    annotate("text", x = 430, y = .45, label = "Globally Dispersed", size = 7/.pt) + 
    annotate("text", x = 80, y = .3, label = "Locally Dense", size = 7/.pt)
  networks
  # save figure in different file formats
  ggsave(networks, file = "1_icefishing_social_foraging/output/networks.png", 
         width = 9, height = 9,bg = "white", dpi = 600, units = "cm")
  ggsave(networks, file = "1_icefishing_social_foraging/output/networks.svg", 
         width = 9, height = 9,bg = "white", dpi = 600, units = "cm",device = svglite::svglite)
  ggsave(networks, file = "1_icefishing_social_foraging/output/networks.pdf", 
         width = 9, height = 9,bg = "white", dpi = 600, units = "cm")
  ggsave(networks, file = "1_icefishing_social_foraging/output/networks.tiff", 
         width = 9, height = 9,bg = "white", dpi = 600, units = "cm")
  
  # how much time spent at angling spots
  spot_selection_data %>% 
    filter(!exclude) %>% 
    summarize(mean(angling_time),sd(angling_time))
  as.data.frame(spot_selection_data %>% 
                  filter(!exclude) %>% 
                  group_by(Return) %>% 
                  summarize(mean(angling_time),sd(angling_time)))
  
  
  # variation in catch success across lakes
  success = catch_data %>% 
    group_by(day, year) %>% 
    summarize(c = mean(catch), sd(catch)) 
  
  # percentage successful spots
  perc_success = spot_selection_data %>% 
    filter(!exclude) %>%
    group_by(day, year) %>% 
    summarize(r = mean(Return))
  
  # how many successful spots with time series information
  patch_leaving_data %>% 
    filter(!exclude) %>% 
    group_by(camera_id, day, year, sequence) %>% 
    summarize(n = sum(reward)) %>%
    filter(n > 0)
  
  # data for both models
  spot_selection_data %>% 
    group_by(year, day) %>% 
    summarize(n = n(),
              missing = sum(exclude),
              percent_missing = sum(exclude)/n())
  
  spot_selection_data %>% 
    group_by(year, day, camera_id) %>% 
    summarize(n = n(),
              missing = sum(exclude),
              percent_missing = sum(exclude)/n(),
              n_missing = sum(exclude),
              n_successful = mean(Return)) %>% 
    group_by(year, day) %>% 
    summarize(mean_n = mean(n),
              median_n = median(n),
              sd_n = sd(n),
              median_missing = median(percent_missing),
              mean_missing = mean(percent_missing),
              n_missing_mean = mean(n_missing),
              sd_n_missing = sd(n_missing),
              sd_perc_missing = sd(percent_missing),
              n_missing_median = median(n_missing),
              percent_complete = mean(percent_missing < .05),
              percent_successful = mean(n_successful),
              median_successful = median(n_successful),
              sd_successful = sd(n_successful))
  
  spot_selection_data %>% 
    filter(!exclude) %>% 
    group_by(year, day, camera_id) %>% 
    summarize(n = n(),
              n_successful = mean(Return)) %>% 
    group_by(year, day) %>% 
    summarize(mean_n = mean(n),
              median_n = median(n),
              sd_n = sd(n),
              percent_successful = mean(n_successful),
              median_successful = median(n_successful),
              sd_successful = sd(n_successful))
  
  spot_selection_data %>% 
    group_by(year, day) %>% 
    mutate(n_id = length(unique(camera_id))) %>% 
    group_by(year, day) %>% 
    summarize(n = n() + mean(n_id),
              n_id = mean(n_id),
              percent_missing = sum(exclude) / (n() + n_id),
              n_missing = sum(exclude),
              n_analysed = n - n_id - n_missing)
  
  patch_leaving_data %>% 
    group_by(year, day, camera_id, sequence) %>% 
    summarize(angling_time = n()*10,
              total_catch = sum(reward),
              missing = mean(as.numeric((exclude)>0))) %>% 
    group_by(year, day) %>% 
    summarize(n_missing = sum(missing),
              percent_missing = mean(missing),
              n_spots = n())
  
  patch_leaving_data %>% 
    group_by(year, day, camera_id, sequence) %>% 
    slice(1) %>% 
    group_by(year, day) %>% 
    mutate(n_id = length(unique(camera_id))) %>% 
    group_by(year, day) %>% 
    summarize(n = n() + mean(n_id),
              n_id = mean(n_id),
              percent_missing = sum(exclude) / (n() + n_id),
              n_missing = sum(exclude),
              n_analysed = n - n_id - n_missing)
  
  
  patch_leaving_data %>% 
    group_by(year, day, camera_id, sequence) %>% 
    summarize(angling_time = n()*10,
              total_catch = sum(reward),
              missing = mean(as.numeric((exclude)>0))) %>% 
    group_by(year, day, camera_id) %>% 
    summarize(n_missing = sum(missing),
              percent_missing = mean(missing),
              n_spots = n()) %>% 
    group_by(year, day) %>% 
    summarize(mean_locations_trip = mean(n_spots),
              sd_locations_trip = sd(n_spots),
              median_locations_trip = median(n_spots),
              mean_missing_trip = mean(n_missing),
              sd_missing_trip = sd(n_missing),
              median_missing_trip = median(n_missing),
              mean_perc_missing_trip = mean(percent_missing),
              sd_perc_missing_trip = sd(percent_missing),
              median_perc_missing_trip = median(percent_missing),
              trips_with_less_than_5 = mean(percent_missing<.05))
  
  patch_leaving_data %>% 
    group_by(year, day, camera_id, sequence) %>% 
    summarize(angling_time = n()*10,
              total_catch = sum(reward),
              missing = mean(as.numeric((exclude)>0))) %>%
    filter(!missing) %>% 
    group_by(year, day, camera_id) %>% 
    summarize(n_spots = n(),
              percent_successful = mean(total_catch >0)) %>%
    group_by(year, day) %>% 
    summarize(perc_successful = mean(percent_successful),
              sd_successful = sd(percent_successful),
              median_successful = median(percent_successful))
  
  patch_leaving_data %>% 
    group_by(year, day, camera_id, sequence) %>% 
    summarize(angling_time = n()*10,
              total_catch = sum(reward),
              missing = mean(as.numeric((exclude)>0))) %>%
    filter(!missing) %>% 
    group_by(year, day, camera_id, x=total_catch > 0) %>% 
    summarize(mean_angling_time = mean(angling_time),
              mean_catch = mean(total_catch)) %>%
    group_by(year, day, x) %>% 
    filter(x) %>% 
    summarize(mean_angling_time_trip = mean(mean_angling_time),
              sd_angling_time_trip = sd(mean_angling_time),
              median_angling_time_trip = median(mean_angling_time),
              mean_catch_trip = mean(mean_catch),
              sd_catch_trip = sd(mean_catch),
              median_catch_trip = median(mean_catch))
  
  patch_leaving_data %>% 
    group_by(year, day, camera_id, sequence) %>% 
    summarize(angling_time = n()*10,
              total_catch = sum(reward),
              missing = mean(as.numeric((exclude)>0))) %>%
    filter(!missing) %>% 
    group_by(year, day, camera_id, x=total_catch > 0) %>% 
    summarize(mean_angling_time = mean(angling_time),
              mean_catch = mean(total_catch)) %>%
    group_by(year, day, x) %>% 
    filter(!x) %>% 
    summarize(mean_angling_time_trip = mean(mean_angling_time),
              sd_angling_time_trip = sd(mean_angling_time),
              median_angling_time_trip = median(mean_angling_time),
              mean_catch_trip = mean(mean_catch),
              sd_catch_trip = sd(mean_catch),
              median_catch_trip = median(mean_catch))
  
  # Supplementary figures
  catch_day = catch_data %>% 
    group_by(day, year, camera_id) %>% 
    mutate(x=mean(catch), group= interaction(day, year)) %>% 
    ggplot(aes(y=x/1000, x = as.factor(as.numeric(interaction(day, year))))) + 
    geom_boxplot(outlier.shape = NA) + 
    ylab("Catch (kg)")+
    labs(tag = "A") +
    xlab("Lake")+
    geom_jitter(alpha = .4, size = .6) + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  
  ratio = spot_selection_data %>% 
    group_by(day, year) %>% 
    summarize(mean(Return))
  min(ratio$`mean(Return)`)
  max(ratio$`mean(Return)`)
  median(ratio$`mean(Return)`)
  
  ratios = spot_selection_data %>% 
    group_by(day, year, camera_id) %>% 
    summarize(r = mean(Return)) %>% 
    ggplot(aes(y=r*100, x = as.factor(as.numeric(interaction(day, year))))) + 
    geom_boxplot(outlier.shape = NA) + 
    ylab("% Successful Spots")+
    xlab("Lake")+
    labs(tag = "B") +
    geom_jitter(alpha = .4,size = .6) + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  
  spot_selection_data$camera_id = as.character(spot_selection_data$camera_id)
  spot_selection_data = left_join(spot_selection_data, catch_data %>% dplyr::select(camera_id, day, year, catch))
  test = spot_selection_data %>% 
    group_by(day, year) %>% 
    summarize(r = mean(Return)*100,
              r_sd = sd(Return)*100,
              c = mean(catch)/1000,
              c_sd = sd(catch/1000),
              n = n())
  
  test = test %>% 
    ungroup() %>% 
    mutate(scale_r = (r - mean(r))/sd(r),
           scale_c = (c - mean(c))/sd(c))
  
  correl = brms::brm(data = test, formula = scale(r)~1+scale(c), 
                     chains = 4,
                     cores = 4)
  
  correl_day = test %>% 
    ggplot(aes(x = r, y = c)) + 
    geom_point(size = .6) + 
    geom_errorbar(aes(ymin = c - 2*(c_sd/sqrt(n)), ymax = c + 2*(c_sd/sqrt(n)))) + 
    geom_errorbarh(aes(xmin = r - 2*(r_sd/sqrt(n)), xmax = r + 2*(r_sd/sqrt(n)))) + 
    labs(tag = "C") +
    coord_cartesian(xlim = c(0, 35), ylim = c(0, 2.3))+
    scale_x_continuous(breaks = seq(from = 0, to = 35, by = 5))+
    xlab("Avg. % Successful Spots") + 
    ylab("Avg. Catch (kg)") +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.background = element_rect(colour = "black", fill = "white"))
  
  # catch in n and kg correlated
  test = spot_selection_data
  
  correl_time = brm(data = test, formula = scale(angling_time)~1+scale(total_catch), 
                    chains = 4,
                    cores = 4)
  correl_time2 = brm(data = test, formula = (angling_time)~1+(total_catch), 
                     chains = 4,
                     cores = 4)
  
  get_variables(correl_time)
  prediction_data = correl_time2 %>% 
    spread_draws(b_Intercept, 
                 b_total_catch) %>% 
    slice(rep(1:n(), each = 151)) %>% 
    group_by(.draw) %>%
    mutate(catch_pred = seq(from = 0, to = 150, by = 1)) %>% 
    ungroup() %>% 
    mutate(posterior_prediction = b_Intercept + b_total_catch * catch_pred)
  
  cor.test(test$total_catch, test$angling_time)
  catch_time = 
    ggplot() + 
    stat_lineribbon(data = prediction_data, aes(x = catch_pred, y = posterior_prediction/60),.width = c(.66, .95), alpha = .4, size = .6)+
    geom_point(data = spot_selection_data, aes(x = total_catch, y = angling_time / 60), alpha = .2, size = .6)+
    scale_fill_manual(values = c("darkgrey","lightgrey"), guide = "none")+
    labs(tag = "D") +
    xlab("Catch (n fish)") + 
    ylab("Angling Time (min)") + 
    #geom_smooth(method = "lm", color = "black") + 
    #stat_summary(color = "cyan3")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.background = element_rect(colour = "black", fill = "white"))# +
  #  annotate("text", x = 125, y = 25, label = "r = .715 [.670, .733]", size = 6)
  catch_time
  
  # catch in kg and n fish per day
  n_dat = spot_selection_data %>% 
    group_by(day, year, camera_id) %>% 
    summarize(n = sum(total_catch, na.rm = T))
  n_dat$camera_id = as.character(n_dat$camera_id)
  n_dat = left_join(n_dat, catch_data)
  
  
  correl_n = brm(data = n_dat, formula = scale(catch)~1+scale(n), 
                 chains = 4,
                 cores = 4)
  correl_n2 = brm(data = n_dat, formula = (catch)~1+(n), 
                  chains = 4,
                  cores = 4)
  
  prediction_data_n = correl_n2 %>% 
    spread_draws(b_Intercept, 
                 b_n) %>% 
    slice(rep(1:n(), each = 275)) %>% 
    group_by(.draw) %>%
    mutate(n_pred = seq(from = 0, to = 274, by = 1)) %>% 
    ungroup() %>% 
    mutate(posterior_prediction = b_Intercept + b_n * n_pred)
  
  
  catch_n_kg =  
    ggplot() + 
    stat_lineribbon(data = prediction_data_n, aes(x = n_pred, y = posterior_prediction/1000),.width = c(.66, .95), alpha = .4, size = .6)+
    geom_point(data = n_dat, aes(x = n, y = catch/1000), alpha = .2, size = .6)+
    scale_fill_manual(values = c("darkgrey","lightgrey"), guide = "none")+
    labs(tag = "E") +
    xlab("Catch (n fish)") + 
    ylab("Catch (kg)") + 
    scale_x_continuous(breaks = seq(from = 0, to = 250, by = 50))+
    scale_y_continuous(breaks = seq(from = 0, to = 9, by = 1))+
    #geom_smooth(method = "lm", color = "black") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          legend.title = element_blank(),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.background = element_rect(colour = "black", fill = "white"))
  catch_n_kg
  
  # make panel for supplement
  panel_lake_characteristics = arrangeGrob(
    grobs = list(catch_day, 
                 ratios,
                 correl_day,
                 catch_time,
                 catch_n_kg),
    widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    heights = c(1, 1, 2),
    layout_matrix = rbind(c(1, 1, 1, 1, 1, 1, 3, 3, 3, 3),
                          c(2, 2, 2, 2, 2, 2, 3, 3, 3, 3),
                          c(4, 4, 4, 4, 4, 5, 5, 5, 5, 5))
  )
  # save figure in different file formats
  ggsave(panel_lake_characteristics, file = "1_icefishing_social_foraging/output/panel_lake_characteristics.png", 
         width = 18.4, height = 14,bg = "white", dpi = 600, units = "cm")
  ggsave(panel_lake_characteristics, file = "1_icefishing_social_foraging/output/panel_lake_characteristics.svg", 
         width = 18.4, height = 14,bg = "white", dpi = 600, units = "cm",device = svglite::svglite)
  ggsave(panel_lake_characteristics, file = "1_icefishing_social_foraging/output/panel_lake_characteristics.pdf", 
         width = 18.4, height = 14,bg = "white", dpi = 600, units = "cm")
  ggsave(panel_lake_characteristics, file = "1_icefishing_social_foraging/output/panel_lake_characteristics.tiff", 
         width = 18.4, height = 14,bg = "white", dpi = 600, units = "cm")
  
  # save supplementary figures
  time_allocation = video_data %>% 
    group_by(day, year, camera_id) %>% 
    summarize(angling = sum(angling, na.rm = T),
              relocating = sum(relocating, na.rm = T),
              drilling = sum(drilling, na.rm = T))
  
  write_csv(time_allocation, file = "utils/data/raw_data/time_allocation_data.csv")
  time_allocation = read.csv(file = "utils/data/raw_data/time_allocation_data.csv")
  n_spots_id = spot_selection_data %>% group_by(camera_id, day, year) %>% summarize(n_spots = n())
  time_allocation = left_join(time_allocation, n_spots_id)
  time_allocation_plot = time_allocation %>%
    group_by(day,year) %>% 
    arrange(angling / (angling + relocating + drilling)) %>% 
    mutate(id = 1:n()) %>% 
    pivot_longer(angling:drilling) %>% 
    ggplot(aes(fill=name, y=value, x=id)) + 
    geom_bar(position="fill", stat="identity", alpha = .7) +#, color = "lightgrey", alpha = .2 
    scale_fill_brewer(palette = "Set2")+
    facet_wrap(~day*year, nrow = 2) + 
    labs(fill = "Activity")+ 
    ylab("% Time Allocation")+
    xlab("Participant ID")+
    scale_y_continuous(labels = scales::percent)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 5),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "bottom",
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
          legend.box.margin = margin(1, 1, 1, 1),
          legend.key.size = unit(.1, "cm"))
  time_allocation_plot
  
  # save figure in different file formats
  ggsave(time_allocation_plot, file = "1_icefishing_social_foraging/output/time_allocation_plot.png", 
         width = 18, height = 9,bg = "white", dpi = 600, units = "cm")
  ggsave(time_allocation_plot, file = "1_icefishing_social_foraging/output/time_allocation_plot.svg", 
         width = 18, height = 9,bg = "white", dpi = 600, units = "cm",device = svglite::svglite)
  ggsave(time_allocation_plot, file = "1_icefishing_social_foraging/output/time_allocation_plot.pdf", 
         width = 18, height = 9,bg = "white", dpi = 600, units = "cm")
  ggsave(time_allocation_plot, file = "1_icefishing_social_foraging/output/time_allocation_plot.tiff", 
         width = 18, height = 9,bg = "white", dpi = 600, units = "cm")
  
  # how much time relocating per spot?
  time_allocation %>% 
    summarize(mean(relocating / n_spots, na.rm = T),
              sd(relocating/n_spots, na.rm = T))
  
  avg_relocation_time = 
    time_allocation %>%  
    group_by(camera_id, day, year) %>% 
    summarize(time = mean(relocating / n_spots, na.rm = T)) %>% 
    ggplot(aes(x = time)) + 
    geom_histogram(color = "black")+
    ylab("ID count")+
    xlab("Avg. Relocation Time (s)")+
    facet_wrap(~day*year, nrow = 2) +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  
  # save figure in different file formats
  ggsave(avg_relocation_time, file = "1_icefishing_social_foraging/output/avg_relocation_time.png", 
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm")
  ggsave(avg_relocation_time, file = "1_icefishing_social_foraging/output/avg_relocation_time.svg", 
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm",device = svglite::svglite)
  ggsave(avg_relocation_time, file = "1_icefishing_social_foraging/output/avg_relocation_time.pdf", 
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm")
  ggsave(avg_relocation_time, file = "1_icefishing_social_foraging/output/avg_relocation_time.tiff",
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm")
  
  
  ################################################################################
  # REVISION: Distribution of patches within and across lakes
  # estimate overdispersion of patches with negative binomial model
  neg_bin_data = fread("utils/data/raw_data/spot_selection_data.csv")
  neg_bin_data = neg_bin_data %>% filter(!exclude)
  
  # fit negative binomial model with zero inflation
  if (!file.exists("utils/data/processed_data/dispersion_fit.rds")){
    ngzi = brms::brm(data = neg_bin_data, 
                     formula = bf(total_catch ~ (1|lake_id),
                                  shape ~ (1|lake_id),
                                  zi ~ (1|lake_id)),
                     family = zero_inflated_negbinomial(),
                     chains = 4,
                     cores = 4)
    saveRDS(ngzi, file = "utils/data/processed_data/dispersion_fit.rds")
  }
  
  ngzi = readRDS(file = "utils/data/processed_data/dispersion_fit.rds")
  
  avg = ngzi %>% 
    spread_draws(b_Intercept, b_shape_Intercept) %>% 
    mutate(additional_variance = exp(b_Intercept)^2 / exp(b_shape_Intercept)) %>% 
    ggplot(aes(x = additional_variance)) + 
    stat_pointinterval()
  lake_effects = ngzi %>% 
    spread_draws(b_Intercept, b_shape_Intercept, r_lake_id[lake, parameter],
                 r_lake_id__shape[lake, parameter]) %>% 
    mutate(additional_variance = 1 / exp(b_shape_Intercept + r_lake_id__shape)) %>% 
    ggplot(aes(x = lake, y = additional_variance)) + 
    stat_pointinterval()
  lake_effects
  
  # visualize
  ngzi_plot_data=ngzi %>% 
    spread_draws(b_Intercept, b_shape_Intercept, r_lake_id[lake, parameter],
                 r_lake_id__shape[lake, parameter]) %>% 
    mutate(additional_variance = 1 / exp(b_shape_Intercept + r_lake_id__shape))
  main_effects = ngzi %>% 
    spread_draws(b_Intercept, b_shape_Intercept) %>% 
    mutate(additional_variance = 1 / exp(b_shape_Intercept))
  main_effects %>% 
    summarize(mean = mean(additional_variance),
              lower = quantile(additional_variance, .025),
              upper = quantile(additional_variance, .975))
  
  overdispersion_plot = ggplot(data = ngzi_plot_data, aes(x = lake, y = additional_variance)) + 
    annotate("rect", ymin = quantile(main_effects$additional_variance, .025), ymax =  quantile(main_effects$additional_variance, .975), xmin = .75, xmax = 10.25,
             alpha = .5,fill = "lightgrey")+
    annotate("segment", x = .75, y = quantile(main_effects$additional_variance, .5), xend = 10.25, yend = quantile(main_effects$additional_variance, .5)) + 
    stat_pointinterval() + 
    xlab("Lake") +
    labs(tag = "B")+
    ylab(expression("Overdispersion"~phi^{-1})) +
    scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  
  # visualize
  ngzi_plot_data_mean=ngzi %>% 
    spread_draws(b_Intercept, b_shape_Intercept, b_zi_Intercept, r_lake_id[lake, parameter],
                 r_lake_id__shape[lake, parameter],
                 r_lake_id__zi[lake, parameter]) %>% 
    mutate(mean = exp(b_Intercept + r_lake_id),
           p_zero = inv_logit(b_zi_Intercept + r_lake_id__zi))
  main_effects_mu = ngzi %>% 
    spread_draws(b_Intercept, b_shape_Intercept, b_zi_Intercept) %>% 
    mutate(mean = exp(b_Intercept),
           p_zero = inv_logit(b_zi_Intercept))
  
  overdispersion_plot_mean = ggplot(data = ngzi_plot_data_mean, aes(x = lake, y = mean)) + 
    annotate("rect", ymin = quantile(main_effects_mu$mean, .025), ymax =  quantile(main_effects_mu$mean, .975), xmin = .75, xmax = 10.25,
             alpha = .5,fill = "lightgrey")+
    annotate("segment", x = .75, y = quantile(main_effects_mu$mean, .5), xend = 10.25, yend = quantile(main_effects_mu$mean, .5)) + 
    stat_pointinterval() + 
    xlab("Lake") +
    labs(tag = "A")+
    ylab(expression("Expectation"~mu)) +
    scale_x_continuous(breaks = seq(from = 1, to = 10, by = 1))+
    scale_y_continuous(breaks = seq(from = 0, to = 10, by = 1))+
    #geom_hline(yintercept = 0, linetype = "dashed") +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  overdispersion_plot_mean
  
  # make panel for supplement
  panel_dispersion = arrangeGrob(
    grobs = list(overdispersion_plot_mean,
                 overdispersion_plot),
    layout_matrix = rbind(c(1),
                          c(2))
  )
  
  
  ggsave(panel_dispersion, file = "1_icefishing_social_foraging/output/overdispersion_plot.png", 
         width = 12, height = 10,bg = "white", dpi = 600, units = "cm")
  ggsave(panel_dispersion, file = "1_icefishing_social_foraging/output/overdispersion_plot.svg", 
         width = 12, height = 10,bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  ggsave(panel_dispersion, file = "1_icefishing_social_foraging/output/overdispersion_plot.pdf", 
         width = 12, height = 10,bg = "white", dpi = 600, units = "cm")
  ggsave(panel_dispersion, file = "1_icefishing_social_foraging/output/overdispersion_plot.tiff",
         width = 12, height = 10,bg = "white", dpi = 600, units = "cm")
  
  # mean
  neg_bin_data %>% 
    group_by(lake_id) %>% 
    summarize(m = mean(total_catch),
              sd = sd(total_catch))
  
  ################################################################################
  # how much time drilling holes per spot?
  n_holes_drilled = video_data %>% 
    group_by(day, year, camera_id) %>% 
    mutate(drilling_seq = rleid(drilling)) %>% 
    filter(drilling == 1) %>% 
    group_by(day, year, camera_id) %>% 
    summarize(n_holes_drilled = length(unique(drilling_seq)))
  
  time_allocation = left_join(time_allocation, n_holes_drilled)
  time_allocation %>% 
    summarize(mean(drilling / n_holes_drilled, na.rm = T),
              sd(drilling/n_holes_drilled, na.rm = T))
  
  avg_drilling_time = 
    time_allocation %>% 
    group_by(camera_id, day, year) %>% 
    summarize(time = mean(drilling / n_holes_drilled, na.rm = T)) %>% 
    ggplot(aes(x = time)) + 
    geom_histogram(color = "black")+
    ylab("ID count")+
    xlab("Avg. Drilling Time (s)")+
    scale_x_continuous(breaks = seq(from = 0, to = 200, by = 25)) + 
    facet_wrap(~day*year, nrow = 2) +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  # save figure in different file formats
  ggsave(avg_drilling_time, file = "1_icefishing_social_foraging/output/avg_drilling_time.png", 
         width = 18, height = 8,bg = "white", dpi = 300, units = "cm")
  ggsave(avg_drilling_time, file = "1_icefishing_social_foraging/output/avg_drilling_time.svg", 
         width = 18, height = 8,bg = "white", dpi = 300, units = "cm", device = svglite::svglite)
  ggsave(avg_drilling_time, file = "1_icefishing_social_foraging/output/avg_drilling_time.pdf", 
         width = 18, height = 8,bg = "white", dpi = 300, units = "cm")
  ggsave(avg_drilling_time, file = "1_icefishing_social_foraging/output/avg_drilling_time.tiff",
         width = 18, height = 8,bg = "white", dpi = 300, units = "cm")
  
  # avg relocation distance
  avg_relocation = 
    spot_selection_data %>% filter(step>0) %>% 
    group_by(camera_id, day, year) %>% 
    summarize(step = mean(step)) %>% 
    ggplot(aes(x = step)) + 
    geom_histogram(color = "black")+
    ylab("ID count")+
    xlab("Avg. Relocation Distance (m)")+
    facet_wrap(~day*year, nrow = 2) +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  
  # save figure in different file formats
  ggsave(avg_relocation, file = "1_icefishing_social_foraging/output/avg_relocation.png", 
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm")
  ggsave(avg_relocation, file = "1_icefishing_social_foraging/output/avg_relocation.svg", 
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  ggsave(avg_relocation, file = "1_icefishing_social_foraging/output/avg_relocation.pdf", 
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm")
  ggsave(avg_relocation, file = "1_icefishing_social_foraging/output/avg_relocation.tiff",
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm")
  
  spot_selection_data %>% filter(step>0) %>% 
    group_by(camera_id, day, year) %>% 
    summarize(step = sum(step)) %>% 
    ungroup() %>% 
    summarize(m_step = mean(step), 
              sd = sd(step, na.rm = T))
  
  total_relocation = 
    spot_selection_data %>% filter(step>0) %>% 
    group_by(camera_id, day, year) %>% 
    summarize(step = sum(step)) %>% 
    ggplot(aes(x = step/1000)) + 
    geom_histogram(color = "black")+
    ylab("ID count")+
    xlab("Total Relocation Distance (km)")+
    facet_wrap(~day*year, nrow = 2) +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  
  # save figure in different file formats
  ggsave(total_relocation, file = "1_icefishing_social_foraging/output/total_relocation.png", 
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm")
  ggsave(total_relocation, file = "1_icefishing_social_foraging/output/total_relocation.svg", 
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  ggsave(total_relocation, file = "1_icefishing_social_foraging/output/total_relocation.pdf", 
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm")
  ggsave(total_relocation, file = "1_icefishing_social_foraging/output/total_relocation.tiff",
         width = 18, height = 8,bg = "white", dpi = 600, units = "cm")
  
  # average angling time successful/unsuccessful spots
  exact_angling_time = video_data %>% 
    group_by(day, year, camera_id) %>% 
    mutate(angling_seq = rleid(angling)) %>% 
    filter(angling == 1) %>%
    group_by(day, year, camera_id, angling_seq) %>% 
    summarize(Return = ifelse(sum(fish_catch)>0, 1, 0),
              time = n()) %>% 
    group_by(day, year, camera_id, Return) %>% 
    summarize(time = mean(time))
  
  avg_angling_time = 
    exact_angling_time %>% 
    filter(time < 1800) %>%
    mutate(Return = ifelse(Return == 0, "Unsuccessful Spots", "Successful Spots")) %>% 
    ggplot(aes(x = time/60, fill = as.factor(Return), group = as.factor(Return))) + 
    geom_histogram(color = "black", position = position_dodge(.5), alpha = .6, bins = 20)+
    ylab("ID count")+
    scale_fill_brewer(palette = "Set2")+
    xlab("Avg. Angling Time (min)")+
    labs(fill = "") +
    facet_wrap(~day*year, nrow = 2) +
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "bottom",
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  avg_angling_time
  # n = 27
  exact_angling_time %>% 
    filter(time >= 1800)
  # save figure in different file formats
  ggsave(avg_angling_time, file = "1_icefishing_social_foraging/output/avg_angling_time.png", 
         width = 18, height = 10,bg = "white", dpi = 600, units = "cm")
  ggsave(avg_angling_time, file = "1_icefishing_social_foraging/output/avg_angling_time.svg", 
         width = 18, height = 10,bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  ggsave(avg_angling_time, file = "1_icefishing_social_foraging/output/avg_angling_time.pdf", 
         width = 18, height = 10,bg = "white", dpi = 600, units = "cm")
  ggsave(avg_angling_time, file = "1_icefishing_social_foraging/output/avg_angling_time.tiff",
         width = 18, height = 10,bg = "white", dpi = 600, units = "cm")
  
  
  ################################################################################
  # create figure 1
  ################################################################################
  # data containing depth contours
  lake_data = readRDS("utils/data/raw_data/lake_data/depth_profiles.rds")
  
  # data containing boundary of competition area
  competition_area = readRDS("utils/data/raw_data/lake_data/study_lakes.rds")
  
  # get in correct order
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
  
  ################################################################################
  # Create Plot with Trajectories on Depth Profile
  ################################################################################
  day5 = gps_data %>% filter(day == 5 & year == 2022)
  saveRDS(day5, file = "utils/data/raw_data/gps_data_day5.rds")
  
  LongLatToUTM<-function(x,y,zone){
    xy <- data.frame(camera_id = 1:length(x), X = x, Y = y)
    coordinates(xy) <- c("X", "Y")
    proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example# +datum=WGS84
    res <- spTransform(xy, CRS(paste("+proj=utm +datum=WGS84 +zone=",zone," ellps=WGS84",sep='')))
    return(as.data.frame(res))
  }
  
  # project coordinates to XY
  xy = LongLatToUTM(day5$position_long, day5$position_lat, zone = 35)
  day5$X = xy[, 2]
  day5$Y = xy[, 3] 
  day5 = st_as_sf(day5, coords = c("X", "Y"))
  st_crs(day5) <- "epsg:3067"
  day5$X = xy[, 2]
  day5$Y = xy[, 3]
  
  # project depth profile in common crs
  a = competition_area[[9]]
  d = st_transform(depth_profiles[[9]], crs = "epsg:3067")
  
  # example traj
  subset = spot_selection_data %>% 
    filter(day == 5 & year == 2022) %>% 
    filter(camera_id == 5)
  xy = subset[, c("x", "y")]
  subset = st_as_sf(subset, coords = c("x", "y"))
  st_crs(subset) <- "epsg:3067"
  subset$x = xy$x
  subset$y = xy$y
  rectangle = raster::extent(subset)
  
  library(ggspatial)
  trajectories = ggplot() + 
    geom_sf(data = depth_profiles[[9]], color ="black", alpha = .3) +
    geom_path(data = day5, aes(x = X, y = Y, group = camera_id), color = "darkgreen", alpha = .4) + 
    geom_path(data = day5 %>% filter(camera_id == 5), aes(x = X, y = Y, group = camera_id), color = "black", alpha = .8, size = .6) + 
    geom_rect(aes(xmin = rectangle[1]-15, xmax = rectangle[2]+15, ymin = rectangle[3]-15, ymax = rectangle[4]+15), fill = NA, color = "black", linetype = "dashed", size = .6)+
    ylab("Northing")+
    xlab("Easting")+
    labs(tag = "D")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white")) + 
    annotation_scale()
  
  # zoom into map
  profile_cropped = st_crop(depth_profiles[[9]], xmin = min(subset$x) - 15, xmax = max(subset$x) + 15,
                            ymin = min(subset$y) - 15, ymax = max(subset$y) + 15)
  zoom = raster::extent(profile_cropped)
  
  # assign color values
  colors <- c("Successful Spot" = "cyan3", "Unsuccessful Spot  " = "darkgreen")
  
  # for fish icons
  #devtools::install_github("mattiaghilardi/fishualize", ref = "fix-fishapes")
  #library(fishualize)
  
  # read images to display on plot
  relocation <- readPNG("utils/data/raw_data/plot_images/relocating.png")
  angling <- readPNG("utils/data/raw_data/plot_images/angling.png")
  
  # coordinates for arrows in plot
  coords_arrows = subset %>% filter(Return == 0) 
  coords_arrows = coords_arrows[30, c("x", "y")] 
  
  # depth
  # sample points from lakebed
  pts = list()
  for (i in 1:nrow(profile_cropped)){
    pts[[i]] = st_line_sample(profile_cropped[i,] %>% st_cast("LINESTRING"), density = 1/10)
    pts[[i]] = pts[[i]] %>% st_cast("POINT")
    pts[[i]] = st_coordinates(pts[[i]])
    pts[[i]] = as.data.frame(pts[[i]])
    pts[[i]]$depth = pull(profile_cropped[i, "depth"] %>% st_drop_geometry())
  }
  pts = do.call(rbind, pts)
  colnames(pts) <- c("x", "y", "depth")
  pts = pts[!is.na(pts$x) & !is.na(pts$y) & !is.na(pts$depth),]
  
  # depth model linear interpolation: 
  gs <- gstat(formula = depth ~ 1, locations= ~x+y, data = pts, nmax = 50, set = list(idp = 0))
  
  # sample evaluation points
  eval_points = st_sample(st_bbox(profile_cropped), size = 100000, type = "regular")
  eval_points = as.data.frame(st_coordinates(eval_points))
  colnames(eval_points) <- c("x", "y")
  
  # predict depth for new points
  pred = predict(gs, eval_points)
  
  example_trajectory = ggplot() +
    geom_raster(data = pred, aes(x = x, y = y, fill = (-var1.pred)), alpha = .4, interpolate = T) +
    scale_fill_distiller(type = "seq",
                         direction = -1,
                         palette = "Greys",
                         guide = "none")+
    geom_sf(data = profile_cropped, color ="black", alpha = .2) +
    geom_path(data = subset, aes(x = x, y = y), color = "black", alpha = .4) +
    geom_point(data = subset %>% filter(Return == 0), aes(x = x, y = y, color = "Unsuccessful Spot  "), stroke = 1, size = 2, alpha = .7) +
    geom_point(data = subset %>% filter(Return == 1), aes(x = x, y = y, color = "Successful Spot"), stroke = 1, size = 2, alpha = .7) +
    ylab("Northing")+
    xlab("Easting")+
    labs(tag = "F")+
    scale_color_manual(values = colors)+
    theme(panel.background = element_rect(fill = NA),
          legend.key.size = unit(.3, "cm"),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linetype = "dashed", size = .4),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", size = .2)) + 
    annotation_raster(relocation, xmin = zoom[1] + 30, xmax = zoom[1] + 130, ymin = zoom[3] + 65, ymax = zoom[3] + 65 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50)) + 
    geom_segment(aes(x = zoom[1] + 30 + 50, y = zoom[3] + 65 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50), xend = zoom[1] + 130 + 10, yend = zoom[3] + 200 + 40),
                 arrow = arrow(length = unit(0.15, "cm")), color = "black", size = .8)+
    annotation_raster(angling, xmin = zoom[1] + 250, xmax = zoom[1] + 350, ymin = zoom[3] + 200, ymax = zoom[3] + 200 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50)) +
    geom_segment(aes(x = zoom[1] + 250 + 50, y = zoom[3] + 200, xend = coords_arrows$x+10, yend = coords_arrows$y+10),
                 arrow = arrow(length = unit(0.15, "cm")), color = "black", size = .8) + 
    annotation_scale()
  example_trajectory
  
  
  example_trajectory_no_shade = ggplot() +
    geom_sf(data = profile_cropped, color ="grey", alpha = 1) +
    geom_path(data = subset, aes(x = x, y = y), color = "black", alpha = .4) +
    geom_point(data = subset %>% filter(Return == 0), aes(x = x, y = y, color = "Unsuccessful Spot  "), stroke = 1, size = 2, alpha = .7) +
    geom_point(data = subset %>% filter(Return == 1), aes(x = x, y = y, color = "Successful Spot"), stroke = 1, size = 2, alpha = .7) +
    ylab("Northing")+
    xlab("Easting")+
    labs(tag = "E")+
    scale_color_manual(values = colors)+
    theme(panel.background = element_rect(fill = NA),
          legend.key.size = unit(.3, "cm"),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA, linetype = "dashed", size = .4),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white", size = .2)) + 
    annotation_raster(relocation, xmin = zoom[1] + 30, xmax = zoom[1] + 130, ymin = zoom[3] + 65, ymax = zoom[3] + 65 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50)) + 
    geom_segment(aes(x = zoom[1] + 30 + 50, y = zoom[3] + 65 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50), xend = zoom[1] + 130 + 10, yend = zoom[3] + 200 + 40),
                 arrow = arrow(length = unit(0.15, "cm")), color = "black", size = .8)+
    annotation_raster(angling, xmin = zoom[1] + 250, xmax = zoom[1] + 350, ymin = zoom[3] + 200, ymax = zoom[3] + 200 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50)) +
    geom_segment(aes(x = zoom[1] + 250 + 50, y = zoom[3] + 200, xend = coords_arrows$x+10, yend = coords_arrows$y+10),
                 arrow = arrow(length = unit(0.15, "cm")), color = "black", size = .8) + 
    annotation_scale()
  example_trajectory_no_shade
  
  # angling plot
  angling_sub = patch_leaving_data %>% 
    filter(day == 5 & year == 2022 & camera_id == 5) %>% 
    filter(sequence == 50) %>% 
    mutate(cumulative_r = cumsum(reward)) %>% 
    mutate(bin = bin - 1) %>% 
    mutate(bin = 10 * bin / 60)
  
  # load catch/start/leave pictures
  start <- readPNG("utils/data/raw_data/plot_images/start.png")
  catch1 <- readPNG("utils/data/raw_data/plot_images/catch1.png")
  catch2 <- readPNG("utils/data/raw_data/plot_images/catch2.png")
  leave <- readPNG("utils/data/raw_data/plot_images/leave.png")
  
  angling_plot = ggplot() + 
    geom_point(data = angling_sub %>% filter(reward == 1), aes(x = bin, y = cumulative_r)) + 
    geom_path(data = angling_sub, aes(x = bin, y = cumulative_r)) + 
    geom_vline(data = angling_sub %>% filter(reward == 1), aes(xintercept = bin), linetype = "dotted", size = .6) + 
    geom_vline(data = angling_sub %>% filter(bin == 0 | bin == max(bin)), aes(xintercept = bin), size = .6) + 
    ylab("Cumulative Catch (n)") + 
    xlab("Angling Time (min)") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 5),
          axis.text.y = element_text(size = 5),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white")) + 
    coord_fixed(ratio = .75)+
    scale_x_continuous(breaks = seq(from = 0, to = 10, by = 1))+
    scale_y_continuous(breaks = seq(from = 0, to = 10, by = 1))+
    labs(tag = "F")+
    geom_rect(aes(xmin = 0-.2, xmax = 0+.2, ymin = 7-.2*1/.75, ymax = 7+.2*1/.75), color = "black",fill = "white", size = .6)+
    geom_rect(aes(xmin = 8.5-.2, xmax = 8.5+.2, ymin = 7-.2*1/.75, ymax = 7+.2*1/.75), color = "black",fill = "white", size = .6)+
    annotate("text", x = 0, y=7, label = "S", size = 7/.pt) +
    add_fishape(family = "Percidae",option = "Perca_flavescens",
                xmin = .6,xmax = 1.4,
                ymin = 6.5, fill = "cadetblue",alpha = .9)+
    add_fishape(family = "Percidae",option = "Perca_flavescens",
                xmin = .93,xmax = 1.73,
                ymin = 6.5, fill = "cadetblue",alpha = .9)+
    add_fishape(family = "Percidae",option = "Perca_flavescens",
                xmin = 2.17-.4,xmax = 2.17+.4,
                ymin = 6.5, fill = "cadetblue",alpha = .9)+
    add_fishape(family = "Percidae",option = "Perca_flavescens",
                xmin = 3.83-.4,xmax = 3.83+.4,
                ymin = 6.5, fill = "cadetblue",alpha = .9)+
    add_fishape(family = "Percidae",option = "Perca_flavescens",
                xmin = 4.67-.4,xmax = 4.67+.4,
                ymin = 6.5, fill = "cadetblue",alpha = .9)+
    add_fishape(family = "Percidae",option = "Perca_flavescens",
                xmin = 7.33-.4,xmax = 7.33+.4,
                ymin = 6.5, fill = "cadetblue",alpha = .9)+
    annotate("text", x = 8.5, y=7, label = "L", size = 7/.pt) +
    annotation_raster(start, xmin = .1, xmax = .9, ymin = 1.8, ymax = 1.8 + 1/.75*dim(start)[1]/dim(start)[2] * (.9-.1)) + 
    geom_segment(aes(x = .1 + (.9-.1) / 2, y = 1.8, xend = 0, yend = 0),
                 arrow = arrow(length = unit(0.15, "cm")), color = "black", size = .6)+
    annotation_raster(catch1, xmin = 2.4, xmax = 3.2, ymin = 3.8, ymax = 3.8 + 1/.75*dim(start)[1]/dim(start)[2] * (.9-.1)) + 
    geom_segment(aes(x = 2.4 + (.9-.1) / 2, y = 3.8, xend = 2.17, yend = 3),
                 arrow = arrow(length = unit(0.15, "cm")), color = "black", size = .6)+
    annotation_raster(catch2, xmin = 5, xmax = 5.8, ymin = 3, ymax = 3 + 1/.75*dim(start)[1]/dim(start)[2] * (.9-.1)) + 
    geom_segment(aes(x = 5 + (.9-.1) / 2, y = 2.5 + 2*dim(start)[1]/dim(start)[2] * (.9-.1), xend = 4.67, yend = 5),
                 arrow = arrow(length = unit(0.15, "cm")), color = "black", size = .6)+
    annotation_raster(leave, xmin = 7.5, xmax = 8.3, ymin = 4, ymax = 4 + 1/.75*dim(start)[1]/dim(start)[2] * (.9-.1)) + 
    geom_segment(aes(x = 7.5 + (.9-.1) / 2, y = 3.5 + 2*dim(start)[1]/dim(start)[2] * (.9-.1), xend = 8.5, yend = 6),
                 arrow = arrow(length = unit(0.15, "cm")), color = "black", size = .6) 
  angling_plot
  
  
  left <- readPNG("utils/data/raw_data/plot_images/left.png")
  middle <- readPNG("utils/data/raw_data/plot_images/middle.png")
  right <- readPNG("utils/data/raw_data/plot_images/right.png")
  
  # plot with picture as layer
  ggplot(mapping = aes(1:10, 1:10)) +
    annotation_raster(left, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    geom_point()
  
  # make facet from each picture
  left_im <- ggplot(mapping = aes(x = 1:max(dim(left)), y = 1:max(dim(left)))) +
    annotation_raster(left, xmin = 1, xmax = dim(left)[2], ymin = 1, ymax = dim(left)[1]) +
    geom_point(alpha = 0) + 
    coord_fixed(ratio = 1, xlim = c(1, dim(left)[2]),ylim = c(1, dim(left)[1]), expand = 0) + 
    labs(tag = "A") + 
    theme(plot.tag = element_text(size = 7, face = "bold"),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x =element_blank(),
          axis.ticks.y =element_blank(),
          panel.border = element_rect(colour = "black", fill=NA))
  middle_im <- ggplot(mapping = aes(x = 1:max(dim(middle)), y = 1:max(dim(middle)))) +
    annotation_raster(middle, xmin = 1, xmax = dim(middle)[2], ymin = 1, ymax = dim(middle)[1]) +
    geom_point(alpha = 0) + 
    coord_fixed(ratio = 1, xlim = c(1, dim(middle)[2]),ylim = c(1, dim(middle)[1]), expand = 0) + 
    labs(tag = "B") + 
    theme(plot.tag = element_text(size = 7, face = "bold"),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x =element_blank(),
          axis.ticks.y =element_blank(),
          panel.border = element_rect(colour = "black", fill=NA))
  right_im <- ggplot(mapping = aes(x = 1:max(dim(right)), y = 1:max(dim(right)))) +
    annotation_raster(right, xmin = 1, xmax = dim(right)[2], ymin = 1, ymax = dim(right)[1]) +
    geom_point(alpha = 0) + 
    coord_fixed(ratio = 1, xlim = c(1, dim(right)[2]),ylim = c(1, dim(right)[1]), expand = 0) + 
    labs(tag = "C") + 
    theme(plot.tag = element_text(size = 7, face = "bold"),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x =element_blank(),
          axis.ticks.y =element_blank(),
          panel.border = element_rect(colour = "black", fill=NA))
  left_im
  middle_im
  right_im
  
  ################################################################################
  # Make Plots of movement trajectories for 3-4 lakes w/ one highlighted foraging trip
  day5 = gps_data %>% filter(day == 5 & year == 2022)
  
  # project coordinates to XY
  xy = LongLatToUTM(day5$position_long, day5$position_lat, zone = 35)
  day5$X = xy[, 2]
  day5$Y = xy[, 3] 
  day5 = st_as_sf(day5, coords = c("X", "Y"))
  st_crs(day5) <- "epsg:3067"
  day5$X = xy[, 2]
  day5$Y = xy[, 3]
  
  # project depth profile in common crs
  a = competition_area[[9]]
  d = st_transform(depth_profiles[[9]], crs = "epsg:3067")
  
  profile_cropped1 = st_crop(depth_profiles[[9]], xmin = min(day5$X) - 15, xmax = max(day5$X) + 15,
                             ymin = min(day5$Y) - 15, ymax = max(day5$Y) + 15)
  area_cropped1 = st_crop(a, xmin = min(day5$X) - 15, xmax = max(day5$X) + 15,
                          ymin = min(day5$Y) - 15, ymax = max(day5$Y) + 15)
  profile_cropped1 = st_intersection(st_transform(profile_cropped1, crs = crs(area_cropped1)),area_cropped1)
  
  # sample points from lakebed
  pts = list()
  profile_cropped1 = profile_cropped1 %>% st_cast("LINESTRING")
  for (i in 1:nrow(profile_cropped1)){
    pts[[i]] = st_line_sample(profile_cropped1[i,], density = 1/10)
    pts[[i]] = pts[[i]] %>% st_cast("POINT")
    pts[[i]] = st_coordinates(pts[[i]])
    pts[[i]] = as.data.frame(pts[[i]])
    pts[[i]]$depth = pull(profile_cropped1[i, "depth"] %>% st_drop_geometry())
  }
  pts = do.call(rbind, pts)
  colnames(pts) <- c("x", "y", "depth")
  pts = pts[!is.na(pts$x) & !is.na(pts$y) & !is.na(pts$depth),]
  
  # depth model linear interpolation: 
  library(gstat)
  gs1 <- gstat(formula = depth ~ 1, locations= ~x+y, data = pts, nmax = 50, set = list(idp = 0))
  
  # sample evaluation points
  eval_points1 = st_sample(area_cropped1, size = 100000, type = "regular")
  eval_points1 = as.data.frame(st_coordinates(eval_points1))
  colnames(eval_points1) <- c("x", "y")
  
  # predict depth for new points
  pred1 = predict(gs1, eval_points1)
  
  # plot
  lake1 = ggplot() + 
    geom_raster(data = pred1, aes(x = x, y = y, fill = (-var1.pred)), alpha = .8, interpolate = T) +
    #scale_fill_distiller(type = "seq",
    #                     direction = -1,
    #                     palette = "Greys")+
    scale_fill_gradient(low = "darkgrey", high = "white")+
    geom_sf(data =profile_cropped1, color ="black", alpha = .3) +
    #labs(tag = "iii)")+
    geom_path(data = day5, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .7) + 
    geom_path(data = day5 %>% filter(camera_id == 5), aes(x = X, y = Y, group = camera_id), color = "black", alpha = .7, size = .8) + 
    geom_rect(aes(xmin = rectangle[1]-15, xmax = rectangle[2]+15, ymin = rectangle[3]-15, ymax = rectangle[4]+15), fill = NA, color = "black", linetype = "dashed", size = .4)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  lake1
  
  lake1_no_shade = ggplot() + 
    #geom_raster(data = pred1, aes(x = x, y = y, fill = (-var1.pred)), alpha = .8, interpolate = T) +
    #scale_fill_distiller(type = "seq",
    #                     direction = -1,
    #                     palette = "Greys")+
    #scale_fill_gradient(low = "darkgrey", high = "white")+
    geom_sf(data =profile_cropped1, color ="grey", alpha = 1) +
    #labs(tag = "iii)")+
    geom_path(data = day5, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .7) + 
    geom_path(data = day5 %>% filter(camera_id == 5), aes(x = X, y = Y, group = camera_id), color = "black", alpha = .7, size = .8) + 
    geom_rect(aes(xmin = rectangle[1]-15, xmax = rectangle[2]+15, ymin = rectangle[3]-15, ymax = rectangle[4]+15), fill = NA, color = "black", linetype = "dashed", size = .4)+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  lake1_no_shade
  
  
  lake1_no_shade_alt = ggplot() + 
    #geom_raster(data = pred1, aes(x = x, y = y, fill = (-var1.pred)), alpha = .8, interpolate = T) +
    #scale_fill_distiller(type = "seq",
    #                     direction = -1,
    #                     palette = "Greys")+
    #scale_fill_gradient(low = "darkgrey", high = "white")+
    geom_sf(data =profile_cropped1, color ="grey", alpha = 1) +
    #labs(tag = "iii)")+
    geom_path(data = day5, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .7) + 
    geom_path(data = day5 %>% filter(camera_id == 5), aes(x = X, y = Y, group = camera_id), color = "black", alpha = .7, size = .8) + 
    #geom_rect(aes(xmin = rectangle[1]-15, xmax = rectangle[2]+15, ymin = rectangle[3]-15, ymax = rectangle[4]+15), fill = NA, color = "black", linetype = "dashed", size = .4)+
    theme(legend.key.size = unit(.3, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'))
  lake1_no_shade_alt
  
  # spots
  subset = spot_selection_data %>% 
    filter(day == 5 & year == 2022) %>% 
    filter(camera_id == 5)
  xy = subset[, c("x", "y")]
  subset = st_as_sf(subset, coords = c("x", "y"))
  st_crs(subset) <- "epsg:3067"
  subset$x = xy$x
  subset$y = xy$y
  
  example_trajectory_no_shade_alt = ggplot() +
    geom_sf(data = profile_cropped1, color ="grey", alpha = 1) +
    geom_path(data = day5, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .6, size = .3) + 
    geom_path(data = day5 %>% filter(camera_id == 5), aes(x = X, y = Y, group = camera_id), color = "black", alpha = .6, size = .6) + 
    geom_path(data = subset, aes(x = x, y = y), color = "black", alpha = .2) +
    geom_point(data = subset, aes(x = x, y = y), stroke = 1, size = 2, alpha = .7, color = "darkgreen") +
    theme(legend.key.size = unit(.3, "cm"),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = c(.8, .2),
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          panel.background = element_rect(fill='transparent'), #transparent panel bg
          plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
          panel.grid.major = element_blank(), #remove major gridlines
          panel.grid.minor = element_blank(), #remove minor gridlines
          legend.background = element_rect(fill='transparent'), #transparent legend bg
          legend.box.background = element_rect(fill='transparent'))
  example_trajectory_no_shade_alt
  
  ggsave(example_trajectory_no_shade_alt, file = "1_icefishing_social_foraging/output/lake_trajectory.png", 
         width = 12, height = 12,bg = "transparent", dpi = 700, units = "cm")
  ggsave(lake1_no_shade_alt, file = "1_icefishing_social_foraging/output/lake_trajectories.png", 
         width = 12, height = 12,bg = "transparent", dpi = 700, units = "cm")
  
  ################################################################################
  
  ################################################################################
  day4 = gps_data %>% filter(day == 4 & year == 2022)
  
  # project coordinates to XY
  xy = LongLatToUTM(day4$position_long, day4$position_lat, zone = 35)
  day4$X = xy[, 2]
  day4$Y = xy[, 3] 
  day4 = st_as_sf(day4, coords = c("X", "Y"))
  st_crs(day4) <- "epsg:3067"
  day4$X = xy[, 2]
  day4$Y = xy[, 3]
  
  # project depth profile in common crs
  a = competition_area[[10]]
  d = st_transform(depth_profiles[[10]], crs = "epsg:3067")
  
  # example traj
  profile_cropped2 = st_crop(depth_profiles[[10]], xmin = min(day4$X) - 15, xmax = max(day4$X) + 15,
                             ymin = min(day4$Y) - 15, ymax = max(day4$Y) + 15)
  area_cropped2 = st_crop(a, xmin = min(day4$X) - 15, xmax = max(day4$X) + 15,
                          ymin = min(day4$Y) - 15, ymax = max(day4$Y) + 15)
  
  profile_cropped2 = st_intersection(st_transform(profile_cropped2, crs = crs(area_cropped2)),area_cropped2)
  
  # sample points from lakebed
  # sample points from lakebed
  pts = list()
  for (i in 1:nrow(profile_cropped2)){
    pts[[i]] = st_line_sample(profile_cropped2[i,] %>% st_cast("LINESTRING"), density = 1/10)
    pts[[i]] = pts[[i]] %>% st_cast("POINT")
    pts[[i]] = st_coordinates(pts[[i]])
    pts[[i]] = as.data.frame(pts[[i]])
    pts[[i]]$depth = pull(profile_cropped2[i, "depth"] %>% st_drop_geometry())
  }
  pts = do.call(rbind, pts)
  colnames(pts) <- c("x", "y", "depth")
  pts = pts[!is.na(pts$x) & !is.na(pts$y) & !is.na(pts$depth),]
  
  # depth model linear interpolation: 
  gs2 <- gstat(formula = depth ~ 1, locations= ~x+y, data = pts, nmax = 50, set = list(idp = 0))
  
  # sample evaluation points
  eval_points2 = st_sample(area_cropped2, size = 100000, type = "regular")
  eval_points2 = as.data.frame(st_coordinates(eval_points2))
  colnames(eval_points2) <- c("x", "y")
  
  # predict depth for new points
  pred2 = predict(gs2, eval_points2)
  
  # visualize
  lake2 = ggplot() +
    geom_raster(data = pred2, aes(x = x, y = y, fill = (-var1.pred)), alpha = .8, interpolate = T) +
    #scale_fill_distiller(type = "seq",
    #                     direction = -1,
    #                     palette = "Greys")+
    scale_fill_gradient(low = "darkgrey", high = "white")+
    geom_sf(data =profile_cropped2, color ="black", alpha = .3) +
    #labs(tag = "iv)")+
    geom_path(data = day4, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .7, size = .8) + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  lake2
  
  lake2_no_shade = ggplot() +
    #geom_raster(data = pred2, aes(x = x, y = y, fill = (-var1.pred)), alpha = .8, interpolate = T) +
    #scale_fill_distiller(type = "seq",
    #                     direction = -1,
    #                     palette = "Greys")+
    #scale_fill_gradient(low = "darkgrey", high = "white")+
    geom_sf(data =profile_cropped2, color ="grey", alpha = 1) +
    #labs(tag = "iv)")+
    geom_path(data = day4, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .7, size = .8) + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  lake2_no_shade
  ################################################################################
  day1 = gps_data %>% filter(day == 1 & year == 2022)
  
  # project coordinates to XY
  xy = LongLatToUTM(day1$position_long, day1$position_lat, zone = 35)
  day1$X = xy[, 2]
  day1$Y = xy[, 3] 
  day1 = st_as_sf(day1, coords = c("X", "Y"))
  st_crs(day1) <- "epsg:3067"
  day1$X = xy[, 2]
  day1$Y = xy[, 3]
  
  # project depth profile in common crs
  a = competition_area[[6]]
  d = st_transform(depth_profiles[[6]], crs = "epsg:3067")
  
  # example traj
  profile_cropped3 = st_crop(d, xmin = min(day1$X) - 15, xmax = max(day1$X) + 15,
                             ymin = min(day1$Y) - 15, ymax = max(day1$Y) + 15)
  area_cropped3 = st_crop(a, xmin = min(day1$X) - 15, xmax = max(day1$X) + 15,
                          ymin = min(day1$Y) - 15, ymax = max(day1$Y) + 15)
  profile_cropped3 = st_intersection(st_transform(profile_cropped3, crs = crs(area_cropped3)),area_cropped3)
  
  # sample points from lakebed
  pts = list()
  for (i in 1:nrow(profile_cropped3)){
    pts[[i]] = st_line_sample(profile_cropped3[i,] %>% st_cast("LINESTRING"), density = 1/10)
    pts[[i]] = pts[[i]] %>% st_cast("POINT")
    pts[[i]] = st_coordinates(pts[[i]])
    pts[[i]] = as.data.frame(pts[[i]])
    pts[[i]]$depth = pull(profile_cropped3[i, "depth"] %>% st_drop_geometry())
  }
  pts = do.call(rbind, pts)
  colnames(pts) <- c("x", "y", "depth")
  pts = pts[!is.na(pts$x) & !is.na(pts$y) & !is.na(pts$depth),]
  
  # depth model linear interpolation: 
  gs3 <- gstat(formula = depth ~ 1, locations= ~x+y, data = pts, nmax = 50, set = list(idp = 0))
  
  # sample evaluation points
  eval_points3 = st_sample(area_cropped3, size = 100000, type = "regular")
  eval_points3 = as.data.frame(st_coordinates(eval_points3))
  colnames(eval_points3) <- c("x", "y")
  
  # predict depth for new points
  pred3 = predict(gs3, eval_points3)
  
  # visualize
  lake3 = ggplot() + 
    geom_raster(data = pred3, aes(x = x, y = y, fill = (-var1.pred)), alpha = .8, interpolate = T) +
    #scale_fill_distiller(type = "seq",
    #                     direction = -1,
    #                     palette = "Greys")+
    scale_fill_gradient(low = "darkgrey", high = "white")+
    #labs(tag = "ii)")+
    geom_sf(data = profile_cropped3, color ="black", alpha = .3) +
    geom_path(data = day1, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .6, size = .8) + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  lake3
  
  lake3_no_shade = ggplot() + 
    #geom_raster(data = pred3, aes(x = x, y = y, fill = (-var1.pred)), alpha = .8, interpolate = T) +
    #scale_fill_distiller(type = "seq",
    #                     direction = -1,
    #                     palette = "Greys")+
    #scale_fill_gradient(low = "darkgrey", high = "white")+
    #labs(tag = "ii)")+
    geom_sf(data = profile_cropped3, color ="grey", alpha = 1) +
    geom_path(data = day1, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .6, size = .8) + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  lake3_no_shade
  ################################################################################
  day10 = gps_data %>% filter(day == 5 & year == 2023)
  
  # project coordinates to XY
  xy = LongLatToUTM(day10$position_long, day10$position_lat, zone = 35)
  day10$X = xy[, 2]
  day10$Y = xy[, 3] 
  day10 = st_as_sf(day10, coords = c("X", "Y"))
  st_crs(day10) <- "epsg:3067"
  day10$X = xy[, 2]
  day10$Y = xy[, 3]
  
  # project depth profile in common crs
  a = competition_area[[5]]
  d = st_transform(depth_profiles[[5]], crs = "epsg:3067")
  #plot(depth_profiles[[5]])
  
  # example traj
  profile_cropped4 = st_crop(d, xmin = min(day10$X) - 15, xmax = max(day10$X) + 15,
                             ymin = min(day10$Y) - 15, ymax = max(day10$Y) + 15)
  area_cropped4 = st_crop(a, xmin = min(day10$X) - 15, xmax = max(day10$X) + 15,
                          ymin = min(day10$Y) - 15, ymax = max(day10$Y) + 15)
  profile_cropped4 = st_intersection(st_transform(profile_cropped4, crs = crs(area_cropped4)),area_cropped4)
  
  # sample points from lakebed
  pts = list()
  for (i in 1:nrow(profile_cropped4)){
    pts[[i]] = st_line_sample(profile_cropped4[i,] %>% st_cast("LINESTRING"), density = 1/10)
    pts[[i]] = pts[[i]] %>% st_cast("POINT")
    pts[[i]] = st_coordinates(pts[[i]])
    pts[[i]] = as.data.frame(pts[[i]])
    pts[[i]]$depth = pull(profile_cropped4[i, "depth"] %>% st_drop_geometry())
  }
  pts = do.call(rbind, pts)
  colnames(pts) <- c("x", "y", "depth")
  pts = pts[!is.na(pts$x) & !is.na(pts$y) & !is.na(pts$depth),]
  
  # depth model linear interpolation: 
  gs4 <- gstat(formula = depth ~ 1, locations= ~x+y, data = pts, nmax = 50, set = list(idp = 0))
  
  # sample evaluation points
  eval_points4 = st_sample(area_cropped4, size = 100000, type = "regular")
  eval_points4 = as.data.frame(st_coordinates(eval_points4))
  colnames(eval_points4) <- c("x", "y")
  
  # predict depth for new points
  pred4 = predict(gs4, eval_points4)
  
  # visualize
  lake4 = ggplot() + 
    geom_raster(data = pred4, aes(x = x, y = y, fill = (-var1.pred)), alpha = .8, interpolate = T) +
    #scale_fill_distiller(type = "seq",
    #                     direction = -1,
    #                     palette = "Greys")+
    scale_fill_gradient(low = "darkgrey", high = "white")+
    labs(tag="i)")+
    geom_sf(data = profile_cropped4, color ="black", alpha = .3) +
    geom_path(data = day10, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .7, size = .8) + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  lake4
  
  lake4_no_shade = ggplot() + 
    #geom_raster(data = pred4, aes(x = x, y = y, fill = (-var1.pred)), alpha = .8, interpolate = T) +
    #scale_fill_distiller(type = "seq",
    #                     direction = -1,
    #                     palette = "Greys")+
    #scale_fill_gradient(low = "darkgrey", high = "white")+
    #labs(tag="i)")+
    geom_sf(data = profile_cropped4, color ="grey", alpha = 1) +
    geom_path(data = day10, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .6, size = .8) + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "none",
          legend.text = element_text(size = 5),
          plot.tag = element_text(size = 7, face = "bold"),
          legend.title = element_blank(),
          legend.background = element_rect(colour = "black", fill = "white"))
  lake4_no_shade
  
  ################################################################################
  # example plot
  ################################################################################
  examples = readPNG("utils/data/raw_data/plot_images/lake_plots.png")
  example_plot <- ggplot(mapping = aes(x = 1:max(dim(examples)), y = 1:max(dim(examples)))) +
    annotation_raster(examples, xmin = 1, xmax = dim(examples)[2], ymin = 1, ymax = dim(examples)[1]) +
    #geom_point() + 
    coord_fixed(ratio = 1, xlim = c(1.1, dim(examples)[2]-.1),ylim = c(1.1, dim(examples)[1]-.1), expand = 0) + 
    labs(tag = "D") + 
    theme(plot.tag = element_text(size = 7, face = "bold"),
          plot.margin = margin(t = 0,  # Top margin
                               r = 0,  # Right margin
                               b = 0,  # Bottom margin
                               l = 0),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x =element_blank(),
          axis.ticks.y =element_blank())
  example_plot
  
  ################################################################################
  panel_figure1 = arrangeGrob(
    grobs = list(left_im,
                 middle_im,
                 right_im,
                 example_plot,
                 example_trajectory_no_shade,
                 angling_plot),#,
    #model_im),
    #widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    widths = c(4.3, 6.7, 7),
    layout_matrix = rbind(c(1, 4, 8),
                          c(1, 4, 8),
                          c(1, 4, 8),
                          c(1, 4, 8),
                          c(1, 4, 8),
                          c(2, 4, 8),
                          c(2, 4, 8),
                          c(2, 4, 8),
                          c(2, 4, 9),
                          c(2, 4, 9),
                          c(3, 4, 9),
                          c(3, 4, 9),
                          c(3, 4, 9),
                          c(3, 4, 9),
                          c(3, 4, 9))
  )
  
  
  # save figure in different file formats
  ggsave(panel_figure1, file = "1_icefishing_social_foraging/output/panel_figure1.png", 
         width = 18, height = 10, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_figure1, file = "1_icefishing_social_foraging/output/panel_figure1.svg", 
         width = 18, height = 10,bg = "white", dpi = 600, units = "cm",device = svglite::svglite)
  ggsave(panel_figure1, file = "1_icefishing_social_foraging/output/panel_figure1.pdf", 
         width = 18, height = 10,bg = "white", dpi = 600, units = "cm")
  ggsave(panel_figure1, file = "1_icefishing_social_foraging/output/panel_figure1.tiff",
         width = 18, height = 10,bg = "white", dpi = 600, units = "cm")
  ################################################################################
  # supplementary figure: lake locations (only works on windows bc. of OSM package and java)
  ################################################################################
  # load model data
  library(data.table)
  library(dplyr)
  
  # lake locations
  day = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
  year = c(2022, 2023, 2022, 2023, 2022, 2023, 2022, 2023, 2022, 2023)
  lat = c(62.62746, 62.59883, 62.56772, 62.36665, 62.66712, 62.35356, 62.75662, 62.62056, 62.68050, 62.84421)
  long = c(29.39998, 29.70948, 29.77903, 29.94114, 29.20515, 29.94607, 29.83124, 29.62019, 29.53324, 29.76650)
  site = c("Kuorinka (2022)", 
           "PyhÃ¤selka - Linnunlahti (2023)",
           "PyhÃ¤selka - Koivuniemi (2022)", 
           "Pieni-Onkamo (2023)", 
           "Viiniranta (2022)", 
           "Suuri-Onkamo (2023)", 
           "HÃ¶ytiÃ¤inen - Kontiolahti (2022)", 
           "PyhÃ¤selka - Lautasuo (2023)",
           "Keretinlahti (2022)", 
           "HÃ¶ytiÃ¤inen - Varparanta (2023)")
  lakes = data.frame(day = day, 
                     year = year,
                     lat = lat, 
                     long = long,
                     site = site)
  lakes$lat[8] = lakes$lat[8] - .01
  lakes$lat[6] = lakes$lat[6] - .01
  
  ################################################################################
  library(sp)
  lakes_sp = SpatialPointsDataFrame(coords = lakes[, c("long", "lat")], data = lakes,
                                    proj4string = CRS("+proj=longlat +zone=35 +ellps=WGS84 +units=m +no_defs"))
  plot(lakes_sp)
  
  bbox = bbox(lakes_sp)
  lakes_sp <- spTransform(lakes_sp,  CRS("+proj=merc +datum=WGS84 +units=m +no_defs"))
  
  Sys.setenv(JAVA_HOME="C:/Program Files/Java/jdk-24/")
  # java settings
  library(rJava)
  .jinit()  # Initialize Java VM with the new headless setting
  .jcall("java/awt/GraphicsEnvironment", "Z", "isHeadless")
  
  library(OpenStreetMap)
  map <- openmap(c(bbox[2,2]+.05, bbox[1,1]-.21), c(bbox[2,1]-.05, bbox[1,2]+.1),
                 type = "esri-imagery", mergeTiles = TRUE)
  # reproject onto WGS84
  map <- openproj(map, projection = "+proj=merc +datum=WGS84 +units=m +no_defs")# 
  
  ################################################################################
  library(OSMscale)
  library(ggspatial)
  library(ggplot2)
  plot_lakes <- OpenStreetMap::autoplot.OpenStreetMap(map) + 
    geom_point(data = data.frame(coordinates(lakes_sp)),
               aes(x = coords.x1, y = coords.x2), # slightly shift the points
               colour = "red", size =  1, shape = 1, stroke = 1) +
    xlab("") + ylab("")+
    annotation_scale(location = "bl", 
                     width_hint = 0.2, 
                     text_col = "white", 
                     pad_x = unit(1, "cm"), 
                     pad_y = unit(1, "cm"), 
                     height = unit(.3, "cm"), 
                     text_face = "bold", 
                     text_cex = .5) +
    coord_sf(crs = CRS("+proj=merc +datum=WGS84 +units=m +no_defs"))
  plot_lakes2 <- plot_lakes +
    geom_text(data = data.frame(coordinates(lakes_sp), site = lakes$site), # Choose dataframe
              aes(x = coords.x1, y = coords.x2, label = site), # Set aesthetics
              hjust = 1.1, vjust = 0.5, fontface = "bold",# Adjust vertical and horizontal
              size = 6 * 5/14, colour = "white") + # Adjust appearance
    theme_minimal() + 
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          plot.title = element_text(size = 7, 
                                    hjust = -0.0, 
                                    vjust=2.12, 
                                    face = "bold"),
          #aspect.ratio = 1/1,
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.spacing = unit(2, "lines"),
          legend.position = c(.6, .6),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7))
  ggsave(plot_lakes2, file = "1_icefishing_social_foraging/output/lake_locations.png", 
         width = 15, bg = "white", dpi = 600, units = "cm")
  
  ################################################################################
  # REVISION: Figure for structured summary
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
    filter(name %in% c("+Ecological", "+Personal", "+Conditional\nSocial Info")) %>% 
    mutate(improvement = elpd_waic-lead(elpd_waic)) %>% 
    filter(name != "+Ecological") %>% 
    ggplot(aes(x = fct_reorder(name, -model_id), y = improvement)) + 
    geom_col(color = "black", width = .6, aes(fill = name)) +
    scale_fill_manual(values = c("turquoise2", "cadetblue4"))+
    scale_x_discrete(labels = c("Personal\nInformation", "Social\nInformation"))+
    ggtitle("Spot-Selection") + 
    ylab("Importance") + 
    xlab("") + 
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  comparison_spot_selection
  
  comparison_basic = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/comparison/comparison_basic_strategies.rds")
  comparison_basic = as.data.frame(comparison_basic)
  comparison_basic$name = c("GUT","Fixed-Time","Fixed-n","Incremental", "Baseline")
  comparison_basic$model_id = c(2, 5, 4, 3, 1)
  comparison_extended = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/comparison/comparison_extended_strategies.rds")
  comparison_extended = as.data.frame(comparison_extended)
  comparison_extended$name = c("+Spatial\nFeatures", "+Updating", "+Fish\nDiscovery", "GUT")
  comparison_extended$model_id = 1:4
  comp_tmp = rbind(comparison_basic, comparison_extended)
  comp_tmp = comp_tmp[-1,]
  comparison_spot_leaving = comp_tmp %>%
    filter(name %in% c("Baseline", "+Updating", "+Spatial\nFeatures")) %>%
    arrange(-elpd_waic) %>% 
    mutate(improvement = elpd_waic-lead(elpd_waic)) %>% 
    filter(name != "Baseline") %>% 
    ggplot(aes(x = fct_reorder(name, -model_id), y = improvement)) + 
    geom_col(color = "black", width = .6, aes(fill = name)) +
    scale_fill_manual(values = c("turquoise2", "cadetblue4"))+
    scale_x_discrete(labels = c("Personal\nInformation", "Social\nInformation"))+
    ggtitle("Spot-Leaving") + 
    ylab("Importance") + 
    xlab("") + 
    #labs(tag = "a")+
    theme(panel.background = element_rect(fill = NA),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          axis.text.x = element_text(size = 6),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(colour = "black", fill=NA),
          legend.position = "none",
          plot.tag = element_text(size = 7, face = "bold"),
          title = element_text(size = 7),
          legend.text = element_text(size = 7),
          legend.title = element_text(size = 7),
          legend.background = element_rect(colour = "white", fill = "white"),
          text = element_text(family = "sans"))
  comparison_spot_leaving
  
  # make panel
  panel_tmp_comparison = arrangeGrob(
    grobs = list(comparison_spot_selection,
                 comparison_spot_leaving),
    widths = c(1),
    layout_matrix = rbind(c(1),
                          c(2))
  )
  
  # save figure
  ggsave(panel_tmp_comparison, file = "1_icefishing_social_foraging/output/model_comparison_fig0.png", 
         width = 4, height = 7.6, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_tmp_comparison, file = "1_icefishing_social_foraging/output/model_comparison_fig0.svg", 
         width = 4, height = 7, bg = "white", dpi = 600, units = "cm", device = svglite::svglite)
  ggsave(panel_tmp_comparison, file = "1_icefishing_social_foraging/output/model_comparison_fig0.pdf", 
         width = 4, height = 7, bg = "white", dpi = 600, units = "cm")
  ggsave(panel_tmp_comparison, file = "1_icefishing_social_foraging/output/model_comparison_fig0.tiff", 
         width = 4, height = 7, bg = "white", dpi = 600, units = "cm")
  
  ################################################################################
  # Survey analysis
  ################################################################################
  source("1_icefishing_social_foraging/scripts/survey_analysis.r")
  visualize_survey_data(survey_data, catch_data)
  
  ################################################################################
  # GPS error analysis
  ################################################################################
  source("1_icefishing_social_foraging/scripts/gps_error_analysis.r")
  visualize_gps_error()
  
  ################################################################################
  # Coding reliability analysis
  ################################################################################
  source("1_icefishing_social_foraging/scripts/video_reliability_analysis.r")
  analyse_video_reliability()
  
  ################################################################################
  # END
  ################################################################################
}
