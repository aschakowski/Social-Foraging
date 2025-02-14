################################################################################
#
# Title: Visualize spot selection simulation
#
# Description: Visualizes agend based model simulations.
#
# Authors: Schakowski, A.
#
# Last updated: 14/02/2025 (DD/MM/YYYY)
#
################################################################################
visualize_spot_selection_simulation <- function(){

# load ids
identifiers = readRDS(file = "utils/data/raw_data/identifiers_spot_selection.rds")
spot_selection_data = fread(file = "utils/data/raw_data/spot_selection_data.csv")

# add model ids
identifiers$lake = as.numeric(as.factor(identifiers$lake))
identifiers$trip = as.numeric(as.factor(identifiers$trip_id))
identifiers$id = as.numeric(as.factor(identifiers$participant_id))
  
# convert time
identifiers$time_since_start = as.numeric(identifiers$time_since_start) / max(as.numeric(identifiers$time_since_start))
  
# load sim results
random = readRDS(file = "utils/data/raw_data/spot_selection_simulations/sim_random.rds")
predictions = readRDS(file = "utils/data/raw_data/spot_selection_simulations/sim_predictions.rds")
nonsocial = readRDS(file = "utils/data/raw_data/spot_selection_simulations/sim_non_social.rds")
nonloss = readRDS(file = "utils/data/raw_data/spot_selection_simulations/sim_non_loss.rds")
nonsuccess = readRDS(file = "utils/data/raw_data/spot_selection_simulations/sim_non_success.rds")

# add time since last success to original data
identifiers$spot_id = as.numeric(identifiers$spot_id)
identifiers$ID = identifiers$camera_id
identifiers = left_join(identifiers, random %>% group_by(ID, day, year, spot_id) %>% slice(1) %>% dplyr::select(ID, day, year, spot_id, time_since_last_success))
identifiers$time_since_last_success = ifelse(is.na(identifiers$time_since_last_success), -1, identifiers$time_since_last_success)
identifiers$time_since_last_success = ifelse(identifiers$time_since_last_success > 0, identifiers$time_since_last_success - 1, identifiers$time_since_last_success)
identifiers$time_since_last_success2 = ifelse(identifiers$time_since_last_success > 10, 11, identifiers$time_since_last_success)

identifiers$day = as.numeric(identifiers$day)
identifiers$year = as.numeric(identifiers$year)
identifiers = left_join(identifiers, spot_selection_data %>% dplyr::select(ID, day, year, spot_id, step, angle))

random$time_since_last_success = ifelse(is.na(random$time_since_last_success), -1, random$time_since_last_success)
random$time_since_last_success = ifelse(random$time_since_last_success > 0, random$time_since_last_success - 1, random$time_since_last_success)
random$time_since_last_success2 = ifelse(random$time_since_last_success > 10, 11, random$time_since_last_success)

predictions$time_since_last_success = ifelse(is.na(predictions$time_since_last_success), -1, predictions$time_since_last_success)
predictions$time_since_last_success = ifelse(predictions$time_since_last_success > 0, predictions$time_since_last_success - 1, predictions$time_since_last_success)
predictions$time_since_last_success2 = ifelse(predictions$time_since_last_success > 10, 11, predictions$time_since_last_success)


nonsocial$time_since_last_success = ifelse(is.na(nonsocial$time_since_last_success), -1, nonsocial$time_since_last_success)
nonsocial$time_since_last_success = ifelse(nonsocial$time_since_last_success > 0, nonsocial$time_since_last_success - 1, nonsocial$time_since_last_success)
nonsocial$time_since_last_success2 = ifelse(nonsocial$time_since_last_success > 10, 11, nonsocial$time_since_last_success)

nonloss$time_since_last_success = ifelse(is.na(nonloss$time_since_last_success), -1, nonloss$time_since_last_success)
nonloss$time_since_last_success = ifelse(nonloss$time_since_last_success > 0, nonloss$time_since_last_success - 1, nonloss$time_since_last_success)
nonloss$time_since_last_success2 = ifelse(nonloss$time_since_last_success > 10, 11, nonloss$time_since_last_success)

nonsuccess$time_since_last_success = ifelse(is.na(nonsuccess$time_since_last_success), -1, nonsuccess$time_since_last_success)
nonsuccess$time_since_last_success = ifelse(nonsuccess$time_since_last_success > 0, nonsuccess$time_since_last_success - 1, nonsuccess$time_since_last_success)
nonsuccess$time_since_last_success2 = ifelse(nonsuccess$time_since_last_success > 10, 11, nonsuccess$time_since_last_success)

n_fun <- function(x){
  return(data.frame(y = median(x), label = paste0(length(x),"\n\n")))
}

median_se = function(z) {
  data.frame(
    y = median(z),
    ymin = median(z) - 2*(1.253 * (sd(z)/sqrt(length(z)))),
    ymax = median(z) + 2*(1.253 * (sd(z)/sqrt(length(z))))
  )
}


plot_model_data = identifiers %>% 
  filter(step > -99 & angle > -99) %>% 
  group_by(time_since_last_success2) %>% 
  mutate(n = n())

plot_random_data = random %>% 
  filter(step > -99) %>% 
  group_by(iter, time_since_last_success2) %>% 
  summarize(step = median((step), na.rm = T),
            angle = median(abs(angle), na.rm = T))

plot_brw_data = predictions %>% 
  filter(step > -99) %>% 
  group_by(iter,time_since_last_success2) %>% 
  summarize(step = median((step), na.rm = T),
            angle = median(abs(angle), na.rm = T))

# define colors
colors <- c("Model Simulation" = "darkgreen", "Data" = "black", "Random Walk" = "cadetblue")

# steps
steps = ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = .6)+
  stat_lineribbon(data = plot_random_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step), color = "Random Walk"),.width = .95, alpha = .7, fill = "lightblue")+
  stat_lineribbon(data = plot_brw_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step), color = "Model Simulation"), .width = .95, alpha = .7, fill = "darkseagreen2")+
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step), color = "Data", size = n), alpha = .2, fun.data = median_se, show.legend = F) + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step)), alpha = .6, geom = "line", fun.data = median_se, color = "darkgrey") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  scale_color_manual(values = colors) +
  scale_size(range=c(0.1, 2))+
  scale_y_continuous(breaks = seq(from = 0, to = 55, by = 10))+
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 1), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))+
  ylab("Median relocation distance (m)")+
  ggtitle("Area-restricted search")+
  coord_cartesian(xlim = c(-.5,11.5), ylim = c(0,40))+
  labs(tag = "a", color = "") +
  xlab("Successive unsuccessful spots (n)")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_blank(),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
steps
# angles
angles = ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = .6)+
  stat_lineribbon(data = plot_random_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (angle)*180/pi, color = "Random Walk"),.width = .95, alpha = .7, fill = "lightblue")+
  stat_lineribbon(data = plot_brw_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (angle)*180/pi, color = "Model Simulation"), .width = .95, alpha = .7, fill = "darkseagreen2")+
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = abs(angle*180/pi), size = n), alpha = .2, fun.data = median_se, show.legend = F, color = "black") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = abs(angle*180/pi)), alpha = .6, geom = "line", fun.data = median_se, color = "darkgrey") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = abs(angle*180/pi)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 1), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))+
  labs(tag = "b")+
  scale_size(range=c(0.1, 2))+
  ggtitle("Area-restricted search")+
  scale_y_continuous(breaks = seq(from = 0, to = 160, by = 20))+
  coord_cartesian(ylim = c(0, 120), 
                  xlim = c(-.5, 11.5))+
  ylab("Median turning angle (deg.)")+
  xlab("Successive unsuccessful spots (n)")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        title = element_blank(),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "bottom",
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        legend.direction = "horizontal",
        legend.key.size = unit(.3, "cm"),
        legend.margin = margin(5, 7, 5, 7))
angles
##################################################################################angle##############################################################################
# load fitted features
fitted_features = readRDS(file = "utils/data/raw_data/fitted_features.rds")
fitted_features$spot_id = as.numeric(fitted_features$spot_id)

# add to data and simulation
identifiers = left_join(identifiers, fitted_features %>% group_by(trip, spot_id) %>% slice(1) %>% unnest() %>% dplyr::select(day, year, spot_id, ID = camera_id, social_feature))

# median split of observed social densities
cutoff = median(identifiers$social_feature, na.rm = T)

################################################################################
colors_social <- c("High Social Density" = "darkolivegreen", "Low Social Density" = "darkseagreen2")

plot_model_data_social = identifiers %>% 
  group_by(time_since_last_success2, cond = social_feature > cutoff) %>% 
  mutate(n = n()) %>% 
  filter(step > -99 & angle > -99)

plot_random_data_social = random %>% 
  filter(step > -99) %>% 
  group_by(iter, time_since_last_success2, cond = social_feature > cutoff) %>% 
  summarize(step = median((step), na.rm = T),
            angle = median(abs(angle), na.rm = T))

plot_brw_data = predictions %>% 
  filter(step > -99) %>% 
  group_by(iter,time_since_last_success2, cond = social_feature > cutoff) %>% 
  summarize(step = median((step), na.rm = T),
            angle = median(abs(angle), na.rm = T))


steps_social = ggplot() +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = .6)+
  stat_lineribbon(data = plot_brw_data %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y =  step, color = "High Social Density"),fill = "darkseagreen3", .width = .95, alpha = .4)+
  stat_lineribbon(data = plot_brw_data %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y =  step, color = "Low Social Density"),fill = "darkseagreen1", .width = .95, alpha = .7)+
  
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = (step), color = "High Social Density", size = n), alpha = .4, fun.data = median_se, show.legend = F) + 
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = (step), color = "High Social Density"), alpha = .8, geom = "line", fun.data = median_se) + 
  
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2 + .1, y = (step), color = "Low Social Density", size = n), alpha = .6, fun.data = median_se, show.legend = F) + 
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2 + .1, y = (step), color = "Low Social Density"), alpha = .8, geom = "line", fun.data = median_se) + 
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2 + .1, y = (step)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = (step)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  
  
  scale_y_continuous(breaks = seq(from = 0, to = 55, by = 10))+
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 1), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))+
  ylab("Median relocation distance (m)")+
  xlab("Successive unsuccessful spots (n)")+
  coord_cartesian(xlim = c(-.5,11.5), ylim = c(0,40))+
  ggtitle("Social Area-Restricted Search")+
  scale_color_manual(values = colors_social)+
  labs(tag = "c", color = "")+
  scale_size(range=c(0.1, 2))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_blank(),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"))
#steps_social
angles_social = ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = .6)+
  stat_lineribbon(data = plot_brw_data %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = angle * (180/pi), color = "High Social Density"),fill = "darkseagreen3", .width = .95, alpha = .4)+
  stat_lineribbon(data = plot_brw_data %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = angle * (180/pi), color = "Low Social Density"),fill = "darkseagreen1", .width = .95, alpha = .7)+
  
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = abs(angle) * (180/pi), color = "High Social Density", size = n), alpha = .4, fun.data = median_se, show.legend = F) + 
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = abs(angle) * (180/pi), color = "High Social Density"), alpha = .8, geom = "line", fun.data = median_se) + 
  
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2 + .1, y = abs(angle) * (180/pi), color = "Low Social Density", size = n), alpha = .6, fun.data = median_se, show.legend = F) + 
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2 + .1, y = abs(angle) * (180/pi), color = "Low Social Density"), alpha = .8, geom = "line", fun.data = median_se) + 
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2 + .1, y = abs(angle) * (180/pi)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = abs(angle) * (180/pi)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  
  scale_color_manual(values = colors_social) +
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 1), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))+
  labs(tag = "d")+
  scale_size(range=c(0.1, 2))+
  ylab("Median turning angle (deg.)")+
  xlab("Successive unsuccessful spots (n)")+
  ggtitle("Social area-restricted search")+
  scale_y_continuous(breaks = seq(from = 0, to = 160, by = 20))+
  coord_cartesian(ylim = c(0, 120), 
                  xlim = c(-.5, 11.5))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        title = element_blank(),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "bottom",
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white", linewidth = .2),
        legend.direction = "horizontal",
        legend.key.size = unit(.3, "cm"),
        legend.margin = margin(5, 7, 5, 7))

# make panel
library(gridExtra)
panel_simulations = arrangeGrob(
  grobs = list(steps,
               angles,
               steps_social,
               angles_social),
  widths = c(.013, 1,.013,  1),
  heights = c(1, 1.2),
  layout_matrix = rbind(c(NA, 1,NA, 3),
                        c(2,2, 4, 4))
)

ggsave(panel_simulations, file = "2_3_socio_ecological_features/output/panel_simulations.png", 
       width = 16, height = 13, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_simulations, file = "2_3_socio_ecological_features/output/panel_simulations.svg", 
       width = 18, height = 13, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_simulations, file = "2_3_socio_ecological_features/output/panel_simulations.pdf", 
       width = 18, height = 13, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_simulations, file = "2_3_socio_ecological_features/output/panel_simulations.tiff", 
       width = 18, height = 13, bg = "white", dpi = 1200, units = "cm")

################################################################################
# supplementary figures: lesion plots
################################################################################
plot_brw_data_success_lesion = nonsuccess %>% 
  filter(step > -99) %>% 
  group_by(iter,time_since_last_success2) %>% 
  summarize(step = median((step), na.rm = T),
            angle = median(abs(angle), na.rm = T))

plot_brw_data_loss_lesion = nonloss %>% 
  filter(step > -99) %>% 
  group_by(iter,time_since_last_success2) %>% 
  summarize(step = median((step), na.rm = T),
            angle = median(abs(angle), na.rm = T))

# steps
success_lesion_steps = ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = .6)+
  stat_lineribbon(data = plot_random_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step), color = "Random Walk"),.width = .95, alpha = .7, fill = "lightblue")+
  stat_lineribbon(data = plot_brw_data_success_lesion %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step), color = "Model Simulation"), .width = .95, alpha = .7, fill = "darkseagreen2")+
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step), size = n), alpha = .2, fun.data = median_se, show.legend = F, color = "black") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step)), alpha = .6, geom = "line", fun.data = median_se, color = "darkgrey") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  scale_color_manual(values = colors) +
  scale_y_continuous(breaks = seq(from = 0, to = 55, by = 10))+
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 1), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))+
  ylab("Median relocation distance (m)")+
  ggtitle("Lesion Success")+
  scale_size(range=c(0.1, 2))+
  coord_cartesian(xlim = c(-.5,11.5), ylim = c(0,40))+
  labs(tag = "a", color = "") +
  xlab("Successive unsuccessful spots (n)")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_blank(),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.direction = "horizontal",
        legend.background = element_rect(colour = "black", fill = "white"))
success_lesion_steps

loss_lesion_steps = ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = .6)+
  stat_lineribbon(data = plot_random_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step), color = "Random Walk"),.width = .95, alpha = .7, fill = "lightblue")+
  stat_lineribbon(data = plot_brw_data_loss_lesion %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step), color = "Model Simulation"), .width = .95, alpha = .7, fill = "darkseagreen2")+
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step), size = n), alpha = .2, fun.data = median_se, show.legend = F, color = "black") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step)), alpha = .6, geom = "line", fun.data = median_se, color = "darkgrey") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (step)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  scale_color_manual(values = colors) +
  scale_y_continuous(breaks = seq(from = 0, to = 55, by = 10))+
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 1), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))+
  ylab("Median relocation distance (m)")+
  ggtitle("Lesion Unsuccessful")+
  scale_size(range=c(0.1, 2))+
  coord_cartesian(xlim = c(-.5,11.5), ylim = c(0,40))+
  labs(tag = "c", color = "") +
  xlab("Successive unsuccessful spots [n]")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_blank(),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
loss_lesion_steps

success_lesion_angles = ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = .6)+
  stat_lineribbon(data = plot_random_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (angle)*180/pi, color = "Random Walk"),.width = .95, alpha = .7, fill = "lightblue")+
  stat_lineribbon(data = plot_brw_data_success_lesion %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (angle)*180/pi, color = "Model Simulation"), .width = .95, alpha = .7, fill = "darkseagreen2")+
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = abs(angle*180/pi), size = n), alpha = .2, fun.data = median_se, show.legend = F, color = "black") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = abs(angle*180/pi)), alpha = .6, geom = "line", fun.data = median_se, color = "darkgrey") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = abs(angle*180/pi)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 1), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))+
  labs(tag = "b")+
  ggtitle("")+
  scale_size(range=c(0.1, 2))+
  scale_y_continuous(breaks = seq(from = 0, to = 160, by = 20))+
  coord_cartesian(ylim = c(0, 120), 
                  xlim = c(-.5, 11.5))+
  ylab("Median turning angle (deg.)")+
  xlab("Successive unsuccessful spots (n)")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_blank(),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))

loss_lesion_angles = ggplot() +
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = .6)+
  stat_lineribbon(data = plot_random_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (angle)*180/pi, color = "Random Walk"),.width = .95, alpha = .7, fill = "lightblue")+
  stat_lineribbon(data = plot_brw_data_loss_lesion %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = (angle)*180/pi, color = "Model Simulation"), .width = .95, alpha = .7, fill = "darkseagreen2")+
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = abs(angle*180/pi), size = n), alpha = .2, fun.data = median_se, show.legend = F, color = "black") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = abs(angle*180/pi)), alpha = .6, geom = "line", fun.data = median_se, color = "darkgrey") + 
  stat_summary(data = plot_model_data %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = abs(angle*180/pi)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  scale_color_manual(values = colors) +
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 1), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))+
  labs(tag = "d")+
  ggtitle("")+
  scale_size(range=c(0.1, 2))+
  scale_y_continuous(breaks = seq(from = 0, to = 160, by = 20))+
  coord_cartesian(ylim = c(0, 120), 
                  xlim = c(-.5, 11.5))+
  ylab("Median turning angle (deg.)")+
  xlab("Successive unsuccessful spots (n)")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_blank(),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.5, .15),
        legend.direction = "horizontal",
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"))

################################################################################
plot_brw_data_social_lesion = nonsocial %>% 
  filter(step > -99) %>% 
  group_by(iter,time_since_last_success2, cond = social_feature > cutoff) %>% 
  summarize(step = median((step), na.rm = T),
            angle = median(abs(angle), na.rm = T))

social_lesion_steps = ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = .6)+
  stat_lineribbon(data = plot_brw_data_social_lesion %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y =  step, color = "High Social Density"),fill = "darkseagreen3", .width = .95, alpha = .4)+
  stat_lineribbon(data = plot_brw_data_social_lesion %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y =  step, color = "Low Social Density"),fill = "darkseagreen1", .width = .95, alpha = .7)+
  
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = (step), color = "High Social Density", size = n), alpha = .4, fun.data = median_se, show.legend = F) + 
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = (step), color = "High Social Density"), alpha = .8, geom = "line", fun.data = median_se) + 
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = (step)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2+.1, y = (step), color = "Low Social Density", size = n), alpha = .6, fun.data = median_se, show.legend = F) + 
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2+.1, y = (step), color = "Low Social Density"), alpha = .8, geom = "line", fun.data = median_se) + 
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2+.1, y = (step)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  
  
  scale_y_continuous(breaks = seq(from = 0, to = 55, by = 10))+
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 1), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))+
  ylab("Median relocation distance (m)")+
  scale_size(range=c(0.1, 2))+
  xlab("Successive unsuccessful spots (n)")+
  coord_cartesian(xlim = c(-.5,11.5), ylim = c(0,40))+
  ggtitle("Lesion Social")+
  scale_color_manual(values = colors_social)+
  labs(tag = "e", color = "")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        legend.text = element_text(size = 6),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))

social_lesion_angles = ggplot() + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = .6)+
  stat_lineribbon(data = plot_brw_data_social_lesion %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = angle * (180/pi), color = "High Social Density"),fill = "darkseagreen3", .width = .95, alpha = .4)+
  stat_lineribbon(data = plot_brw_data_social_lesion %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2, y = angle * (180/pi), color = "Low Social Density"),fill = "darkseagreen2", .width = .95, alpha = .7)+
  
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = abs(angle) * (180/pi), color = "High Social Density", size = n), alpha = .4, fun.data = median_se, show.legend = F) + 
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = abs(angle) * (180/pi), color = "High Social Density"), alpha = .8, geom = "line", fun.data = median_se) + 
  stat_summary(data = plot_model_data_social %>% filter(cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2-.1, y = abs(angle) * (180/pi)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2+.1, y = abs(angle) * (180/pi), color = "Low Social Density", size = n), alpha = .6, fun.data = median_se, show.legend = F) + 
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2+.1, y = abs(angle) * (180/pi), color = "Low Social Density"), alpha = .8, geom = "line", fun.data = median_se) + 
  stat_summary(data = plot_model_data_social %>% filter(!cond) %>% filter(time_since_last_success2 >= 0), aes(x = time_since_last_success2+.1, y = abs(angle) * (180/pi)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  
  scale_color_manual(values = colors_social) +
  scale_x_continuous(breaks = seq(from = 0, to = 11, by = 1), labels = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10"))+
  labs(tag = "f")+
  ylab("Median turning angle (deg.)")+
  scale_size(range=c(0.1, 2))+
  xlab("Successive unsuccessful spots (n)")+
  ggtitle("")+
  scale_y_continuous(breaks = seq(from = 0, to = 160, by = 20))+
  coord_cartesian(ylim = c(0, 120), 
                  xlim = c(-.5, 11.5))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.5, .15),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"),
        legend.direction = "horizontal")

library(gridExtra)
panel_lesions = arrangeGrob(
  grobs = list(success_lesion_steps,
               loss_lesion_steps,
               social_lesion_steps,
               success_lesion_angles,
               loss_lesion_angles,
               social_lesion_angles),
  widths = c(1, 1),
  layout_matrix = rbind(c(1, 4),
                        c(2, 5), 
                        c(3, 6))
)
ggsave(panel_lesions, file = "2_3_socio_ecological_features/output/panel_lesions.png", 
       width = 18, height = 17, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_lesions, file = "2_3_socio_ecological_features/output/panel_lesions.svg", 
       width = 18, height = 17, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_lesions, file = "2_3_socio_ecological_features/output/panel_lesions.pdf", 
       width = 18, height = 17, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_lesions, file = "2_3_socio_ecological_features/output/panel_lesions.tiff", 
       width = 18, height = 17, bg = "white", dpi = 1200, units = "cm")

################################################################################
# explored area
# compute convex hull for each foraging trip
areas_data = spot_selection_data %>% 
  group_by(ID, day, year) %>% 
  mutate(n_spots = n()) %>% 
  filter(n_spots > 3) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            ratio = mean(total_catch > 0)) %>% 
  group_by(ID, day, year) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            ratio = mean(ratio)) %>% 
  mutate(area = area / (1000 * 1000)) %>% 
  mutate(ratio_binned = cut(ratio*100, breaks = c(0, 5, 10, 15, 20, 25, 30, Inf),
                            include.lowest = TRUE, right = TRUE)) %>% 
  group_by(ratio_binned) %>% 
  mutate(n = n())

areas_brw = predictions %>% 
  group_by(trip, iter) %>% 
  mutate(n_spots = n()) %>% 
  filter(n_spots > 3) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            ratio = mean(Return > 0)) %>% 
  group_by(trip, iter) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            ratio = mean(ratio)) %>% 
  mutate(area = area / (1000 * 1000)) %>% 
  mutate(ratio_binned = cut(ratio*100, breaks = c(0, 5, 10, 15, 20, 25, 30, Inf),
                            include.lowest = TRUE, right = TRUE)) %>% 
  group_by(iter, ratio_binned) %>% 
  summarize(area = median(area))

areas_random = random %>% 
  group_by(trip, iter) %>% 
  mutate(n_spots = n()) %>% 
  #filter(n_spots > 3) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            ratio = mean(Return > 0)) %>% 
  group_by(trip, iter) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            ratio = mean(ratio)) %>% 
  mutate(area = area / (1000 * 1000)) %>% 
  mutate(ratio_binned = cut(ratio*100, breaks = c(0, 5, 10, 15, 20, 25, 30, Inf),
                            include.lowest = TRUE, right = TRUE)) %>% 
  group_by(iter, ratio_binned) %>% 
  summarize(area = median(area))

areas_success_lesion = nonsuccess %>% 
  group_by(trip, iter) %>% 
  mutate(n_spots = n()) %>% 
  #filter(n_spots > 3) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            ratio = mean(Return > 0)) %>% 
  group_by(trip, iter) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            ratio = mean(ratio)) %>% 
  mutate(area = area / (1000 * 1000)) %>% 
  mutate(ratio_binned = cut(ratio*100, breaks = c(0, 5, 10, 15, 20, 25, 30, Inf),
                            include.lowest = TRUE, right = TRUE)) %>% 
  group_by(iter, ratio_binned) %>% 
  summarize(area = median(area))

areas_loss_lesion = nonloss %>% 
  group_by(trip, iter) %>% 
  mutate(n_spots = n()) %>% 
  #filter(n_spots > 3) %>% 
  summarize(area = (cbind(x[chull(x,y)], y[chull(x, y)])),
            ratio = mean(Return > 0)) %>% 
  group_by(trip, iter) %>% 
  mutate(n = n()) %>% 
  filter(n >= 4) %>% 
  summarize(area = Polygon(area)@area,
            ratio = mean(ratio)) %>% 
  mutate(area = area / (1000 * 1000)) %>% 
  mutate(ratio_binned = cut(ratio*100, breaks = c(0, 5, 10, 15, 20, 25, 30, Inf),
                            include.lowest = TRUE, right = TRUE)) %>% 
  group_by(iter, ratio_binned) %>% 
  summarize(area = median(area))


explored_area = ggplot() + 
  stat_lineribbon(data = areas_random, aes(x = ratio_binned, y = (area), color = "Random Walk"),.width = .95, alpha = .5, fill = "lightblue")+
  stat_lineribbon(data = areas_brw, aes(x = ratio_binned, y = (area), color = "Model Simulation"), .width = .95, alpha = .5, fill = "lightgreen")+
  stat_summary(data = areas_data, aes(x = ratio_binned, y = (area), size = n), alpha = .2, fun.data = median_se, show.legend = F, color = "black") + 
  stat_summary(data = areas_data, aes(x = ratio_binned, y = (area)), alpha = .6, geom = "line", fun.data = median_se, color = "darkgrey") + 
  stat_summary(data = areas_data, aes(x = ratio_binned, y = (area)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  scale_color_manual(values = colors) +
  scale_y_continuous(breaks = seq(from = 0, to = .3, by = .05))+
  ylab("Median explored area (sqkm)")+
  ggtitle("Full model")+
  scale_size(range=c(0.1, 2))+
  labs(tag = "a", color = "") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "lightgrey", size = 1)+
  xlab("% Successful spots")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 6, angle = 55, vjust = .5),
        plot.title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"))
explored_area

# lesion plots
explored_area_nosuccess = ggplot() + 
  stat_lineribbon(data = areas_random, aes(x = ratio_binned, y = (area), color = "Random Walk"),.width = .95, alpha = .5, fill = "lightblue")+
  stat_lineribbon(data = areas_success_lesion, aes(x = ratio_binned, y = (area), color = "Model Simulation"), .width = .95, alpha = .5, fill = "lightgreen")+
  stat_summary(data = areas_data, aes(x = ratio_binned, y = (area), size = n), alpha = .2, fun.data = median_se, show.legend = F, color = "black") + 
  stat_summary(data = areas_data, aes(x = ratio_binned, y = (area)), alpha = .6, geom = "line", fun.data = median_se, color = "darkgrey") + 
  stat_summary(data = areas_data, aes(x = ratio_binned, y = (area)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  scale_color_manual(values = colors) +
  scale_y_continuous(breaks = seq(from = 0, to = .3, by = .05))+
  ylab("Median explored area (sqkm)")+
  ggtitle("Lesion successful")+
  scale_size(range=c(0.1, 2))+
  labs(tag = "b", color = "") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "lightgrey", size = 1)+
  xlab("% Successful spots")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6, angle = 55, vjust = .5),
        plot.title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"))

explored_area_loss = ggplot() + 
  stat_lineribbon(data = areas_random, aes(x = ratio_binned, y = (area), color = "Random Walk"),.width = .95, alpha = .5, fill = "lightblue")+
  stat_lineribbon(data = areas_loss_lesion, aes(x = ratio_binned, y = (area), color = "Model Simulation"), .width = .95, alpha = .5, fill = "lightgreen")+
  stat_summary(data = areas_data, aes(x = ratio_binned, y = (area), size = n), alpha = .2, fun.data = median_se, show.legend = F, color = "black") + 
  stat_summary(data = areas_data, aes(x = ratio_binned, y = (area)), alpha = .6, geom = "line", fun.data = median_se, color = "darkgrey") + 
  stat_summary(data = areas_data, aes(x = ratio_binned, y = (area)), fun.data = n_fun, geom = "text", size = 5/.pt) + 
  scale_color_manual(values = colors) +
  scale_y_continuous(breaks = seq(from = 0, to = .3, by = .05))+
  ylab("Median explored area (sqkm)")+
  ggtitle("Lesion unsuccessful")+
  labs(tag = "c", color = "") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "lightgrey", size = 1)+
  xlab("% Successful spots")+
  scale_size(range=c(0.1, 2))+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x =  element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 6, angle = 55, vjust = .5),
        plot.title = element_text(size = 7),
        plot.tag = element_text(size = 7, face = "bold"),
        axis.text.y = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.65, .7),
        legend.text = element_text(size = 5),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"))

library(gridExtra)
panel_expl_area_lesions = arrangeGrob(
  grobs = list(explored_area,
               explored_area_nosuccess,
               explored_area_loss),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3))
)
ggsave(panel_expl_area_lesions, file = "2_3_socio_ecological_features/output/panel_expl_area_lesions.png", 
       width = 18, height = 6, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_expl_area_lesions, file = "2_3_socio_ecological_features/output/panel_expl_area_lesions.svg", 
       width = 18, height = 6, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_expl_area_lesions, file = "2_3_socio_ecological_features/output/panel_expl_area_lesions.pdf", 
       width = 18, height = 6, bg = "white", dpi = 1200, units = "cm")
ggsave(panel_expl_area_lesions, file = "2_3_socio_ecological_features/output/panel_expl_area_lesions.tiff", 
       width = 18, height = 6, bg = "white", dpi = 1200, units = "cm")

}
################################################################################
# END
################################################################################
