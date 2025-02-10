################################################################################
#
# Title: Descriptive Analyses and Figures
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 29/11/2024 (DD/MM/YYYY)
#
################################################################################
run_1_icefishing_social_foraging <- function(wd){

# set working directory
setwd(wd)

# check
getwd()

################################################################################
# load data
################################################################################

# dataset contains basic demographics, catch success for each competition and camera/watch ids
catch_data <- read_delim("utils/data/raw_data/catch_data.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
length(unique(catch_data$participant_id))

# data contains survey responses for all participants (also for participants in group competitions, not reported here)
survey_data <- read_delim("utils/data/raw_data/survey_data.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)
length(unique(survey_data$participant_id))

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
# how many participants/foraging trips per day
catch_data %>% 
  group_by(day, year) %>% 
  summarize(n = n())
length(unique(catch_data$participant_id))
nrow(catch_data)

# how many individuals fished before?
catch_data %>% 
  filter(!is.na(fished_before) & fished_before != "NA") %>% 
  mutate(fished_before = as.numeric(as.factor(fished_before))) %>%
  group_by(year) %>% 
  summarize(f = mean(as.numeric(as.factor(fished_before))-1, na.rm = T),
            n = n())

# how many of these on all days?
catch_data %>% 
  group_by(participant_id) %>% 
  summarize(n = n()) %>% 
  filter(n == 10)

# how many more than 5 days?
catch_data %>% 
  group_by(participant_id) %>% 
  summarize(n = n()) %>% 
  filter(n >= 5)

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
#participation
# save figure in different file formats
ggsave(participation, file = "1_icefishing_social_foraging/output/participation_plot.png", 
       width = 9, height = 5,bg = "white", dpi = 1200, units = "cm")
ggsave(participation, file = "1_icefishing_social_foraging/output/participation_plot.svg", 
       width = 9, height = 5,bg = "white", dpi = 1200, units = "cm")
ggsave(participation, file = "1_icefishing_social_foraging/output/participation_plot.pdf", 
       width = 9, height = 5,bg = "white", dpi = 1200, units = "cm")
ggsave(participation, file = "1_icefishing_social_foraging/output/participation_plot.tiff", 
       width = 9, height = 5,bg = "white", dpi = 1200, units = "cm")

# avg experience
mean(survey_data$years_icefishing, na.rm= T)
sd(survey_data$years_icefishing, na.rm = T)
sum(is.na(survey_data$years_icefishing))

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
  summarize(mean_sex = mean(as.numeric(as.factor(sex))-1, na.rm = T),
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
  group_by(day, year, ID) %>% 
  summarize(rel_angling_time = sum(angling_time) / max(time_end))
mean(rel_angling_time$rel_angling_time)
sd(rel_angling_time$rel_angling_time)

# how many spots visited during a competition?
n_spots_visited = spot_selection_data %>% 
  filter(!exclude) %>% 
  group_by(day, year, ID) %>% 
  summarize(n_spots = max(spot_id))
mean(n_spots_visited$n_spots)
sd(n_spots_visited$n_spots)

# total number of foraging locations analysed
nrow(spot_selection_data %>% filter(!exclude))

# network analysis with dbscan
radii = seq(from = 0, to = 500, by = 10)
n_iter = 50
source("1_icefishing_social_foraging/scripts/create_network_data.r")
create_network_data(gps_data, radii, n_iter)
iterated_results = readRDS(file = "utils/data/processed_data/network_data.rds")

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
  xlab(expression(epsilon~(m)))+
  ylab("Relative Community Size")+
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
        axis.ticks.y = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white")) + 
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
       width = 9, height = 9,bg = "white", dpi = 1200, units = "cm")
ggsave(networks, file = "1_icefishing_social_foraging/output/networks.svg", 
       width = 9, height = 9,bg = "white", dpi = 1200, units = "cm")
ggsave(networks, file = "1_icefishing_social_foraging/output/networks.pdf", 
       width = 9, height = 9,bg = "white", dpi = 1200, units = "cm")
ggsave(networks, file = "1_icefishing_social_foraging/output/networks.tiff", 
       width = 9, height = 9,bg = "white", dpi = 1200, units = "cm")

# how much time spent at angling spots
spot_selection_data %>% 
  filter(!exclude) %>% 
  summarize(mean(angling_time),sd(angling_time))

cor.test(spot_selection_data$total_catch[spot_selection_data$exclude == 0],
         spot_selection_data$angling_time[spot_selection_data$exclude == 0])

# variation in catch success across lakes
success = catch_data %>% 
  group_by(day, year) %>% 
  summarize(c = mean(catch), sd(catch)) 

# percentage successful spots
perc_success = spot_selection_data %>% 
  filter(!exclude) %>%
  group_by(day, year) %>% 
  summarize(r = mean(Return))

test = left_join(perc_success, success)
cor.test(test$r, test$c)
plot(test$r, test$c)

# how many successful spots with time series information
patch_leaving_data %>% 
  filter(!exclude2) %>% 
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
  group_by(year, day, ID) %>% 
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
  group_by(year, day, ID) %>% 
  summarize(n = n(),
            n_successful = mean(Return)) %>% 
  group_by(year, day) %>% 
  summarize(mean_n = mean(n),
            median_n = median(n),
            sd_n = sd(n),
            percent_successful = mean(n_successful),
            median_successful = median(n_successful),
            sd_successful = sd(n_successful))


patch_leaving_data %>% 
  group_by(year, day, camera_id, sequence) %>% 
  summarize(angling_time = n()*10,
            total_catch = sum(reward),
            missing = mean(as.numeric((exclude + exclude2)>0))) %>% 
  group_by(year, day) %>% 
  summarize(n_missing = sum(missing),
            percent_missing = mean(missing),
            n_spots = n())


patch_leaving_data %>% 
  group_by(year, day, camera_id, sequence) %>% 
  summarize(angling_time = n()*10,
            total_catch = sum(reward),
            missing = mean(as.numeric((exclude + exclude2)>0))) %>% 
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
            missing = mean(as.numeric((exclude + exclude2)>0))) %>%
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
            missing = mean(as.numeric((exclude + exclude2)>0))) %>%
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
            missing = mean(as.numeric((exclude + exclude2)>0))) %>%
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
  labs(tag = "a") +
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
catch_day
data_a = catch_data %>% 
  group_by(day, year, camera_id) %>% 
  mutate(catch=mean(catch/1000), group= interaction(day, year)) %>% 
  dplyr::select(catch, group)

ratio = spot_selection_data %>% 
  group_by(day, year) %>% 
  summarize(mean(Return))
min(ratio$`mean(Return)`)
max(ratio$`mean(Return)`)
median(ratio$`mean(Return)`)

ratios = spot_selection_data %>% 
  group_by(day, year, ID) %>% 
  summarize(r = mean(Return)) %>% 
  ggplot(aes(y=r*100, x = as.factor(as.numeric(interaction(day, year))))) + 
  geom_boxplot(outlier.shape = NA) + 
  ylab("% Successful Spots")+
  xlab("Lake")+
  labs(tag = "b") +
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
ratios
data_b = spot_selection_data %>% 
  group_by(day, year, ID) %>% 
  summarize(r = mean(Return)*100) 

spot_selection_data$ID = as.character(spot_selection_data$ID)
spot_selection_data = left_join(spot_selection_data, catch_data %>% dplyr::select(ID=camera_id, day, year, catch))
test = spot_selection_data %>% 
  group_by(day, year) %>% 
  summarize(r = mean(Return)*100,
            r_sd = sd(Return)*100,
            c = mean(catch)/1000,
            c_sd = sd(catch/1000),
            n = n())
cor.test(test$r, test$c)

test = test %>% 
  ungroup() %>% 
  mutate(scale_r = (r - mean(r))/sd(r),
         scale_c = (c - mean(c))/sd(c))

correl = brm(data = test, formula = scale(r)~1+scale(c), 
             chains =4,
             cores = 4)

correl_day = test %>% 
  ggplot(aes(x = r, y = c)) + 
  geom_point(size = .6) + 
  geom_errorbar(aes(ymin = c - 2*(c_sd/sqrt(n)), ymax = c + 2*(c_sd/sqrt(n)))) + 
  geom_errorbarh(aes(xmin = r - 2*(r_sd/sqrt(n)), xmax = r + 2*(r_sd/sqrt(n)))) + 
  #theme_classic() +
  labs(tag = "c") +
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
correl_day
data_c = test

# catch in n and kg correlated
test = spot_selection_data# %>% 
#  filter(total_catch > 0)

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
  labs(tag = "d") +
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
data_d1 = prediction_data
data_d2 = spot_selection_data

# catch in kg and n fish per day
n_dat = spot_selection_data %>% 
  group_by(day, year, camera_id = ID) %>% 
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


cor.test(n_dat$catch, n_dat$n)
catch_n_kg =  
  ggplot() + 
  stat_lineribbon(data = prediction_data_n, aes(x = n_pred, y = posterior_prediction/1000),.width = c(.66, .95), alpha = .4, size = .6)+
  geom_point(data = n_dat, aes(x = n, y = catch/1000), alpha = .2, size = .6)+
  scale_fill_manual(values = c("darkgrey","lightgrey"), guide = "none")+
  labs(tag = "e") +
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
        legend.background = element_rect(colour = "black", fill = "white"))# +
  #annotate("text", x = 200, y = 1, label = "r = .743 [.696, .783]", size = 6)
catch_n_kg
data_e1 = prediction_data_n
data_e2 = n_dat

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
       width = 18, height = 14,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_lake_characteristics, file = "1_icefishing_social_foraging/output/panel_lake_characteristics.svg", 
       width = 18, height = 14,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_lake_characteristics, file = "1_icefishing_social_foraging/output/panel_lake_characteristics.pdf", 
       width = 18, height = 14,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_lake_characteristics, file = "1_icefishing_social_foraging/output/panel_lake_characteristics.tiff", 
       width = 18, height = 14,bg = "white", dpi = 1200, units = "cm")

# save supplementary figures
time_allocation = video_data %>% 
  group_by(day, year, camera_id) %>% 
  summarize(angling = sum(angling, na.rm = T),
            relocating = sum(relocating, na.rm = T),
            drilling = sum(drilling, na.rm = T))

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
       width = 18, height = 9,bg = "white", dpi = 1200, units = "cm")
ggsave(time_allocation_plot, file = "1_icefishing_social_foraging/output/time_allocation_plot.svg", 
       width = 18, height = 9,bg = "white", dpi = 1200, units = "cm")
ggsave(time_allocation_plot, file = "1_icefishing_social_foraging/output/time_allocation_plot.pdf", 
       width = 18, height = 9,bg = "white", dpi = 1200, units = "cm")
ggsave(time_allocation_plot, file = "1_icefishing_social_foraging/output/time_allocation_plot.tiff", 
       width = 18, height = 9,bg = "white", dpi = 1200, units = "cm")

spots_visited = spot_selection_data %>% 
  group_by(ID, day, year) %>% 
  summarize(n = n()) %>% 
  ggplot(aes(x = n)) + 
  geom_histogram(color = "black") + 
  ylab("ID count")+
  xlab("Visited Foraging Locations (n)")+
  facet_wrap(~day*year, nrow = 2)+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))

# save figure in different file formats
ggsave(spots_visited, file = "1_icefishing_social_foraging/output/spots_visited.png", 
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")
ggsave(spots_visited, file = "1_icefishing_social_foraging/output/spots_visited.svg", 
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")
ggsave(spots_visited, file = "1_icefishing_social_foraging/output/spots_visited.pdf", 
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")
ggsave(spots_visited, file = "1_icefishing_social_foraging/output/spots_visited.tiff", 
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")

avg_relocation = 
  spot_selection_data %>% filter(step>0) %>% 
  group_by(ID, day, year) %>% 
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
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")
ggsave(avg_relocation, file = "1_icefishing_social_foraging/output/avg_relocation.svg", 
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")
ggsave(avg_relocation, file = "1_icefishing_social_foraging/output/avg_relocation.pdf", 
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")
ggsave(avg_relocation, file = "1_icefishing_social_foraging/output/avg_relocation.tiff",
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")

spot_selection_data %>% filter(step>0) %>% 
  group_by(ID, day, year) %>% 
  summarize(step = sum(step)) %>% 
  ungroup() %>% 
  summarize(m_step = mean(step), 
            sd = sd(step, na.rm = T))


total_relocation = 
  spot_selection_data %>% filter(step>0) %>% 
  group_by(ID, day, year) %>% 
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
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")
ggsave(total_relocation, file = "1_icefishing_social_foraging/output/total_relocation.svg", 
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")
ggsave(total_relocation, file = "1_icefishing_social_foraging/output/total_relocation.pdf", 
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")
ggsave(total_relocation, file = "1_icefishing_social_foraging/output/total_relocation.tiff",
       width = 18, height = 8,bg = "white", dpi = 1200, units = "cm")

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

# compute x and y coordinates#Function: https://stackoverflow.com/questions/18639967/converting-latitude-and-longitude-points-to-utm
LongLatToUTM<-function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
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
  filter(ID == 5)
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
  geom_path(data = day5 %>% filter(camera_id == 5), aes(x = X, y = Y, group = camera_id), color = "black", alpha = 1, size = .8) + 
  geom_rect(aes(xmin = rectangle[1]-15, xmax = rectangle[2]+15, ymin = rectangle[3]-15, ymax = rectangle[4]+15), fill = NA, color = "black", linetype = "dashed", size = 1)+
  ylab("Northing")+
  xlab("Easting")+
  labs(tag = "d")+
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
trajectories

# zoom into map
profile_cropped = st_crop(depth_profiles[[9]], xmin = min(subset$x) - 15, xmax = max(subset$x) + 15,
                          ymin = min(subset$y) - 15, ymax = max(subset$y) + 15)

# assign color values
colors <- c("Successful Spots" = "cyan3", "Unsuccessful Spots     " = "darkgreen")

# for fish icons
#devtools::install_github("mattiaghilardi/fishualize", ref = "fix-fishapes")
#library(fishualize)

# read images to display on plot
relocation <- readPNG("utils/data/raw_data/plot_images/relocating.png")
angling <- readPNG("utils/data/raw_data/plot_images/angling.png")

# coordinates for arrows in plot
coords_arrows = subset %>% filter(Return == 0) 
coords_arrows = coords_arrows[30, c("x", "y")] 

zoom = raster::extent(profile_cropped)
example_trajectory = ggplot() + 
  geom_sf(data = profile_cropped, color ="black", alpha = .2) +
  geom_path(data = subset, aes(x = x, y = y), color = "black", alpha = .4) +
  geom_point(data = subset %>% filter(Return == 1), aes(x = x, y = y, color = "Successful Spots"), stroke = 3, size = 3, alpha = .7) +
  geom_point(data = subset %>% filter(Return == 0), aes(x = x, y = y, color = "Unsuccessful Spots     "), stroke = 1, size = 2, alpha = .7) +
  ylab("Northing")+
  xlab("Easting")+
  labs(tag = "e")+
  scale_color_manual(values = colors)+
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
  annotation_raster(relocation, xmin = zoom[1] + 30, xmax = zoom[1] + 130, ymin = zoom[3] + 65, ymax = zoom[3] + 65 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50)) + 
  geom_segment(aes(x = zoom[1] + 30 + 50, y = zoom[3] + 65 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50), xend = zoom[1] + 130 + 10, yend = zoom[3] + 200 + 40),
               arrow = arrow(length = unit(0.15, "cm")), color = "black", size = .8)+
  annotation_raster(angling, xmin = zoom[1] + 250, xmax = zoom[1] + 350, ymin = zoom[3] + 200, ymax = zoom[3] + 200 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50)) +
  geom_segment(aes(x = zoom[1] + 250 + 50, y = zoom[3] + 200, xend = coords_arrows$x+10, yend = coords_arrows$y+10),
               arrow = arrow(length = unit(0.15, "cm")), color = "black", size = .8) + 
  annotation_scale()
example_trajectory

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
  labs(tag = "f")+
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
  labs(tag = "a") + 
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
  labs(tag = "b") + 
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
  labs(tag = "c") + 
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

# 3 example lakes for plot
################################################################################
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

lake1 = ggplot() + 
  geom_sf(data =profile_cropped1, color ="black", alpha = .3) +
  geom_path(data = day5, aes(x = X, y = Y, group = camera_id), color = "green4", alpha = .4) + 
  geom_path(data = day5 %>% filter(camera_id == 5), aes(x = X, y = Y, group = camera_id), color = "black", alpha = 1, size = .8) + 
  geom_rect(aes(xmin = rectangle[1]-15, xmax = rectangle[2]+15, ymin = rectangle[3]-15, ymax = rectangle[4]+15), fill = NA, color = "black", linetype = "dashed", size = 1)+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
lake1
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
lake2 = ggplot() + 
  geom_sf(data = profile_cropped2, color ="black", alpha = .3) +
  geom_path(data = day4, aes(x = X, y = Y, group = camera_id), color = "cadetblue", alpha = .4) + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
lake2
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
profile_cropped3 = st_crop(depth_profiles[[6]], xmin = min(day1$X) - 15, xmax = max(day1$X) + 15,
                          ymin = min(day1$Y) - 15, ymax = max(day1$Y) + 15)
lake3 = ggplot() + 
  geom_sf(data = profile_cropped3, color ="black", alpha = .3) +
  geom_path(data = day1, aes(x = X, y = Y, group = camera_id), color = "cyan3", alpha = .4) + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
lake3

################################################################################
#ggsave(lake1, file = "1_icefishing_social_foraging/output/example_lake1.png", 
#       width = 5, height = 5,bg = "white", dpi = 1200, units = "cm")
#ggsave(lake2, file = "1_icefishing_social_foraging/output/example_lake2.png", 
#       width = 5, height = 5,bg = "white", dpi = 1200, units = "cm")
#ggsave(lake10, file = "1_icefishing_social_foraging/output/example_lake3.png", 
#       width = 5, height = 5,bg = "white", dpi = 1200, units = "cm")
################################################################################
# make fake variation plots
colors_lakes = rev(c("cyan3", "green4", "cadetblue"))
colors_people = rev(c("red2", "burlywood", "orange3", "orange1", "maroon"))

lake_distribution1 = rnorm(50, 10, 10)
lake_distribution2 = rnorm(50, 5, 10)
lake_distribution3 = rnorm(50, 15, 10)

id_distribution1 = rnorm(10, 3, 3)
id_distribution2 = rnorm(10, 7, 4)
id_distribution3 = rnorm(10, 14, 2)
id_distribution4 = rnorm(10, 1, 4)
id_distribution5 = rnorm(10, 16, 5)


lake_data_frame = data.frame(lake1 = lake_distribution1,
                             lake2 = lake_distribution2,
                             lake3 = lake_distribution3) %>% 
  pivot_longer(1:3, names_to = "lake", values_to = "dist")
lake_variation = lake_data_frame %>% 
  ggplot(aes(x = dist, y = lake, group = lake, fill = lake, color = lake)) + 
  geom_dots(size = 1) + 
  scale_fill_manual(values = colors_lakes) + 
  scale_color_manual(values = colors_lakes) + 
  xlab("Behaviour") + 
  ylab("Lake")+ 
  labs(tag = "")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
#ggsave(lake_variation, file = "1_icefishing_social_foraging/output/lake_variation.png", 
#       width = 5, height = 5,bg = "white", dpi = 1200, units = "cm")

id_data_frame = data.frame(id1 = id_distribution1,
                             id2 = id_distribution2,
                             id3 = id_distribution3,
                           id4 = id_distribution4,
                           id5 = id_distribution5) %>% 
  pivot_longer(1:5, names_to = "id", values_to = "dist")
id_variation = id_data_frame %>% 
  ggplot(aes(x = dist, y = id, group = id, fill = id, color = id)) + 
  geom_dots(size = 1) + 
  scale_fill_manual(values = colors_people) +
  scale_color_manual(values = colors_people) +
  xlab("Behaviour") + 
  ylab("Participant ID")+ 
  labs(tag = "g")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill=NA),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
id_variation
#ggsave(id_variation, file = "1_icefishing_social_foraging/output/id_variation.png", 
#       width = 5, height = 5,bg = "white", dpi = 1200, units = "cm")
################################################################################
#illustration = readPNG("utils/data/raw_data/plot_images/variation.png")
#illustration_plot <- ggplot(mapping = aes(x = 1:max(dim(illustration)), y = 1:max(dim(illustration)))) +
#  annotation_raster(illustration, xmin = 1, xmax = dim(illustration)[2], ymin = 1, ymax = dim(illustration)[1]) +
#  geom_point(alpha = 0) + 
#  coord_fixed(ratio = 1, xlim = c(1, dim(illustration)[2]),ylim = c(1, dim(illustration)[1]), expand = 0) + 
#  labs(tag = "d") + 
#  theme(plot.tag = element_text(size = 7, face = "bold"),
#        plot.margin = margin(t = 0,  # Top margin
#                             r = 0,  # Right margin
#                             b = 0,  # Bottom margin
#                             l = 0),
#        axis.title.x = element_blank(),
#        axis.title.y = element_blank(),
#        axis.text.x = element_blank(),
#        axis.text.y = element_blank(),
#        axis.ticks.x =element_blank(),
#        axis.ticks.y =element_blank())
#illustration_plot

method = readPNG("utils/data/raw_data/plot_images/method.png")
method_plot <- ggplot(mapping = aes(x = 1:max(dim(method)), y = 1:max(dim(method)))) +
  annotation_raster(method, xmin = 1, xmax = dim(method)[2], ymin = 1, ymax = dim(method)[1]) +
  geom_point(alpha = 0) + 
  coord_fixed(ratio = 1, xlim = c(1, dim(method)[2]),ylim = c(1, dim(method)[1]), expand = 0) + 
  labs(tag = "d") + 
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
method_plot

#panel_figure1_v2 = arrangeGrob(
#  grobs = list(left_im,
#               middle_im,
#               right_im,
#               trajectories,
#               example_trajectory,
#               angling_plot,
#               illustration_plot),#,
#  #model_im),
#  widths = c(1, 1, 1, 1, 1, 1, 1, 1),
#  layout_matrix = rbind(c(1, 1, 1, 2, 2, 2, 3, 3, 3),
#                        c(1, 1, 1, 2, 2, 2, 3, 3, 3),
#                        c(1, 1, 1, 2, 2, 2, 3, 3, 3),
#                        c(4, 4, 4, 4, 5, 5, 5, 7, 7),
#                        c(4, 4, 4, 4, 5, 5, 5, 7, 7),
#                        c(4, 4, 4, 4, 5, 5, 5, 7, 7),
#                        c(4, 4, 4, 4, 5, 5, 5, 7, 7),
#                        c(4, 4, 4, 4, 6, 6, 6, 7, 7),
#                        c(4, 4, 4, 4, 6, 6, 6, 7, 7),
#                        c(4, 4, 4, 4, 6, 6, 6, 7, 7),
#                        c(4, 4, 4, 4, 6, 6, 6, 7, 7))
#)

panel_figure1 = arrangeGrob(
  grobs = list(left_im,
               middle_im,
               right_im,
               method_plot,
               example_trajectory,
               angling_plot,
               id_variation,
               lake_variation),#,
  #model_im),
  widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
  layout_matrix = rbind(c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3),
                        c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3),
                        c(1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3),
                        c(4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5),
                        c(4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5),
                        c(4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5),
                        c(4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5),
                        c(4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5),
                        c(4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6),
                        c(4, 4, 4, 4, 4, 4, 6, 6, 6, 6, 6, 6),
                        c(7, 7, 7, 8, 8, 8, 6, 6, 6, 6, 6, 6),
                        c(7, 7, 7, 8, 8, 8, 6, 6, 6, 6, 6, 6))
)


# save figure in different file formats
ggsave(panel_figure1, file = "1_icefishing_social_foraging/output/panel_figure1.png", 
       width = 18, height = 21,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_figure1, file = "1_icefishing_social_foraging/output/panel_figure1.svg", 
       width = 18, height = 21,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_figure1, file = "1_icefishing_social_foraging/output/panel_figure1.pdf", 
       width = 18, height = 21,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_figure1, file = "1_icefishing_social_foraging/output/panel_figure1.tiff",
       width = 18, height = 21,bg = "white", dpi = 1200, units = "cm")

################################################################################
# Survey analysis
################################################################################
source("1_icefishing_social_foraging/scripts/survey_analysis.r")
visualize_survey_data(survey_data, catch_data)

################################################################################
# GPS accuracy analysis
################################################################################
source("1_icefishing_social_foraging/scripts/survey_analysis.r")
visualize_survey_data(survey_data, catch_data)

################################################################################
# Video reliability analysis
################################################################################
source("1_icefishing_social_foraging/scripts/survey_analysis.r")
visualize_survey_data(survey_data, catch_data)


################################################################################
# END
################################################################################
}
