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

# load library and functions
source("utils/library/library.R")

################################################################################
# load data
################################################################################

# dataset contains basic demographics, catch success for each competition and camera/watch ids
catch_data <- read_delim("utils/data/catch_data.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

# data contains survey responses for all participants
survey_data <- read_delim("utils/data/survey_data.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

# data contains information for all foraging locations (ID equivalent to camera_id in other datasets)
spot_selection_data = fread(file = "utils/data/spot_selection_data.csv")

# data contains information for each angling sequence
patch_leaving_data = fread(file = "utils/data/patch_leaving_data.csv")

# data contains video labels and gps for those with working recording
video_data = fread(file = "utils/data/video_data.csv")

# data only contains gps coordinates for all participants (also for those with missing recordings)
gps_data = fread(file = "utils/data/gps_data.csv")

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
  xlab("d.yyyy") + 
  ylab("Participant ID") +
  #theme_classic()+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 20),
        plot.tag = element_text(size = 24, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.5, .15),
        legend.text = element_text(size = 16),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "white", fill = "white"))
participation
ggsave(participation, file = "1_icefishing_social_foraging/output/participation_plot.png", 
       width = 10, height = 8, bg = "white", dpi = 300)

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
radii = seq(from = 0, to = 500, length.out = 100)
n_iter = 50
source("utils/library/create_network_data.r")
create_network_data(gps_data, radii, n_iter)
iterated_results = readRDS(file = "utils/data/network_data.rds")

# plot results
colors <- c("Random" = "cadetblue", "Data" = "black")
networks = ggplot() +
  stat_summary(data = iterated_results %>% dplyr::select(data, iteration, r) %>% filter(iteration == 1), aes(x = r, y = data, color = "Data"), alpha = .5) + 
  stat_lineribbon(data = iterated_results %>% 
                    dplyr::select(random, iteration, r, timestamp) %>% 
                    group_by(r, iteration) %>% summarize(random = mean(random)), 
                  aes(x = r, y = random, color = "Random"),.width = .95, alpha = .4, fill = "darkgreen")+
  scale_color_manual(values = colors)+
  theme_classic()+
  xlab(expression(epsilon~(italic(m))))+
  ylab("Relative Community Size")+
  theme(axis.title.x = element_text(size = 25),
        axis.title.y = element_text(size = 25),
        axis.text.x = element_text(size = 25),
        axis.text.y = element_text(size = 25),
        plot.title = element_text(size = 25, 
                                  hjust = -0.0, 
                                  vjust=2.12, 
                                  face = "bold"),
        #aspect.ratio = 1/1.5,
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.spacing = unit(2, "lines"),
        legend.position = c(.2, .8),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.background = element_rect(color = "black", fill = NA)) + 
  geom_vline(xintercept = 210, linetype = "dotted")+
  annotate("segment", x = 80, y = .25, xend = 100, yend = .1,
           arrow = arrow(type = "closed", length = unit(0.03, "npc"))) + 
  annotate("segment", x = 430, y = .5, xend = 410, yend = .65,
           arrow = arrow(type = "closed", length = unit(0.03, "npc"))) + 
  annotate("text", x = 430, y = .45, label = "Globally Dispersed", size = 5) + 
  annotate("text", x = 80, y = .3, label = "Locally Dense", size = 5)

ggsave(networks, file = "1_icefishing_social_foraging/output/networks.png", 
       width = 8, height = 8, bg = "white", dpi = 300)

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
  ggplot(aes(y=x, x = as.factor(as.numeric(interaction(day, year))))) + 
  geom_boxplot() + 
  ylab("Catch [kg]")+
  labs(tag = "a") +
  xlab("Competition Site")+
  geom_jitter(alpha = .4) + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))

ratio = spot_selection_data %>% 
  group_by(day, year) %>% 
  summarize(mean(Return))
min(ratio$`mean(Return)`)
max(ratio$`mean(Return)`)
median(ratio$`mean(Return)`)

ratios = spot_selection_data %>% 
  group_by(day, year, ID) %>% 
  summarize(r = mean(Return)) %>% 
  ggplot(aes(y=r, x = as.factor(as.numeric(interaction(day, year))))) + 
  geom_boxplot() + 
  ylab("Patch discoveries [%]")+
  xlab("Competition Site")+
  labs(tag = "b") +
  geom_jitter(alpha = .4) + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
ratios

spot_selection_data$ID = as.character(spot_selection_data$ID)
spot_selection_data = left_join(spot_selection_data, catch_data %>% dplyr::select(ID=camera_id, day, year, catch))
test = spot_selection_data %>% 
  group_by(day, year) %>% 
  summarize(r = mean(Return),
            r_sd = sd(Return),
            c = mean(catch),
            c_sd = sd(catch),
            n = n())
cor.test(test$r, test$c)
correl_day = test %>% 
  ggplot(aes(x = r, y = c)) + 
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = c - 2*(c_sd/sqrt(n)), ymax = c + 2*(c_sd/sqrt(n)))) + 
  geom_errorbarh(aes(xmin = r - 2*(r_sd/sqrt(n)), xmax = r + 2*(r_sd/sqrt(n)))) + 
  theme_classic() +
  labs(tag = "c") +
  xlab("Avg. Patch discoveries [%]") + 
  ylab("Avg. Catch [kg]") +
  annotate("text", x = .25, y = 500, label = "r = .851 [.478, .964]", size = 6) + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white"))
correl_day


# catch in n and kg correlated
test = spot_selection_data %>% 
  filter(total_catch > 0)
cor.test(test$total_catch, test$angling_time)
catch_time = spot_selection_data %>% 
  filter(total_catch > 0) %>% 
  ggplot(aes(x = total_catch, y = angling_time / 60)) + 
  geom_point(alpha = .2)+
  labs(tag = "d") +
  xlab("Total Catch [n fish]") + 
  ylab("Angling Time [min]") + 
  geom_smooth(method = "lm", color = "black") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white")) +
  annotate("text", x = 125, y = 25, label = "r = .715 [.670, .733]", size = 6)
catch_time

# catch in kg and n fish per day
n_dat = spot_selection_data %>% 
  group_by(day, year, camera_id = ID) %>% 
  summarize(n = sum(total_catch, na.rm = T))
n_dat$camera_id = as.character(n_dat$camera_id)
n_dat = left_join(n_dat, catch_data)

cor.test(n_dat$catch, n_dat$n)
catch_n_kg = n_dat %>% 
  ggplot(aes(x = n, y = catch/1000)) + 
  geom_point(alpha = .2)+
  labs(tag = "e") +
  xlab("Total Catch [n fish]") + 
  ylab("Catch [kg]") + 
  scale_x_continuous(breaks = seq(from = 0, to = 250, by = 50))+
  scale_y_continuous(breaks = seq(from = 0, to = 9, by = 1))+
  geom_smooth(method = "lm", color = "black") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        legend.title = element_blank(),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.background = element_rect(colour = "black", fill = "white")) +
  annotate("text", x = 200, y = 1, label = "r = .743 [.696, .783]", size = 6)
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
ggsave(panel_lake_characteristics, file = "1_icefishing_social_foraging/output/panel_lake_characteristics.png", 
       width = 17, height = 17, bg = "white", dpi = 300)

# save supplementary figures
time_allocation = video_data %>% 
  group_by(day, year, camera_id) %>% 
  summarize(angling = sum(angling, na.rm = T),
            relocating = sum(relocating, na.rm = T),
            drilling = sum(drilling, na.rm = T))

time_allocation_plot = time_allocation %>%
  group_by(day,year) %>% 
  arrange(angling) %>% 
  mutate(id = 1:n()) %>% 
  pivot_longer(angling:drilling) %>% 
  ggplot(aes(fill=name, y=value, x=id)) + 
  geom_bar(position="fill", stat="identity", color = "black") + 
  facet_wrap(~day*year, nrow = 2) + 
  labs(fill = "Activity")+ 
  ylab("Time Allocation [%]")+
  xlab("ID")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
time_allocation_plot
ggsave(time_allocation_plot, file = "1_icefishing_social_foraging/output/time_allocation_plot.png", 
       width = 18, height = 8, bg = "white", dpi = 300)

spots_visited = spot_selection_data %>% 
  group_by(ID, day, year) %>% 
  summarize(n = n()) %>% 
  ggplot(aes(x = n)) + 
  geom_histogram(color = "black") + 
  ylab("ID count")+
  xlab("Visited Foraging Locations [n]")+
  facet_wrap(~day*year, nrow = 2)+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
ggsave(spots_visited, file = "1_icefishing_social_foraging/output/spots_visited.png", 
       width = 18, height = 8, bg = "white", dpi = 300)

avg_relocation = 
  spot_selection_data %>% filter(step>0) %>% 
  group_by(ID, day, year) %>% 
  summarize(step = mean(step)) %>% 
  ggplot(aes(x = step)) + 
  geom_histogram(color = "black")+
  ylab("ID count")+
  xlab("Avg. Relocation Distance [m]")+
  facet_wrap(~day*year, nrow = 2) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
ggsave(avg_relocation, file = "1_icefishing_social_foraging/output/avg_relocation.png", 
       width = 18, height = 8, bg = "white", dpi = 300)


total_relocation = 
  spot_selection_data %>% filter(step>0) %>% 
  group_by(ID, day, year) %>% 
  summarize(step = sum(step)) %>% 
  ggplot(aes(x = step/1000)) + 
  geom_histogram(color = "black")+
  ylab("ID count")+
  xlab("Total Relocation Distance [km]")+
  facet_wrap(~day*year, nrow = 2) +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 18, angle = 90, hjust = .5),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
ggsave(total_relocation, file = "1_icefishing_social_foraging/output/total_relocation.png", 
       width = 18, height = 8, bg = "white", dpi = 300)

################################################################################
# create figure 1
################################################################################

# data containing depth contours
lake_data = readRDS("utils/data/depth_profiles.rds")

# data containing boundary of competition area
competition_area = readRDS("utils/data/study_lakes.rds")

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

# project coordinates to XY
xy = LongLatToUTM(day5$position_long, day5$position_lat, zone = 35)
day5$X = xy$coords.x1
day5$Y = xy$coords.x2
day5 = st_as_sf(day5, coords = c("X", "Y"))
st_crs(day5) <- "epsg:3067"
day5$X = xy$coords.x1
day5$Y = xy$coords.x2

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
rectangle = extent(subset)

trajectories = ggplot() + 
  geom_sf(data = depth_profiles[[9]], color ="black", alpha = .3) +
  geom_path(data = day5, aes(x = X, y = Y, group = camera_id), color = "darkgreen", alpha = .4) + 
  geom_path(data = day5 %>% filter(camera_id == 5), aes(x = X, y = Y, group = camera_id), color = "black", alpha = 1, size = .8) + 
  geom_rect(aes(xmin = rectangle[1]-15, xmax = rectangle[2]+15, ymin = rectangle[3]-15, ymax = rectangle[4]+15), fill = NA, color = "black", linetype = "dashed", size = 1)+
  ylab("Northing")+
  xlab("Easting")+
  labs(tag = "d")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white"))
trajectories

# zoom into map
profile_cropped = st_crop(depth_profiles[[9]], xmin = min(subset$x) - 15, xmax = max(subset$x) + 15,
                          ymin = min(subset$y) - 15, ymax = max(subset$y) + 15)

# assign color values
colors <- c("successful" = "red2", "unsuccessful   " = "black")

# read images to display on plot
relocation <- readPNG("utils/data/plot_images/relocating.png")
angling <- readPNG("utils/data/plot_images/angling.png")

# coordinates for arrows in plot
coords_arrows = subset %>% filter(Return == 0) 
coords_arrows = coords_arrows[30, c("x", "y")] 

zoom = extent(profile_cropped)
example_trajectory = ggplot() + 
  geom_sf(data = profile_cropped, color ="black", alpha = .2) +
  geom_path(data = subset, aes(x = x, y = y), color = "black", alpha = .4) +
  geom_point(data = subset %>% filter(Return == 1), aes(x = x, y = y, color = "successful"), stroke = 2, size = 2, alpha = .5) +
  geom_point(data = subset %>% filter(Return == 0), aes(x = x, y = y, color = "unsuccessful   "), stroke = 1, size = 1, alpha = .5) +
  ylab("Northing")+
  xlab("Easting")+
  labs(tag = "e")+
  scale_color_manual(values = colors)+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white")) + 
  annotation_raster(relocation, xmin = zoom[1] + 30, xmax = zoom[1] + 130, ymin = zoom[3] + 65, ymax = zoom[3] + 65 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50)) + 
  geom_segment(aes(x = zoom[1] + 30 + 50, y = zoom[3] + 65 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50), xend = zoom[1] + 130 + 10, yend = zoom[3] + 200 + 40),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 2)+
  annotation_raster(angling, xmin = zoom[1] + 250, xmax = zoom[1] + 350, ymin = zoom[3] + 200, ymax = zoom[3] + 200 + 1*dim(relocation)[1]/dim(relocation)[2] * (150-50)) +
  geom_segment(aes(x = zoom[1] + 250 + 50, y = zoom[3] + 200, xend = coords_arrows$x+10, yend = coords_arrows$y+10),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 2)
example_trajectory

# angling plot
angling_sub = patch_leaving_data %>% 
  filter(day == 5 & year == 2022 & camera_id == 5) %>% 
  filter(sequence == 50) %>% 
  mutate(cumulative_r = cumsum(reward)) %>% 
  mutate(bin = bin - 1) %>% 
  mutate(bin = 10 * bin / 60)

# load catch/start/leave pictures
start <- readPNG("utils/data/plot_images/start.png")
catch1 <- readPNG("utils/data/plot_images/catch1.png")
catch2 <- readPNG("utils/data/plot_images/catch2.png")
leave <- readPNG("utils/data/plot_images/leave.png")

angling_plot = ggplot() + 
  geom_point(data = angling_sub %>% filter(reward == 1), aes(x = bin, y = cumulative_r)) + 
  geom_path(data = angling_sub, aes(x = bin, y = cumulative_r)) + 
  geom_vline(data = angling_sub %>% filter(reward == 1), aes(xintercept = bin), linetype = "dotted") + 
  geom_vline(data = angling_sub %>% filter(bin == 0 | bin == max(bin)), aes(xintercept = bin)) + 
  ylab("Cumulative Catch [n]") + 
  xlab("Angling Time [min]") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .2),
        legend.text = element_text(size = 18),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.title = element_blank(),
        legend.background = element_rect(colour = "black", fill = "white")) + 
  coord_fixed(ratio = .75)+
  scale_x_continuous(breaks = seq(from = 0, to = 10, by = 1))+
  scale_y_continuous(breaks = seq(from = 0, to = 10, by = 1))+
  labs(tag = "f")+
  geom_rect(aes(xmin = 0-.15, xmax = 0+.15, ymin = 7-.15*1/.75, ymax = 7+.15*1/.75), color = "black",fill = "white", size = 1)+
  geom_rect(aes(xmin = 1-.15, xmax = 1+.15, ymin = 7-.15*1/.75, ymax = 7+.15*1/.75), color = "black",fill = "white", size = 1)+
  geom_rect(aes(xmin = 1.33-.15, xmax = 1.33+.15, ymin = 7-.15*1/.75, ymax = 7+.15*1/.75), color = "black",fill = "white", size = 1)+
  geom_rect(aes(xmin = 2.17-.15, xmax = 2.17+.15, ymin = 7-.15*1/.75, ymax = 7+.15*1/.75), color = "black",fill = "white", size = 1)+
  geom_rect(aes(xmin = 3.83-.15, xmax = 3.83+.15, ymin = 7-.15*1/.75, ymax = 7+.15*1/.75), color = "black",fill = "white", size = 1)+
  geom_rect(aes(xmin = 4.67-.15, xmax = 4.67+.15, ymin = 7-.15*1/.75, ymax = 7+.15*1/.75), color = "black",fill = "white", size = 1)+
  geom_rect(aes(xmin = 7.33-.15, xmax = 7.33+.15, ymin = 7-.15*1/.75, ymax = 7+.15*1/.75), color = "black",fill = "white", size = 1)+
  geom_rect(aes(xmin = 8.5-.15, xmax = 8.5+.15, ymin = 7-.15*1/.75, ymax = 7+.15*1/.75), color = "black",fill = "white", size = 1)+
  annotate("text", x = 0, y=7, label = "S", size = 7) +
  annotate("text", x = 1, y=7, label = "C", size = 7) +
  annotate("text", x = 1.33, y=7, label = "C", size = 7) +
  annotate("text", x = 2.17, y=7, label = "C", size = 7) +
  annotate("text", x = 3.83, y=7, label = "C", size = 7) +
  annotate("text", x = 4.67, y=7, label = "C", size = 7) +
  annotate("text", x = 7.33, y=7, label = "C", size = 7) +
  annotate("text", x = 8.5, y=7, label = "L", size = 7) +
  annotation_raster(start, xmin = .1, xmax = .9, ymin = 1.8, ymax = 1.8 + 1/.75*dim(start)[1]/dim(start)[2] * (.9-.1)) + 
  geom_segment(aes(x = .1 + (.9-.1) / 2, y = 1.8, xend = 0, yend = 0),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 2)+
  annotation_raster(catch1, xmin = 2.4, xmax = 3.2, ymin = 3.8, ymax = 3.8 + 1/.75*dim(start)[1]/dim(start)[2] * (.9-.1)) + 
  geom_segment(aes(x = 2.4 + (.9-.1) / 2, y = 3.8, xend = 2.17, yend = 3),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 2)+
  annotation_raster(catch2, xmin = 5, xmax = 5.8, ymin = 3, ymax = 3 + 1/.75*dim(start)[1]/dim(start)[2] * (.9-.1)) + 
  geom_segment(aes(x = 5 + (.9-.1) / 2, y = 2.5 + 2*dim(start)[1]/dim(start)[2] * (.9-.1), xend = 4.67, yend = 5),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 2)+
  annotation_raster(leave, xmin = 7.5, xmax = 8.3, ymin = 4, ymax = 4 + 1/.75*dim(start)[1]/dim(start)[2] * (.9-.1)) + 
  geom_segment(aes(x = 7.5 + (.9-.1) / 2, y = 3.5 + 2*dim(start)[1]/dim(start)[2] * (.9-.1), xend = 8.5, yend = 6),
               arrow = arrow(length = unit(0.3, "cm")), color = "black", size = 2) 
angling_plot


left <- readPNG("utils/data/plot_images/left.png")
middle <- readPNG("utils/data/plot_images/middle.png")
right <- readPNG("utils/data/plot_images/right.png")

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
  theme(plot.tag = element_text(size = 24, face = "bold"),
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
middle_im <- ggplot(mapping = aes(x = 1:max(dim(middle)), y = 1:max(dim(middle)))) +
  annotation_raster(middle, xmin = 1, xmax = dim(middle)[2], ymin = 1, ymax = dim(middle)[1]) +
  geom_point(alpha = 0) + 
  coord_fixed(ratio = 1, xlim = c(1, dim(middle)[2]),ylim = c(1, dim(middle)[1]), expand = 0) + 
  labs(tag = "b") + 
  theme(plot.tag = element_text(size = 24, face = "bold"),
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
right_im <- ggplot(mapping = aes(x = 1:max(dim(right)), y = 1:max(dim(right)))) +
  annotation_raster(right, xmin = 1, xmax = dim(right)[2], ymin = 1, ymax = dim(right)[1]) +
  geom_point(alpha = 0) + 
  coord_fixed(ratio = 1, xlim = c(1, dim(right)[2]),ylim = c(1, dim(right)[1]), expand = 0) + 
  labs(tag = "c") + 
  theme(plot.tag = element_text(size = 24, face = "bold"),
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
left_im
middle_im
right_im

################################################################################
# model facet
################################################################################
model_framework = readPNG("utils/data/plot_images/model_framework.png")
model_im <- ggplot(mapping = aes(x = 1:max(dim(model_framework)), y = 1:max(dim(model_framework)))) +
  annotation_raster(model_framework, xmin = 1, xmax = dim(model_framework)[2], ymin = 1, ymax = dim(model_framework)[1]) +
  geom_point(alpha = 0) + 
  coord_fixed(ratio = 1, xlim = c(1, dim(model_framework)[2]),ylim = c(1, dim(model_framework)[1]), expand = 0) + 
  labs(tag = "g") + 
  theme(plot.tag = element_text(size = 24, face = "bold"),
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
model_im

################################################################################
# make panel
################################################################################
# make lesion panel
panel_figure1 = arrangeGrob(
  grobs = list(left_im,
               middle_im,
               right_im,
               trajectories,
               example_trajectory,
               angling_plot,
               model_im),
  widths = c(1, 1, 1, 1, 1, 1),
  layout_matrix = rbind(c(1, 1, 2, 2, 3, 3),
                        c(1, 1, 2, 2, 3, 3),
                        c(4, 4, 4, 5, 5, 5),
                        c(4, 4, 4, 5, 5, 5),
                        c(4, 4, 4, 5, 5, 5),
                        c(6, 6, 6, 7, 7, 7),
                        c(6, 6, 6, 7, 7, 7))
)
ggsave(panel_figure1, file = "1_icefishing_social_foraging/output/panel_figure1.png", 
       width = 24, height = 26, bg = "white", dpi = 300)

################################################################################
# END
################################################################################
}
