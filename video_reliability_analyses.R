################################################################################
#
# Title: Video reliability analyses
#
# Description: Assess reliability of video coding
#
# Authors: Schakowski, A.
#
# Last updated: 03/07/2024 (DD/MM/YYYY)
#
################################################################################
setwd("Projects/IceFishing/collective_foraging_dynamics_git/rproject/")
# library
source("functions/library.R")

# load data rater 1
rater_1 <- read_delim("0_rawdata/reliability_data/video_coding_data/reliability_between_rater_1.csv", delim = ";", escape_double = FALSE, trim_ws = TRUE)

# split id
rater_1$ID = str_split(rater_1$ID, pattern = "_REL", simplify = T)[,1]
rater_1$condition = "REL"

# load data rater 2
rater_2 <- read_csv("2_processed_data/video_data/video_data.csv")
rater_2$condition = "DATA"

# match datasets
rater_2$ID = paste("D", rater_2$day, "_", rater_2$year, "_C", rater_2$camera_id, sep = "")
rater_2 = rater_2 %>% 
  filter(ID %in% unique(rater_1$ID))

# keep only same columns
rater_2 = rater_2 %>% dplyr::select(ID, angling, drilling, end_of_competition, handling_time, missing_data, non_perch, other, relocating, size1, size2, size3, catch, condition)
rater_1 = rater_1 %>% dplyr::select(ID, angling, drilling, end_of_competition, handling_time, missing_data, non_perch, other, relocating, size1, size2, size3, condition)
rater_1$catch = NA
rater_1 = rater_1 %>% 
  filter(ID %in% unique(rater_2$ID))

# merge datasets
raters = rbind(rater_1, rater_2)

# number of labelled catch
n_catch_id = raters %>% 
  group_by(ID, condition) %>% 
  summarize(s = sum(size1) + sum(size2) + sum(size3) + sum(catch, na.rm = T)) %>% 
  pivot_wider(names_from = condition, values_from = s) %>% 
  ggplot(aes(x = DATA, y = REL)) + 
  geom_point() + 
  ggtitle("N Catch/ID") +
  labs(tag = "a") +
  xlab("Rater 1") +
  ylab("Rater 2") +
  geom_abline(slope = 1) + 
  theme_minimal()

# time of labelled angling sequences
raters$angling = round(raters$angling)

raters$unique_id = paste(raters$ID, raters$condition)
raters = setDT(raters)[, sequence := rleid(angling), by = "unique_id"]
raters = as.data.frame(raters)

# number of labelled angling sequences
n_labelled_sequences = raters %>% 
  filter(angling == 1) %>% 
  group_by(ID, condition) %>% 
  summarize(n = length(unique(sequence))) %>% 
  pivot_wider(names_from = condition, values_from = n) %>% 
  ggplot(aes(x = DATA, y = REL)) + 
  geom_point() + 
  ggtitle("N Spots/ID") +
  labs(tag = "b") +
  xlab("Rater 1") +
  ylab("Rater 2") +
  geom_abline(slope = 1) + 
  theme_minimal()

seq = raters %>% 
  filter(angling == 1) %>% 
  group_by(ID, condition) %>% 
  summarize(n = length(unique(sequence))) %>% 
  pivot_wider(names_from = condition, values_from = n) %>% 
  filter(abs(REL-DATA)>0)

angling_time = raters %>% 
  filter(ID != seq$ID) %>% 
  filter(angling == 1) %>%
  group_by(ID, sequence, condition) %>% 
  summarize(time = n()) %>%
  group_by(ID, condition) %>% 
  mutate(sequence = 1:n()) %>% 
  pivot_wider(names_from = condition, values_from = time) %>% 
  ggplot(aes(x = DATA, y = REL)) + 
  geom_point() + 
  ggtitle("Angling Time [s]") +
  labs(tag = "c") +
  xlab("Rater 1") +
  ylab("Rater 2") +
  geom_abline(slope = 1) + 
  theme_minimal()

# time of fish catch
spike_synchronicity = raters %>%
  filter(ID != seq$ID) %>% 
  group_by(ID, sequence,condition) %>% 
  mutate(time = 1:n()) %>% 
  mutate(catch = ifelse(is.na(catch), 0, catch)) %>% 
  mutate(catch_new = ifelse(catch == 1 | size1 == 1 | size2 == 1 | size3 == 1, 1, 0)) %>%
  mutate(total_catch_new = sum(catch_new)) %>% 
  mutate(catch_new = ifelse(condition == "REL", -catch_new, catch_new)) %>% 
  filter(total_catch_new > 0) %>% 
  ggplot(aes(x = time, y = catch_new, color = condition)) + 
  geom_path(alpha = .3) + 
  facet_wrap(~ID*sequence, scales="free") + 
  theme_minimal()+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) + 
  ggtitle("Catches") +
  labs(tag = "d")

# take time index for all labelled catches in data and rel
test = raters %>%
  filter(ID != seq$ID) %>% 
  group_by(ID,condition) %>% 
  mutate(time = 1:n()) %>% 
  mutate(catch = ifelse(is.na(catch), 0, catch)) %>% 
  mutate(catch_new = ifelse(catch == 1 | size1 == 1 | size2 == 1 | size3 == 1, 1, 0)) %>%
  mutate(total_catch_new = sum(catch_new)) %>% 
  mutate(catch_new = ifelse(condition == "REL", -catch_new, catch_new)) %>% 
  filter(total_catch_new > 0) %>% 
  filter(catch_new != 0) %>% 
  group_by(condition) %>% 
  arrange(ID, sequence) %>% 
  group_by(ID, sequence, condition) %>% 
  mutate(n_catch = 1:n())

data = test[test$condition == "DATA",]
rel = test[test$condition == "REL",]

plot(rel$time, data$time[1:902])

# make reliability panel
reliability_suppl = arrangeGrob(
  grobs = list(n_catch_id,
               n_labelled_sequences,
               angling_time,
               spike_synchronicity),
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3),
                        c(4, 4, 4),
                        c(4, 4, 4),
                        c(4, 4, 4))
)
#ggsave(panel_simulations, file = "5_visualizations/SpatialModel/panel_simulations.tiff", 
#       width = 22, height = 12, bg = "white", dpi = 300)
ggsave(reliability_suppl, file = "5_visualizations/reliability_suppl.png", 
       width = 21, height = 30, bg = "white", dpi = 300)


