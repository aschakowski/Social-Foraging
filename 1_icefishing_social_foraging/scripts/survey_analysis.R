################################################################################
# visualize survey data
################################################################################
visualize_survey_data <- function(survey_data, catch_data){
  
# filter participants of correct competition
p_ids = unique(catch_data$participant_id)
survey_data = survey_data %>% filter(participant_id %in% p_ids)

# how many fishing competitions in a year? (P2_Q4)
# older participants more competitions
plot(survey_data$P2_Q4, 2022-survey_data$P2_Q1)
plot(survey_data$P2_Q4, 2022-survey_data$P2_Q3b)

# how many days ice-fishing? (P2_)
plot(survey_data$P2_Q5, 2022-survey_data$P2_Q4)


# self rated skill and experience?
plot(survey_data$angling_skill, survey_data$P2_Q1)

# motivation
hist(survey_data$P4_Q9a_1)
hist(survey_data$P4_Q9a_2)
hist(survey_data$P4_Q9a_3)
hist(survey_data$P4_Q9a_4)
hist(survey_data$P4_Q9a_5)
hist(survey_data$P4_Q9a_6)
hist(survey_data$P4_Q9a_7)
hist(survey_data$P4_Q9a_8)
hist(survey_data$P4_Q9a_9)
hist(survey_data$P4_Q9a_10)
hist(survey_data$P4_Q9a_11)
hist(survey_data$P4_Q9a_12)
hist(survey_data$P4_Q9a_13)
hist(survey_data$P4_Q9a_14)
hist(survey_data$P4_Q9a_15)
hist(survey_data$P4_Q9a_16)
hist(survey_data$P4_Q9a_17)
hist(survey_data$P4_Q9a_18)

# importance aspects participation
survey_data$sex = ifelse(survey_data$P2_Q2 == 2, "f", "m")
motives_participation = survey_data %>%
  group_by(participant_id) %>% 
  slice(1) %>% 
  pivot_longer(P4_Q9a_1:P4_Q9a_18) %>% 
  mutate(name = plyr::revalue(name, c("P4_Q9a_1" = "others think highly of you", 
                                      "P4_Q9a_2" = "show others",
                                      "P4_Q9a_3" = "make good impression",
                                      "P4_Q9a_4" = "become better",
                                      "P4_Q9a_5" = "develop skills",
                                      "P4_Q9a_6" = "test abilities",
                                      "P4_Q9a_7" = "learn what capable",
                                      "P4_Q9a_8" = "feel independence",
                                      "P4_Q9a_9" = "be alone",
                                      "P4_Q9a_10" = "be own boss",
                                      "P4_Q9a_11" = "be free to make own choices",
                                      "P4_Q9a_12" = "be obligated to no one",
                                      "P4_Q9a_13" = "do things own way",
                                      "P4_Q9a_14" = "think for myself",
                                      "P4_Q9a_15" = "be with others who enjoy same thing",
                                      "P4_Q9a_16" = "be with people having similar values",
                                      "P4_Q9a_17" = "observe others closely",
                                      "P4_Q9a_18" = "observe others from distance"))) %>% 
  ggplot(aes(x = fct_reorder(name, value, .fun = "mean"), y = value, color = as.factor(sex), group = as.factor(sex))) + 
  scale_y_continuous(labels=c("Not at all\nimportant",
                              "Slightly\nimportant",
                              "Moderately\nimportant",
                              "Very\nimportant",
                              "Extremely\nimportant"), breaks=1:5, limits=c(1,5)) +
  stat_summary(position = position_dodge(.4), size = .6) +
  labs(color = "Gender", tag = "a")+
  scale_color_brewer(palette = "Set1")+
  geom_hline(yintercept = 3, linetype = "dotted")+
  coord_flip(ylim = c(1, 5)) +
  theme_minimal() + 
  ggtitle("Motives for Participation") + 
  ylab("Avg. Response") + 
  xlab("Item") + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x =  element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 7))
motives_participation

survey_data$AS_social_recognition
survey_data$AS_skill_development
survey_data$AL_independence
survey_data$AL_autonomy
survey_data$SP_similar_people
survey_data$NP_observing

motives_participation_grouped = survey_data %>%
  group_by(participant_id) %>% 
  slice(1) %>%  
  pivot_longer(c(AS_social_recognition,
                 AS_skill_development,
                 AL_independence,
                 AL_autonomy,
                 SP_similar_people,
                 NP_observing)) %>% 
  mutate(name = plyr::revalue(name, c("AS_social_recognition" = "Social Recognition",
                                      "AS_skill_development" = "Skill Development",
                                      "AL_independence" = "Independence",
                                      "AL_autonomy" = "Autonomy",
                                      "SP_similar_people" = "Being with Similar People",
                                      "NP_observing" = "Observing other People"))) %>% 
  ggplot(aes(x = fct_reorder(name, value, .fun = "mean"), y = value, color = as.factor(sex), group = as.factor(sex))) + 
  #ggdist::stat_halfeye(
  #  adjust = .5, 
  #  width = .6, 
  #  .width = 0, 
  #  justification = -.3, 
  #  point_colour = NA) + 
  scale_y_continuous(labels=c("Not at all\nimportant",
                              "Slightly\nimportant",
                              "Moderately\nimportant",
                              "Very\nimportant",
                              "Extremely\nimportant"), breaks=1:5, limits=c(1,5)) +
  stat_summary(position = position_dodge(.4), size = .6) +
  geom_hline(yintercept = 3, linetype = "dotted")+
  #geom_point(
  #  size = 1.5,
  #  alpha = .2,
  #  position = position_jitter(
  #    seed = 1, width = .1
  #  )
  #) + 
  coord_flip(ylim = c(1, 5)) +
  labs(color = "Gender", tag = "b")+
  theme_minimal() + 
  scale_color_brewer(palette = "Set1")+
  ggtitle("Motive Dimensions") + 
  ylab("Avg. Response") + 
  xlab("Dimension") + 
  theme(legend.position = c(.9, .2),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 7),
        axis.text.x =  element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 7))
motives_participation_grouped
motives_participation

# leisure motivation
leisure_motivation = survey_data %>%
  group_by(participant_id) %>% 
  slice(1) %>%  
  pivot_longer(P5_Q9b_1:P5_Q9b_16) %>% 
  mutate(name = plyr::revalue(name, c("P5_Q9b_1" = "develop knowledge", 
                                      "P5_Q9b_2" = "learn about things",
                                      "P5_Q9b_3" = "be close to nature",
                                      "P5_Q9b_4" = "enjoy smells and sounds of nature",
                                      "P5_Q9b_5" = "be where things are natural",
                                      "P5_Q9b_6" = "catch trophy fish",
                                      "P5_Q9b_7" = "master challenges",
                                      "P5_Q9b_8" = "catch as much as possible",
                                      "P5_Q9b_9" = "experience competition",
                                      "P5_Q9b_10" = "outwit difficult to catch fish",
                                      "P5_Q9b_11" = "win a prize",
                                      "P5_Q9b_12" = "catch fish for dinner",
                                      "P5_Q9b_13" = "get exercise",
                                      "P5_Q9b_14" = "get fresh air",
                                      "P5_Q9b_15" = "keep physically fit",
                                      "P5_Q9b_16" = "feel good after being active"))) %>% 
  ggplot(aes(x = fct_reorder(name, value, .fun = "mean"), y = value, group = as.factor(sex), color = as.factor(sex))) + 
  scale_y_continuous(labels=c("Not at all\nimportant",
                              "Slightly\nimportant",
                              "Moderately\nimportant",
                              "Very\nimportant",
                              "Extremely\nimportant"), breaks=1:5, limits=c(1,5)) +
  #ggdist::stat_halfeye(
  #  adjust = .5, 
  #  width = .6, 
  #  .width = 0, 
  #  justification = -.3, 
  #  point_colour = NA) + 
  stat_summary(position = position_dodge(.4), size = .6) +
  labs(color = "Gender", tag = "a")+
  scale_color_brewer(palette = "Set1")+
  geom_hline(yintercept = 3, linetype = "dotted")+
  #geom_point(
  #  size = 1.5,
  #  alpha = .2,
  #  position = position_jitter(
  #    seed = 1, width = .1
  #  )
  #) + 
  coord_flip(ylim = c(1, 5)) +
  theme_minimal() + 
  ggtitle("Leisure Motivation Scale") + 
  ylab("Avg. Response") + 
  xlab("Item") + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x =  element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 7))
leisure_motivation


survey_data$L_exploration
survey_data$EN_general_experience
survey_data$physical_fitness
survey_data$trophy
survey_data$challenge
survey_data$catching
leisure_motivation_grouped = survey_data %>%
  group_by(participant_id) %>% 
  slice(1) %>% 
  pivot_longer(c(L_exploration,
                 EN_general_experience,
                 physical_fitness,
                 trophy,
                 challenge,
                 catching)) %>% 
  mutate(name = plyr::revalue(name, c("L_exploration" = "Exploration",
                                      "EN_general_experience" = "Nature Experience",
                                      "physical_fitness" = "Fitness",
                                      "trophy" = "Trophy",
                                      "challenge" = "Challenge",
                                      "catching" = "Catch/Consumption"))) %>% 
  ggplot(aes(x = fct_reorder(name, value, .fun = "mean"), y = value, color = as.factor(sex), group = as.factor(sex))) + 
  scale_y_continuous(labels=c("Not at all\nimportant",
                              "Slightly\nimportant",
                              "Moderately\nimportant",
                              "Very\nimportant",
                              "Extremely\nimportant"), breaks=1:5, limits=c(1,5)) +
  #ggdist::stat_halfeye(
  #  adjust = .5, 
  #  width = .6, 
  #  .width = 0, 
  #  justification = -.3, 
  #  point_colour = NA) + 
  stat_summary(position = position_dodge(.4), size = .6) +
  scale_color_brewer(palette = "Set1")+
  geom_hline(yintercept = 3, linetype = "dotted")+
  #geom_point(
  #  size = 1.5,
  #  alpha = .2,
  #  position = position_jitter(
  #    seed = 1, width = .1
  #  )
  #) + 
  coord_flip(ylim = c(1, 5)) +
  theme_minimal() + 
  ggtitle("Leisure Motivation Dimensions") + 
  labs(color = "Gender", tag = "b")+
  ylab("Avg. Response") + 
  xlab("Dimension") + 
  theme(legend.position = c(.9, .2),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x =  element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 7))
leisure_motivation_grouped

# involvement
involvement = survey_data %>%
  group_by(participant_id) %>% 
  slice(1) %>%  
  pivot_longer(P6_Q10_1:P6_Q10_15) %>% 
  mutate(name = plyr::revalue(name, c("P6_Q10_1" = "one of the most\nenjoyable things", 
                                      "P6_Q10_2" = "very important",
                                      "P6_Q10_3" = "one of the most\nsatisfying things",
                                      "P6_Q10_4" = "life is organized\naround competitions",
                                      "P6_Q10_5" = "competitions occupy central role",
                                      "P6_Q10_6" = "changing preferences to other\nactivity would require rethinking",
                                      "P6_Q10_7" = "enjoy discussing\ncompetitions with friends",
                                      "P6_Q10_8" = "most friends connected\nto competitions",
                                      "P6_Q10_9" = "competitions provide opportunity\nto connect with friends",
                                      "P6_Q10_10" = "can really be myself\nduring competitions",
                                      "P6_Q10_11" = "identify with people\nand image of competitions",
                                      "P6_Q10_12" = "don't have to be concerned\nwith the way I look",
                                      "P6_Q10_13" = "can tell a lot about a person\nfrom seeing them in a competition",
                                    "P6_Q10_14" = "participating says a\nlot about me personally",
                                      "P6_Q10_15" = "during a competition others\nsee me the way I want to be seen"))) %>% 
  ggplot(aes(x = fct_reorder(name, value, .fun = "mean"), y = value, group = as.factor(sex), color = as.factor(sex))) + 
  scale_y_continuous(labels=c("Strongly\ndisagree",
                              "Disagree",
                              "Neither agree\nor disagree",
                              "Agree",
                              "Strongly\nagree"), breaks=1:5, limits=c(1,5)) +
  #ggdist::stat_halfeye(
  #  adjust = .5, 
  #  width = .6, 
  #  .width = 0, 
  #  justification = -.3, 
  #  point_colour = NA) + 
  stat_summary(position = position_dodge(.4), size = .6) +
  labs(color = "Gender", tag = "a") +
  scale_color_brewer(palette = "Set1")+
  geom_hline(yintercept = 3, linetype = "dotted")+
  #geom_point(
  #  size = 1.5,
  #  alpha = .2,
  #  position = position_jitter(
  #    seed = 1, width = .1
  #  )
  #) + 
  coord_flip(ylim = c(1, 5)) +
  theme_minimal() + 
  ggtitle("Involvement Scale") + 
  ylab("Avg. Response") + 
  xlab("Item") + 
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x =  element_text(size = 6),
        axis.text.y = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 7))
involvement

involvement_grouped = survey_data %>%
  group_by(participant_id) %>% 
  slice(1) %>%  
  pivot_longer(c(attraction,
                 centrality,
                 bonding,
                 affirmation,
                 expression)) %>% 
  mutate(name = plyr::revalue(name, c("attraction" = "Attraction",
                                      "centrality" = "Centrality",
                                      "bonding" = "Social Bonding",
                                      "affirmation" = "Identity Affirmation",
                                      "expression" = "Identity Expression"))) %>% 
  ggplot(aes(x = fct_reorder(name, value, .fun = "mean"), y = value, color = as.factor(sex), group = as.factor(sex))) + 
  scale_y_continuous(labels=c("Strongly\ndisagree",
                              "Disagree",
                              "Neither agree\nor disagree",
                              "Agree",
                              "Strongly\nagree"), breaks=1:5, limits=c(1,5)) +
  #ggdist::stat_halfeye(
  #  adjust = .5, 
  #  width = .6, 
  #  .width = 0, 
  #  justification = -.3, 
  #  point_colour = NA) + 
  labs(color = "Gender", tag = "b")+
  stat_summary(position = position_dodge(.4), size = .6) +
  geom_hline(yintercept = 3, linetype = "dotted")+
  #geom_point(
  #  size = 1.5,
  #  alpha = .2,
  #  position = position_jitter(
  #    seed = 1, width = .1
  #  )
  #) + 
  coord_flip(ylim = c(1, 5)) +
  theme_minimal() + 
  ggtitle("Involvement Scale Dimensions") + 
  ylab("Avg. response") + 
  xlab("Dimension") + 
  scale_color_brewer(palette = "Set1")+
  theme(legend.position = c(.9, .2),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_blank(),
        axis.text.x =  element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.tag = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 7))
involvement_grouped

# supplementary panel
panel_motives1 = arrangeGrob(
  grobs = list(motives_participation,
               motives_participation_grouped),
  widths = c(.065, 1),
  heights = c(1, .5),
  layout_matrix = rbind(c(1, 1),
                        c(NA, 2))
)

panel_motives2 = arrangeGrob(
  grobs = list(leisure_motivation,
               leisure_motivation_grouped),
  widths = c(.09, 1),
  heights = c(1, .5),
  layout_matrix = rbind(c(1, 1),
                        c(NA, 2))
)

panel_motives3 = arrangeGrob(
  grobs = list(involvement,
               involvement_grouped),
  widths = c(.05, 1),
  heights = c(1, .5),
  layout_matrix = rbind(c(1, 1),
                        c(NA, 2))
)

ggsave(panel_motives1, file = "1_icefishing_social_foraging/output/panel_motives1.png", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_motives1, file = "1_icefishing_social_foraging/output/panel_motives1.svg", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_motives1, file = "1_icefishing_social_foraging/output/panel_motives1.pdf", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_motives1, file = "1_icefishing_social_foraging/output/panel_motives1.tiff", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")

ggsave(panel_motives2, file = "1_icefishing_social_foraging/output/panel_motives2.png", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_motives2, file = "1_icefishing_social_foraging/output/panel_motives2.svg", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_motives2, file = "1_icefishing_social_foraging/output/panel_motives2.pdf", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_motives2, file = "1_icefishing_social_foraging/output/panel_motives2.tiff", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")

ggsave(panel_motives3, file = "1_icefishing_social_foraging/output/panel_motives3.png", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_motives3, file = "1_icefishing_social_foraging/output/panel_motives3.svg", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_motives3, file = "1_icefishing_social_foraging/output/panel_motives3.pdf", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")
ggsave(panel_motives3, file = "1_icefishing_social_foraging/output/panel_motives3.tiff", 
       width = 18, height = 17,bg = "white", dpi = 1200, units = "cm")



}
################################################################################
# END
################################################################################