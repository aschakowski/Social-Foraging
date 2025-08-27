################################################################################
# Video reliability analyses
################################################################################
analyse_video_reliability <- function(){

# functions
library(Rcpp)
sourceCpp("utils/library/fast_dist.cpp")
  
# load data
raters = readRDS("utils/data/raw_data/reliability_data/video_coder_reliability.rds")

# number of labelled angling sequences
test = raters %>% 
  filter(angling == 1) %>% 
  group_by(ID, condition) %>% 
  dplyr::summarize(n = length(unique(sequence))) %>% 
  pivot_wider(names_from = condition, values_from = n)
cor.test(test$DATA, test$REL)

n_labelled_sequences = raters %>% 
  filter(angling == 1) %>% 
  group_by(ID, condition) %>% 
  dplyr::summarize(n = length(unique(sequence))) %>% 
  pivot_wider(names_from = condition, values_from = n) %>% 
  ggplot(aes(x = DATA, y = REL)) + 
  geom_point(size = .5) + 
  ggtitle("N = 26 IDs") +
  labs(tag = "A") +
  xlab("N spots/ID Rater 1") +
  ylab("N spots/ID Rater 2") +
  geom_abline(slope = 1) + 
  #theme_classic() + 
  #annotate(label = "r = .999 [.999, .999]", geom = "text", x = 45, y = 25, size = 5) + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "black", fill = "white"))
n_labelled_sequences

seq = raters %>% 
  filter(angling == 1) %>% 
  group_by(ID, condition) %>% 
  dplyr::summarize(n = length(unique(sequence))) %>% 
  pivot_wider(names_from = condition, values_from = n) %>% 
  filter(abs(REL-DATA)>0)

test = raters %>% 
  filter(!ID %in% seq$ID) %>% 
  filter(angling == 1) %>%
  group_by(ID, sequence, condition) %>% 
  dplyr::summarize(time = n()) %>% 
  pivot_wider(names_from = condition, values_from = time)
cor.test(test$DATA, test$REL)

raters %>% 
  filter(!ID %in% seq$ID) %>%
  filter(angling == 1) %>% 
  filter(condition == "DATA") %>% 
  mutate(x = paste(ID, sequence)) %>%
  ungroup() %>% 
  dplyr::summarize(length(unique(x)))


angling_time = raters %>% 
  filter(!ID %in% seq$ID) %>% 
  filter(angling == 1) %>%
  group_by(ID, sequence, condition) %>% 
  dplyr::summarize(time = n()) %>%
  group_by(ID, condition) %>% 
  mutate(sequence = 1:n()) %>% 
  pivot_wider(names_from = condition, values_from = time) %>% 
  ggplot(aes(x = DATA, y = REL)) + 
  geom_point(size = .5) + 
  ggtitle("N = 924 spots") +
  labs(tag = "C") +
  xlab("Angling Time (s) Rater 1") +
  ylab("Angling Time (s) Rater 2") +
  geom_abline(slope = 1) + 
  #theme_classic() + 
  #annotate(label = "r = .999 [.999, .999]", geom = "text", x = 2000, y = 250, size = 5) + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "black", fill = "white"))
angling_time

# number of labelled catch
test = raters %>% 
  filter(!ID %in% seq$ID) %>% 
  group_by(ID, condition, sequence) %>% 
  dplyr::summarize(s = sum(fish_catch)) %>% 
  pivot_wider(names_from = condition, values_from = s)
cor.test(test$DATA, test$REL)

n_catch_spot = raters %>% 
  group_by(ID, condition, sequence) %>%
  filter(angling == 1) %>% 
  filter(!ID %in% seq$ID) %>% 
  dplyr::summarize(s = sum(fish_catch)) %>% 
  pivot_wider(names_from = condition, values_from = s) %>% 
  ggplot(aes(x = DATA, y = REL)) + 
  geom_point(size = .5) + 
  ggtitle("N = 924 spots") +
  labs(tag = "B") +
  xlab("N catch Rater 1") +
  ylab("N catch Rater 2") +
  geom_abline(slope = 1) + 
  #theme_classic()+
  #annotate(label = "r = .999 [.999, .999]", geom = "text", x = 75, y = 25, size = 5) + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "black", fill = "white"))
n_catch_spot

# time of fish catch
# compute synchronicity between raters

# extract sequence
test = raters %>%
  filter(angling == 1) %>% 
  filter(!ID %in% seq$ID) %>% 
  group_by(ID, sequence,condition) %>% 
  mutate(time = 1:n()) %>% 
  #mutate(catch = ifelse(is.na(catch), 0, catch)) %>% 
  #mutate(catch_new = ifelse(catch == 1 | size1 == 1 | size2 == 1 | size3 == 1, 1, 0)) %>%
  mutate(total_catch_new = sum(fish_catch)) %>% 
  mutate(catch_new = ifelse(condition == "REL", -fish_catch, fish_catch)) %>% 
  filter(total_catch_new > 0)

# loop through all sequences
test$unique_sequence_id = paste(test$ID, test$sequence)
results_df = data.frame(sequence = vector(),
                        synch = vector(),
                        delta_catch = vector())
for (i in 1:length(unique(test$unique_sequence_id))){
  
  # seq
  seq = unique(test$unique_sequence_id)[i]
  
  # subset data
  dat = test %>% filter(unique_sequence_id == seq)
  
  # extract both sequences
  x1 = pull(dat[dat$condition == "DATA","fish_catch"])
  x2 = pull(dat[dat$condition == "REL","fish_catch"])
  
  # extract timestamp of catches
  t1 = which(x1==1)
  t2 = which(x2==1)
  
  # compute for each element shortest distance to other elements
  distances = fastPdist2(as.matrix(t1), as.matrix(t2))
  
  # which vector shorter?
  if (length(t1)<length(t2)){
    
    tmp = apply(distances, 1, min)
    tmp = tmp[tmp<20] # assuming that with 20s difference both were coding the same fish catch
    
    synch = mean(tmp)
    
  } else {
    
    tmp = apply(distances, 2, min)
    tmp = tmp[tmp<20]
    synch = mean(tmp)
    
  }
  
  # save n catch
  delta_catch = abs(sum(x1)-sum(x2))
  
  # save to df
  results_df[i, "sequence"] = seq
  results_df[i, "synch"] = synch
  results_df[i, "delta_catch"] = delta_catch
  
}

# raw values
synch = vector()
for (i in 1:length(unique(test$unique_sequence_id))){
  
  # seq
  seq = unique(test$unique_sequence_id)[i]
  
  # subset data
  dat = test %>% filter(unique_sequence_id == seq)
  
  # extract both sequences
  x1 = pull(dat[dat$condition == "DATA","fish_catch"])
  x2 = pull(dat[dat$condition == "REL","fish_catch"])
  
  # extract timestamp of catches
  t1 = which(x1==1)
  t2 = which(x2==1)
  
  # compute for each element shortest distance to other elements
  distances = fastPdist2(as.matrix(t1), as.matrix(t2))
  
  # which vector shorter?
  if (length(t1)<length(t2)){
    
    tmp = apply(distances, 1, min)
    tmp = tmp[tmp<20] # assuming that with 20s difference both were coding the same fish catch
    
    synch = c(synch,tmp)
    
  } else {
    
    tmp = apply(distances, 2, min)
    tmp = tmp[tmp<20]
    synch = c(synch,tmp)
    
  }
  
}


# plot results
catch_coding_s = results_df %>% 
  ggplot(aes(x = synch)) + 
  geom_histogram(color = "black", fill = "darkgrey") + 
  theme_classic() + 
  ggtitle("N = 162 spots") + 
  ylab("count") +
  labs(tag = "F") +
  xlab("Mean abs. time difference catches rater 1 - rater 2 [s]") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "black", fill = "white"))
catch_coding_s

mean(synch)
length(synch)
quantile(synch, prob = seq(from = 0, to = 1, by = .05))
catch_coding_s_noaggregated = data.frame(delta = synch) %>% 
  ggplot(aes(x = delta)) + 
  geom_histogram(color = "black", fill = "darkgrey") + 
  #theme_classic() + 
  ggtitle("N = 1257 catches") + 
  ylab("count") +
  labs(tag = "D") +
  xlab("Abs. time difference\nRater 1 - Rater 2 (s)") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "black", fill = "white"))
catch_coding_s_noaggregated

catch_coding_n = results_df %>% 
  ggplot(aes(x = delta_catch)) + 
  geom_histogram(color = "black", fill = "darkgrey") + 
  #theme_classic() + 
  ggtitle("N = 164 spots") + 
  ylab("count") + 
  labs(tag = "E") +
  xlab("Difference number of fish catches rater 1 - rater 2 (n)") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 7),
        axis.title.y = element_text(size = 7),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 7),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 7, face = "bold"),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        legend.background = element_rect(colour = "black", fill = "white"))
catch_coding_n

# make reliability panel
library(gridExtra)
reliability_suppl = arrangeGrob(
  grobs = list(n_labelled_sequences,
               n_catch_spot,
               angling_time,
               catch_coding_s_noaggregated),
  widths = c(1, 1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3, 4))
)
ggsave(reliability_suppl, file = "1_icefishing_social_foraging/output/panel_interrater.png", 
       width = 18, height = 5,bg = "white", dpi = 600, units = "cm")
ggsave(reliability_suppl, file = "1_icefishing_social_foraging/output/panel_interrater.pdf", 
       width = 18, height = 5,bg = "white", dpi = 600, units = "cm")
ggsave(reliability_suppl, file = "1_icefishing_social_foraging/output/panel_interrater.svg", 
       width = 18, height = 5,bg = "white", dpi = 600, units = "cm")
ggsave(reliability_suppl, file = "1_icefishing_social_foraging/output/panel_interrater.tiff", 
       width = 18, height = 5,bg = "white", dpi = 600, units = "cm")
}
