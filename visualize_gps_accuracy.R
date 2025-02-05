################################################################################
#
# Title: Analyse Accuracy
#
# Description: This script automates the GPS accuracy analysis for fieldwork.
# 
# Last updated: 11/11/2024
#
################################################################################
# set working directory
# on linux server
if (Sys.info()[1] == "Linux"){
  print("Script needs to run on windows!")
  
} else {
  # on windows server
  setwd("P:/IceFishing/collective_foraging_dynamics_git/rproject/")
}

# check
getwd()

# load library and functions
source("functions/library.R")

################################################################################
# helper function:
extract_fit_data <- function(df, # output file of .fit to .csv converter
                             variable #variable of interest
){
  
  # get all rows with data on variable
  positions = which(df == variable, arr.ind=TRUE)
  
  # get all rows that contain timestamps
  positions_of_timestamps = which(df == "timestamp", arr.ind=TRUE)
  
  # now select entries to the right of time and variable entries and keep row index
  positions[,2] = positions[,2] + 1
  positions_of_timestamps[,2] = positions_of_timestamps[,2] + 1
  
  variable_data = data.frame(row_index = vector(), variable_data = vector())
  for (i in 1:length(positions[,1])){
    variable_data[i,] = c(positions[i,1],as.numeric(df[positions[i,1], positions[i,2]]))
  }
  
  time_data = data.frame(row_index = vector(), time_data = vector())
  for (i in 1:length(positions_of_timestamps[,1])){
    time_data[i,] = c(positions_of_timestamps[i,1],as.numeric(df[positions_of_timestamps[i,1], positions_of_timestamps[i,2]]))
  }
  
  matched_data = full_join(variable_data, time_data, by = "row_index")
  
  matched_data[order(matched_data$row_index),]
  
  colnames(matched_data) <- c("row_index",paste(variable), "timestamp")
  return(matched_data)
}


# First, read all the GPS data in. For that first define where to look. The following code works
# if you open this script from the R Project.
wd_path = getwd()

# Path where the Garmin extraction tool is saved
garmin_path = "functions/garmin_tool"

# Path where the GPS data is saved (defined relative to WD)
data_path = "0_rawdata/reliability_data/gps_accuracy_data/raw"
data_path = paste(wd_path, data_path, sep = "/")

# where should converted file be saved
save_path = "0_rawdata/reliability_data/gps_accuracy_data/converted"
save_path = paste(wd_path, save_path, sep = "/")
dir.create(save_path)

# get file names
files_galileo = list.files(paste(data_path, "GALILEO", sep = "/"))
files_glonass = list.files(paste(data_path, "GLONASS", sep = "/"))

# run parser for all these files and save in converted file. First, parse the
# GLONASS files  
for (i in 1:length(files_glonass)){
  shell(paste("cd ", paste(getwd(),garmin_path, sep = "/"), "&& java -jar ./java/FitCSVTool.jar -b ", 
              paste(data_path, "GLONASS", files_glonass[i], sep = "/"),paste(" ",save_path,"/",files_glonass[i], "_glonass.csv", sep = ""), sep = "")) 
}

# Now, read in the converted files
# GLONASS files
converted_glonass_files = list.files(paste(save_path), pattern = "glonass")

glonass_converted_data = list()
for (i in 1:length(converted_glonass_files)){
  glonass_converted_data[[i]] = read.csv(paste(save_path, converted_glonass_files[i], sep = "/"))
}

# And convert these files to readable data
# Glonass files
for (i in 1:length(converted_glonass_files)){
  # Convert first column to readable format
  colnames(glonass_converted_data[[i]])[1]<-"Type"
  
  # Select relevant entries
  glonass_converted_data[[i]] = glonass_converted_data[[i]] %>% 
    filter(Type == "Data" & Message == "record") %>% 
    filter_all(any_vars(str_detect(., pattern = 'position_lat|position_long')))
  
  # Extract variable vectors
  position_long = extract_fit_data(df = glonass_converted_data[[i]], variable = "position_long")
  position_lat = extract_fit_data(df = glonass_converted_data[[i]], variable = "position_lat")
  
  # merge variables to dataframe
  glonass_converted_data[[i]] = full_join(position_long %>% dplyr::select(-timestamp), 
                                          position_lat, 
                                          by = "row_index")
  
  # arrange by timestamp
  glonass_converted_data[[i]] = glonass_converted_data[[i]] %>% arrange(timestamp)
  
  # convert to right format
  glonass_converted_data[[i]]$date = as_datetime(as.numeric(glonass_converted_data[[i]]$timestamp) + 631065600)
  glonass_converted_data[[i]]$time_since_start = as.numeric(glonass_converted_data[[i]]$timestamp) - min(as.numeric(glonass_converted_data[[i]]$timestamp))
  glonass_converted_data[[i]]$position_long = as.numeric(glonass_converted_data[[i]]$position_long) * (180 / 2^31)
  glonass_converted_data[[i]]$position_lat = as.numeric(glonass_converted_data[[i]]$position_lat) * (180 / 2^31)
  
  # add watch ID
  glonass_converted_data[[i]]$watch_id = i
  
  # add setting
  glonass_converted_data[[i]]$setting = "glonass" 
}

data = do.call(rbind.data.frame, glonass_converted_data)

# select same time frame
data = data %>% filter(date > as.POSIXct("2022-02-17 07:45:00", tz="UTC") & date < as.POSIXct("2022-02-17 09:30:00", tz="UTC"))

# now compute pairwise distances
watch_1_glonass = data %>% arrange(date) %>% filter(setting == "glonass" & watch_id == 1) 
watch_2_glonass = data %>% arrange(date) %>% filter(setting == "glonass" & watch_id == 2) 
watch_3_glonass = data %>% arrange(date) %>% filter(setting == "glonass" & watch_id == 3) 
watch_4_glonass = data %>% arrange(date) %>% filter(setting == "glonass" & watch_id == 4) 

# Step 1: Plot all trajectories at once
data %>% 
  ggplot(aes(x = position_long, y = position_lat)) + 
  geom_point()+
  xlab("Longitude")+ylab("Latitude")+
  ggtitle("Accuracy Test")


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
  theme_classic()+
  scale_fill_brewer()+
  xlab("Pairwise distance [m]") + 
  labs(tag = "a")+
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .8),
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.background = element_rect(colour = "black", fill = "white"))
histogram_errors

boxplot_errors = ggplot(data = acc_data %>% filter(setting == "glonass"), aes(x = comparison, y = distance)) + 
  geom_violin(alpha = .2)+
  geom_boxplot(width=.1)+
  theme_classic() + 
  geom_hline(yintercept = 0, linetype = "dotted") +
  xlab("comparison") + 
  labs(tag = "b")+
  ylab("Pairwise distance [m]") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.background = element_rect(colour = "black", fill = "white"))
boxplot_errors

# summary for text
acc_data %>% group_by(setting, comparison) %>% summarize(median_distance = median(distance, na.rm = T),
                                                         min_distance = min(distance, na.rm = T),
                                                         max_distance = max(distance, na.rm = T))


# compute quantiles for text
quantile(acc_data %>% filter(setting != "galileo") %>% dplyr::select(distance), na.rm = T, p = seq(from = 0, to = 1, by = .05))

################################################################################
# compute stationary gps error
################################################################################
# load data
library(data.table)
merged_data <- fread("2_processed_data/merged_data/merged_data.csv")

# compute angling sequence id
merged_data = setDT(merged_data)[, sequence := rleid(angling), by = c("day", "year", "camera_id")]

# filter other
acc_data_angling = merged_data %>% 
  filter(angling == 1) %>% 
  dplyr::select(timestamp, 
                position_long, position_lat, 
                other,
                day, 
                year, 
                watch_id, 
                sequence)
acc_data_angling = acc_data_angling %>% filter(!other)

acc_data_angling = acc_data_angling %>% 
  ungroup() %>% 
  group_by(day, year, watch_id, sequence) %>% 
  mutate(angling_time = max(timestamp)-min(timestamp))

# compute distance to mean position for each angling spot
acc_data_angling = acc_data_angling %>% 
  group_by(day, year, watch_id, sequence) %>% 
  mutate(error = distHaversine(cbind(position_long, position_lat), cbind(mean(position_long), mean(position_lat))))

# aggregate for each spot
acc_data_angling = acc_data_angling %>% 
  group_by(day, year, watch_id, sequence) %>% 
  mutate(mean_error = mean(error),
         angling_time_cont = 1:n())

################################################################################
# visualize
################################################################################
gps_error_angling_histogram = 
  ggplot(data = acc_data_angling, aes(x = error)) + 
  geom_histogram(alpha = .8, position = "dodge", color = "black")+
  theme_classic()+
  scale_fill_brewer()+
  labs(tag = "c")+
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  xlab("distance from center [m]") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .8),
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.background = element_rect(colour = "black", fill = "white"))
gps_error_angling_histogram

gps_error_angling_mean_histogram = 
  ggplot(data = acc_data_angling %>% group_by(day, year, watch_id, sequence) %>% slice(1), aes(x = mean_error)) + 
  geom_histogram(alpha = .8, position = "dodge", color = "black")+
  theme_classic()+
  scale_fill_brewer()+
  labs(tag = "d")+
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  xlab("avg. distance from center [m]") + 
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = c(.8, .8),
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.background = element_rect(colour = "black", fill = "white"))
gps_error_angling_mean_histogram

gps_error_angling_boxplot = ggplot(data = acc_data_angling, aes(x = interaction(day, year), y = error)) + 
  geom_violin()+
  geom_boxplot(width = .1) + 
  theme_classic()+
  #scale_fill_brewer()+
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  xlab("Competition [d.yyyy]") + 
  labs(tag = "e")+
  ylab("Distance from center [m]") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.background = element_rect(colour = "black", fill = "white"))
gps_error_angling_boxplot

gps_error_angling_mean_boxplot = ggplot(data = acc_data_angling %>% group_by(day, year, watch_id, sequence) %>% slice(1), aes(x = interaction(day, year), y = mean_error)) + 
  geom_violin()+
  geom_boxplot(width = .1) + 
  theme_classic()+
  #scale_y_continuous(labels = function(x) format(x, scientific = TRUE))+
  xlab("Competition [d.yyyy]") + 
  labs(tag = "f")+
  ylab("Avg. Distance from center [m]") +
  theme(panel.background = element_rect(fill = NA),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        plot.title = element_text(size = 20),
        panel.border = element_rect(colour = "black", fill=NA),
        legend.position = "none",
        strip.text.y = element_blank(), 
        strip.background = element_blank(),
        plot.tag = element_text(size = 24, face = "bold"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 18),
        legend.background = element_rect(colour = "black", fill = "white"))
gps_error_angling_mean_boxplot

# statistics for text
acc_data_angling %>% group_by(day, year) %>% summarize(median_error = median(error, na.rm = T),
                                                       min_error = min(error, na.rm = T),
                                                       max_error = max(error, na.rm = T))

acc_data_angling %>% group_by(day, year, watch_id, sequence) %>% slice(1) %>% ungroup() %>% summarize(median_error = median(mean_error, na.rm = T),
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
ggsave(gps_error_panel, file = "5_visualizations/gps_error.png", 
       width = 21, height = 30, bg = "white", dpi = 300)

