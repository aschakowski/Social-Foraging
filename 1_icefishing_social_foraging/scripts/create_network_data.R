################################################################################
#
# Title: Descriptive Network Analysis
#
# Description: Simulate random locations and evaluate level of aggregation
#
# Authors: Schakowski, A.
#
# Last updated: 06/11/2024 (DD/MM/YYYY)
#
################################################################################
create_network_data <- function(gps_positions, radii, n_iter){
  
  # downsample data
  gps_positions = gps_positions %>% filter(timestamp %% 60 == 0)
  
  # filter timestamps that do not occur often enough
  gps_positions = gps_positions %>% group_by(timestamp) %>% mutate(n = n()) %>% 
    filter(n > 30)
  
  # source helper functions
  sourceCpp("utils/library/fast_dist.cpp")
  sourceCpp("utils/library/ptinpoly.cpp")
  
  # Import Lake Area
  competition_area = readRDS("utils/data/raw_data/lake_data/study_lakes.rds")
  
  # polygon lakes for faster execution time
  point_lakes = list()
  for (i in 1:10){
    point_lakes[[i]] = competition_area[[i]] %>% 
      st_buffer(10) %>% 
      st_cast("POINT") %>% 
      as.data.frame()
    point_lakes[[i]] = point_lakes[[i]] %>% 
      mutate(x = unlist(map(point_lakes[[i]]$polygons,1)),
             y = unlist(map(point_lakes[[i]]$polygons,2))) %>% 
      dplyr::select(x,y)
  }
  
  # lake key
  lake_key = data.frame(day = c(1:5, 1:3, 5,4), 
                        year = rep(c(2023, 2022), each = 5),
                        lake_index = 1:10)
  gps_positions = left_join(gps_positions, lake_key)
  
  # compute XY from LonLat
  xy = LongLatToUTM(gps_positions$position_long, gps_positions$position_lat, zone = 35)
  gps_positions$X = xy[ ,2]
  gps_positions$Y = xy[ ,3]
  
  # check if any coordinates are outside competition area,
  # compute competition area as unit of drawn area and visited spots
  check = vector()
  c = 1
  competition_area_new=list()
  for (i in 1:10){
    
    area = point_lakes[[i]]
    
    coords = gps_positions %>% ungroup() %>% filter(day == lake_key[i, "day"] & year == lake_key[i, "year"]) %>% dplyr::select(X, Y)
    
    # convex hull of angling spots
    hull = coords[coords[, 1:2] %>% chull(),c("X", "Y")]
    colnames(hull) <- c("x", "y")
    hull2 = rbind(area, hull)
    hull3 = hull2[hull2 %>% chull(),c("x", "y")]
    
    for (j in 1:nrow(coords)){
      check[c] = pnpoly(cbind(coords$X[j], coords$Y[j]), as.matrix(hull3)) 
      c = c+1
    }
    
    competition_area_new[[i]] <- hull3
    
  }
  
  # loop through all lakes and generate n_iter random locations on the lake
  # and evaluate community size
  
  iterated_results = list()
  for (n in 1:n_iter){
  
    lake_results = list()
    for (lake in unique(gps_positions$lake_index)){
      
      # data
      data_lake = gps_positions %>% filter(lake_index == lake)
      
      # loop through all timesteps
      full_results = list()
      for (t in 1:length(unique(data_lake$timestamp))){
      
      # actual locations  
      data_locations = data_lake %>% filter(timestamp == unique(data_lake$timestamp)[t])
      
      # sim
      sim_locations = st_coordinates(st_sample(competition_area[[lake]], size = nrow(data_locations), type = "random"))
      
      # evaluate community structure for different epsilon
      # preallocate results
      avg_community_size = data.frame(data = vector(mode = "numeric", length = length(radii)),
                                      random = vector(mode = "numeric", length = length(radii)),
                                      r = radii,
                                      iteration = n,
                                      timestamp = unique(data_lake$timestamp)[t],
                                      lake = lake)
      
        for (r in 1:length(radii)){
        
          # identify clusters
          clusters = dbscan::dbscan(data_locations[, c("X", "Y")], eps = radii[r], minPts = 1)
        
          # avg cluster size as inverse of number of identified clusters
          avg_community_size[r, "data"] = 1/length(unique(clusters$cluster))
        
          # identify clusters
          clusters = dbscan::dbscan(sim_locations[, c("X", "Y")], eps = radii[r], minPts = 1)
        
          # avg cluster size as inverse of number of identified clusters
          avg_community_size[r, "random"] = 1/length(unique(clusters$cluster))
        
        }#radii
      
        full_results[[t]] <- avg_community_size
  
      } #timestamps
    
      lake_results[[lake]] = do.call(rbind, full_results)
      
    } #lakes
    
    iterated_results[[n]] <- do.call(rbind, lake_results)
    
    print(paste(n, "/", n_iter))
    
  } # iterations
  
  iterated_results = do.call(rbind, iterated_results)
  
  saveRDS(iterated_results, file = "utils/data/processed_data/network_data.rds")
  
}

