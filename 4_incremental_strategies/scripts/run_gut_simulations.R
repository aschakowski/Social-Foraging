################################################################################
#
# Title: Run giving-up time simulations
#
# Description: Run analyses and visualize results
#
# Authors: Schakowski, A.
#
# Last updated: 02/04/2024 (DD/MM/YYYY)
#
################################################################################
run_gut_simulations <- function(){

# load models
e_spatial_features = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/e_patch_leaving_model_spatial_features_fit.rds")
d_globallocal = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/d_patch_leaving_model_globallocal_updating_fit.rds")
c_patch_discovery = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/c_patch_leaving_model_patch_discovery_fit.rds")
b_time_wo_catch = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/b_patch_leaving_model_time_wo_catch_fit.rds")
a_baseline = readRDS(file = "utils/data/processed_data/model_fit/patch_leaving_models/a_patch_leaving_model_baseline_fit.rds")

# load stan data
stan_data = readRDS(file = "utils/data/processed_data/stan_data_patch_leaving.rds")
gc()

# simulation functions
sim_spatial_patch_leaving_model <- function(iter, 
                                            stan_data, 
                                            model){
  
  # extract parameters
  sim_main_effects = model %>% 
    spread_draws(beta[parameter]) %>% 
    filter(parameter < 10) %>% 
    group_by(parameter) %>% 
    summarize(beta = mean(beta))
  
  sim_main_effects_alpha = model %>% 
    spread_draws(alpha[parameter]) %>% 
    group_by(parameter) %>% 
    summarize(alpha = mean(alpha))
  
  sim_lake_effects = model %>% 
    spread_draws(v_lake[lake, parameter]) %>% 
    filter(parameter < 10) %>% 
    group_by(parameter, lake) %>% 
    summarize(lake_offset = mean(v_lake))
  
  sim_lake_effects_alpha = model %>% 
    spread_draws(v_lake[lake, parameter]) %>% 
    filter(parameter >= 10) %>% 
    mutate(parameter = parameter - 9) %>% 
    group_by(parameter, lake) %>% 
    summarize(lake_offset = mean(v_lake))
  
  sim_trip_effects = model %>% 
    spread_draws(v_trip[trip, parameter]) %>% 
    filter(parameter < 10) %>% 
    group_by(parameter, trip) %>% 
    summarize(trip_offset = mean(v_trip))
  
  sim_trip_effects_alpha = model %>% 
    spread_draws(v_trip[trip, parameter]) %>% 
    filter(parameter >= 10) %>%
    mutate(parameter = parameter - 9) %>% 
    group_by(parameter, trip) %>% 
    summarize(trip_offset = mean(v_trip))
  
  sim_id_effects = model %>% 
    spread_draws(v_id[id, parameter]) %>% 
    filter(parameter < 10) %>% 
    group_by(parameter, id) %>% 
    summarize(id_offset = mean(v_id))
  
  sim_id_effects_alpha = model %>% 
    spread_draws(v_id[id, parameter]) %>% 
    filter(parameter >= 10) %>%
    mutate(parameter = parameter - 9) %>% 
    group_by(parameter, id) %>% 
    summarize(id_offset = mean(v_id))
  
  
  
  trip_res = list()
  c = 1
  for (j in stan_data$trips){
    
    # extract trip level data
    start = stan_data$trip_start[j]
    end = stan_data$trip_end[j]
    
    # extract data
    state = stan_data$state[start:end]
    catch = stan_data$catch[start:end]
    cumulative_catch = stan_data$cumulative_catch[start:end]
    angling_time = stan_data$angling_time[start:end]
    time_since_event = stan_data$time_since_event[start:end]
    
    # spatial features
    social_feature = stan_data$social_feature[start:end]
    success_feature = stan_data$success_feature[start:end]
    loss_feature = stan_data$loss_feature[start:end]
    roughness_feature = stan_data$roughness_feature[start:end]
    
    # create spot identifier
    spot = c(1,lag(cumsum(state) + 1)[-1])
    
    # extract parameter identifiers
    lake = stan_data$lakes[start:end][1]
    trip = j
    id = stan_data$ids[start:end][1]
    
    # preallocate results
    res = list()
    
    # loop through all spots
    # counter for trip length
    trip_length_counter = 2
    
    # initialize local and global reward
    global_reward = stan_data$initial_values[j]
    
    # loop through all angling spots    
    for (s in unique(spot)){
      
      # subset predictors
      catch_tmp = catch[spot == s]
      time_since_event_tmp = time_since_event[spot == s]
      angling_time_tmp = angling_time[spot == s]
      
      # subset spatial features
      social_feature_tmp = social_feature[spot == s]
      success_feature_tmp = success_feature[spot == s]
      loss_feature_tmp = loss_feature[spot == s]
      roughness_feature_tmp = roughness_feature[spot == s]
      
      # extend time varying predictors
      catch_tmp = c(catch_tmp, rep(0, 1500 - length(catch_tmp)))
      time_since_event_tmp = c(time_since_event_tmp, seq(from = time_since_event_tmp[length(time_since_event_tmp)]+1, to = 1500, by = 1))[1:1500]
      cumulative_catch_tmp = cumsum(catch_tmp)
      
      # extend spatial features
      social_feature_tmp = c(social_feature_tmp, rep(social_feature_tmp[length(social_feature_tmp)], 1500 - length(social_feature_tmp)))
      success_feature_tmp = c(success_feature_tmp, rep(success_feature_tmp[length(success_feature_tmp)], 1500 - length(success_feature_tmp)))
      loss_feature_tmp = c(loss_feature_tmp, rep(loss_feature_tmp[length(loss_feature_tmp)], 1500 - length(loss_feature_tmp)))
      roughness_feature_tmp = c(roughness_feature_tmp, rep(roughness_feature_tmp[length(roughness_feature_tmp)], 1500 - length(roughness_feature_tmp)))
      
      # compute regression weights
      weights = sim_main_effects$beta + sim_lake_effects$lake_offset[sim_lake_effects$lake == lake] + sim_trip_effects$trip_offset[sim_trip_effects$trip == trip] + sim_id_effects$id_offset[sim_id_effects$id == id]
      
      # compute alphas
      alphas = sim_main_effects_alpha$alpha + sim_lake_effects_alpha$lake_offset[sim_lake_effects_alpha$lake == lake] + sim_trip_effects_alpha$trip_offset[sim_trip_effects_alpha$trip == trip] + sim_id_effects_alpha$id_offset[sim_id_effects_alpha$id == id]
      
      # initialize local reward
      local_reward = 0
      
      # compute leaving probability for each time bin; simulate until leaving
      
      # angling time
      t=2;
      leave = 0;
      while (leave == 0){
        
        # update local reward based on catch
        if (cumulative_catch_tmp[t]>cumulative_catch_tmp[t-1]){
          if (cumulative_catch_tmp[t] == 1){
            # if first fish, set local reward to waiting time for first fish
            local_reward[t] = time_since_event_tmp[t] + 1;
          } else {
            # if 2nd, 3rd, ... fish update local reward according to resc wagner
            local_reward[t] = (1-inv_logit(alphas[2])) * local_reward[t-1] + inv_logit(alphas[2]) * (time_since_event_tmp[t-1] + 1);
          }
          
          # if no fish has been caught carry local reward from last timestamp
        } else {
          local_reward[t] = local_reward[t-1];
        }
        
        # update global reward
        global_reward[trip_length_counter] = global_reward[trip_length_counter-1];
        
        # compute predictors
        if (cumulative_catch_tmp[t] == 0){
          
          predictors = c(1,
                         0,
                         time_since_event_tmp[t],
                         global_reward[trip_length_counter],
                         local_reward[t],
                         success_feature_tmp[t],
                         loss_feature_tmp[t],
                         roughness_feature_tmp[t],
                         social_feature_tmp[t])
          
          
          
        } else {
          
          predictors = c(1,
                         1,
                         time_since_event_tmp[t],
                         global_reward[trip_length_counter],
                         local_reward[t],
                         success_feature_tmp[t],
                         loss_feature_tmp[t],
                         roughness_feature_tmp[t],
                         social_feature_tmp[t])
          
        }
        
        # multiply weights by predictors
        p = inv_logit(weights %*% predictors)
        
        # sample state
        leave = rbinom(1, 1, p)
        
        # angling time can be no longer than 3 hrs
        if (t > 300){
          leave = 1
        }
        
        # add one to trip length
        trip_length_counter = trip_length_counter + 1
        
        # angling time counter
        t = t + 1
        
        # if leave is 1, update global reward rate
        # if no fish has been caught at last spot
        if (cumulative_catch_tmp[t-1] == 0){
          
          # if angling time at last spot exceeds global estimate, update globale estimate according to resc wagner
          if ((t-1)>=global_reward[trip_length_counter-1]){
            global_reward[trip_length_counter] = (1-inv_logit(alphas[1])) * global_reward[trip_length_counter-1] + inv_logit(alphas[1]) * (t-1);
            # otherwise, angling time bears no information
          } else {
            global_reward[trip_length_counter] = global_reward[trip_length_counter-1];
          }
          
          # if fish was caught at last spot, update global estimate with average waiting time for fish at last spot (= angling_time / (cumulative_catch + 1))
        } else if (cumulative_catch_tmp[t-1]>0){
          global_reward[trip_length_counter] = (1-inv_logit(alphas[1])) * global_reward[trip_length_counter-1] + inv_logit(alphas[1]) * ( ((t-1)/(cumulative_catch_tmp[t-1] + 1))); 
        }
        
      }
      
      # save spot info
      res[[s]] = data.frame(iter = iter,
                            trip = j,
                            spot = s,
                            id = id,
                            lake = lake,
                            time_sim = t-1,
                            time_data = max(angling_time_tmp),
                            gut_sim = time_since_event_tmp[t-1],
                            gut_data = time_since_event_tmp[max(angling_time_tmp)],
                            catch = sum(catch_tmp),
                            catch_sim = sum(catch_tmp[1:(t-1)]),
                            iter = iter,
                            social_feature = mean(social_feature_tmp[1:(t-1)]),
                            success_feature = mean(success_feature_tmp[1:(t-1)]),
                            loss_feature = mean(loss_feature_tmp[1:(t-1)]),
                            roughness_feature = mean(roughness_feature_tmp[1:(t-1)]),
                            local_reward = local_reward[t-1],
                            global_reward = global_reward[trip_length_counter],
                            lake_effect_beta1 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][1],
                            lake_effect_beta2 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][2],
                            lake_effect_beta3 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][3],
                            lake_effect_beta4 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][4],
                            lake_effect_beta5 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][5],
                            lake_effect_beta6 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][6],
                            lake_effect_beta7 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][7],
                            lake_effect_beta8 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][8],
                            lake_effect_beta9 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][9],
                            
                            trip_effect_beta1 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][1],
                            trip_effect_beta2 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][2],
                            trip_effect_beta3 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][3],
                            trip_effect_beta4 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][4],
                            trip_effect_beta5 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][5],
                            trip_effect_beta6 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][6],
                            trip_effect_beta7 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][7],
                            trip_effect_beta8 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][8],
                            trip_effect_beta9 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][9],
                            
                            id_effect_beta1 = sim_id_effects$id_offset[sim_id_effects$id == id][1],
                            id_effect_beta2 = sim_id_effects$id_offset[sim_id_effects$id == id][2],
                            id_effect_beta3 = sim_id_effects$id_offset[sim_id_effects$id == id][3],
                            id_effect_beta4 = sim_id_effects$id_offset[sim_id_effects$id == id][4],
                            id_effect_beta5 = sim_id_effects$id_offset[sim_id_effects$id == id][5],
                            id_effect_beta6 = sim_id_effects$id_offset[sim_id_effects$id == id][6],
                            id_effect_beta7 = sim_id_effects$id_offset[sim_id_effects$id == id][7],
                            id_effect_beta8 = sim_id_effects$id_offset[sim_id_effects$id == id][8],
                            id_effect_beta9 = sim_id_effects$id_offset[sim_id_effects$id == id][9],
                            
                            lake_effect_alpha1 = sim_lake_effects_alpha$lake_offset[sim_lake_effects_alpha$lake == lake][1],
                            trip_effect_alpha1 = sim_trip_effects_alpha$trip_offset[sim_trip_effects_alpha$trip == trip][1],
                            id_effect_alpha1 = sim_id_effects_alpha$id_offset[sim_id_effects_alpha$id == id][1],
                            
                            lake_effect_alpha2 = sim_lake_effects_alpha$lake_offset[sim_lake_effects_alpha$lake == lake][2],
                            trip_effect_alpha2 = sim_trip_effects_alpha$trip_offset[sim_trip_effects_alpha$trip == trip][2],
                            id_effect_alpha2 = sim_id_effects_alpha$id_offset[sim_id_effects_alpha$id == id][2])
      
    }
    
    # trip level results
    res = do.call(rbind, res)
    
    # save results for all trips
    trip_res[[c]] = res
    c=c+1
  }
  
  trip_res = do.call(rbind, trip_res)
  
  return(trip_res)
  
}

################################################################################
sim_globallocal_patch_leaving_model <- function(iter, 
                                                stan_data, 
                                                model){
  
  # extract parameters
  sim_main_effects = model %>% 
    spread_draws(beta[parameter]) %>% 
    filter(parameter < 6) %>% 
    group_by(parameter) %>% 
    summarize(beta = mean(beta))
  
  sim_main_effects_alpha = model %>% 
    spread_draws(alpha[parameter]) %>% 
    group_by(parameter) %>% 
    summarize(alpha = mean(alpha))
  
  sim_lake_effects = model %>% 
    spread_draws(v_lake[lake, parameter]) %>% 
    filter(parameter < 6) %>% 
    group_by(parameter, lake) %>% 
    summarize(lake_offset = mean(v_lake))
  
  sim_lake_effects_alpha = model %>% 
    spread_draws(v_lake[lake, parameter]) %>% 
    filter(parameter >= 6) %>% 
    mutate(parameter = parameter - 5) %>% 
    group_by(parameter, lake) %>% 
    summarize(lake_offset = mean(v_lake))
  
  sim_trip_effects = model %>% 
    spread_draws(v_trip[trip, parameter]) %>% 
    filter(parameter < 6) %>% 
    group_by(parameter, trip) %>% 
    summarize(trip_offset = mean(v_trip))
  
  sim_trip_effects_alpha = model %>% 
    spread_draws(v_trip[trip, parameter]) %>% 
    filter(parameter >= 6) %>%
    mutate(parameter = parameter - 5) %>% 
    group_by(parameter, trip) %>% 
    summarize(trip_offset = mean(v_trip))
  
  sim_id_effects = model %>% 
    spread_draws(v_id[id, parameter]) %>% 
    filter(parameter < 6) %>% 
    group_by(parameter, id) %>% 
    summarize(id_offset = mean(v_id))
  
  sim_id_effects_alpha = model %>% 
    spread_draws(v_id[id, parameter]) %>% 
    filter(parameter >= 6) %>%
    mutate(parameter = parameter - 5) %>% 
    group_by(parameter, id) %>% 
    summarize(id_offset = mean(v_id))
  
  
  
  trip_res = list()
  c = 1
  for (j in stan_data$trips){
    
    # extract trip level data
    start = stan_data$trip_start[j]
    end = stan_data$trip_end[j]
    
    # extract data
    state = stan_data$state[start:end]
    catch = stan_data$catch[start:end]
    cumulative_catch = stan_data$cumulative_catch[start:end]
    angling_time = stan_data$angling_time[start:end]
    time_since_event = stan_data$time_since_event[start:end]
    
    # create spot identifier
    spot = c(1,lag(cumsum(state) + 1)[-1])
    
    # extract parameter identifiers
    lake = stan_data$lakes[start:end][1]
    trip = j
    id = stan_data$ids[start:end][1]
    
    # preallocate results
    res = list()
    
    # loop through all spots
    # counter for trip length
    trip_length_counter = 2
    
    # initialize local and global reward
    global_reward = stan_data$initial_values[j]
    
    # loop through all angling spots    
    for (s in unique(spot)){
      
      # subset predictors
      catch_tmp = catch[spot == s]
      time_since_event_tmp = time_since_event[spot == s]
      angling_time_tmp = angling_time[spot == s]
      
      # extend time varying predictors
      catch_tmp = c(catch_tmp, rep(0, 1500 - length(catch_tmp)))
      time_since_event_tmp = c(time_since_event_tmp, seq(from = time_since_event_tmp[length(time_since_event_tmp)]+1, to = 1500, by = 1))[1:1500]
      cumulative_catch_tmp = cumsum(catch_tmp)
      
      # compute regression weights
      weights = sim_main_effects$beta + sim_lake_effects$lake_offset[sim_lake_effects$lake == lake] + sim_trip_effects$trip_offset[sim_trip_effects$trip == trip] + sim_id_effects$id_offset[sim_id_effects$id == id]
      
      # compute alphas
      alphas = sim_main_effects_alpha$alpha + sim_lake_effects_alpha$lake_offset[sim_lake_effects_alpha$lake == lake] + sim_trip_effects_alpha$trip_offset[sim_trip_effects_alpha$trip == trip] + sim_id_effects_alpha$id_offset[sim_id_effects_alpha$id == id]
      
      # initialize local reward
      local_reward = 0
      
      # compute leaving probability for each time bin; simulate until leaving
      
      # angling time
      t=2;
      leave = 0;
      while (leave == 0){
        
        # update local reward based on catch
        if (cumulative_catch_tmp[t]>cumulative_catch_tmp[t-1]){
          if (cumulative_catch_tmp[t] == 1){
            # if first fish, set local reward to waiting time for first fish
            local_reward[t] = time_since_event_tmp[t] + 1;
          } else {
            # if 2nd, 3rd, ... fish update local reward according to resc wagner
            local_reward[t] = (1-inv_logit(alphas[2])) * local_reward[t-1] + inv_logit(alphas[2]) * (time_since_event_tmp[t-1] + 1);
          }
          
          # if no fish has been caught carry local reward from last timestamp
        } else {
          local_reward[t] = local_reward[t-1];
        }
        
        # update global reward
        global_reward[trip_length_counter] = global_reward[trip_length_counter-1];
        
        # compute predictors
        if (cumulative_catch_tmp[t] == 0){
          
          predictors = c(1,
                         0,
                         time_since_event_tmp[t],
                         global_reward[trip_length_counter],
                         local_reward[t])
          
          
          
        } else {
          
          predictors = c(1,
                         1,
                         time_since_event_tmp[t],
                         global_reward[trip_length_counter],
                         local_reward[t])
          
        }
        
        # multiply weights by predictors
        p = inv_logit(weights %*% predictors)
        
        # sample state
        leave = rbinom(1, 1, p)
        
        # angling time can be no longer than 3 hrs
        if (t > 300){
          leave = 1
        }
        
        # add one to trip length
        trip_length_counter = trip_length_counter + 1
        
        # angling time counter
        t = t + 1
        
        # if leave is 1, update global reward rate
        # if no fish has been caught at last spot
        if (cumulative_catch_tmp[t-1] == 0){
          
          # if angling time at last spot exceeds global estimate, update globale estimate according to resc wagner
          if ((t-1)>=global_reward[trip_length_counter-1]){
            global_reward[trip_length_counter] = (1-inv_logit(alphas[1])) * global_reward[trip_length_counter-1] + inv_logit(alphas[1]) * (t-1);
            # otherwise, angling time bears no information
          } else {
            global_reward[trip_length_counter] = global_reward[trip_length_counter-1];
          }
          
          # if fish was caught at last spot, update global estimate with average waiting time for fish at last spot (= angling_time / (cumulative_catch + 1))
        } else if (cumulative_catch_tmp[t-1]>0){
          global_reward[trip_length_counter] = (1-inv_logit(alphas[1])) * global_reward[trip_length_counter-1] + inv_logit(alphas[1]) * ( ((t-1)/(cumulative_catch_tmp[t-1] + 1))); 
        }
        
      }
      
      # save spot info
      res[[s]] = data.frame(iter = iter,
                            trip = j,
                            spot = s,
                            id = id,
                            lake = lake,
                            time_sim = t-1,
                            time_data = max(angling_time_tmp),
                            gut_sim = time_since_event_tmp[t-1],
                            gut_data = time_since_event_tmp[max(angling_time_tmp)],
                            catch = sum(catch_tmp),
                            catch_sim = sum(catch_tmp[1:(t-1)]),
                            iter = iter,
                            local_reward = local_reward[t-1],
                            global_reward = global_reward[trip_length_counter],
                            lake_effect_beta1 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][1],
                            lake_effect_beta2 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][2],
                            lake_effect_beta3 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][3],
                            lake_effect_beta4 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][4],
                            lake_effect_beta5 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][5],
                            
                            trip_effect_beta1 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][1],
                            trip_effect_beta2 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][2],
                            trip_effect_beta3 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][3],
                            trip_effect_beta4 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][4],
                            trip_effect_beta5 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][5],
                            
                            id_effect_beta1 = sim_id_effects$id_offset[sim_id_effects$id == id][1],
                            id_effect_beta2 = sim_id_effects$id_offset[sim_id_effects$id == id][2],
                            id_effect_beta3 = sim_id_effects$id_offset[sim_id_effects$id == id][3],
                            id_effect_beta4 = sim_id_effects$id_offset[sim_id_effects$id == id][4],
                            id_effect_beta5 = sim_id_effects$id_offset[sim_id_effects$id == id][5],
                            
                            lake_effect_alpha1 = sim_lake_effects_alpha$lake_offset[sim_lake_effects_alpha$lake == lake][1],
                            trip_effect_alpha1 = sim_trip_effects_alpha$trip_offset[sim_trip_effects_alpha$trip == trip][1],
                            id_effect_alpha1 = sim_id_effects_alpha$id_offset[sim_id_effects_alpha$id == id][1],
                            
                            lake_effect_alpha2 = sim_lake_effects_alpha$lake_offset[sim_lake_effects_alpha$lake == lake][2],
                            trip_effect_alpha2 = sim_trip_effects_alpha$trip_offset[sim_trip_effects_alpha$trip == trip][2],
                            id_effect_alpha2 = sim_id_effects_alpha$id_offset[sim_id_effects_alpha$id == id][2])
      
    }
    
    # trip level results
    res = do.call(rbind, res)
    
    # save results for all trips
    trip_res[[c]] = res
    c=c+1
  }
  
  trip_res = do.call(rbind, trip_res)
  
  return(trip_res)
  
}

################################################################################
sim_patch_discovery_patch_leaving_model <- function(iter, 
                                                    stan_data, 
                                                    model){
  
  # extract parameters
  sim_main_effects = model %>% 
    spread_draws(beta[parameter]) %>% 
    filter(parameter < 4) %>% 
    group_by(parameter) %>% 
    summarize(beta = mean(beta))
  
  sim_lake_effects = model %>% 
    spread_draws(v_lake[lake, parameter]) %>% 
    filter(parameter < 4) %>% 
    group_by(parameter, lake) %>% 
    summarize(lake_offset = mean(v_lake))
  
  sim_trip_effects = model %>% 
    spread_draws(v_trip[trip, parameter]) %>% 
    filter(parameter < 4) %>% 
    group_by(parameter, trip) %>% 
    summarize(trip_offset = mean(v_trip))
  
  sim_id_effects = model %>% 
    spread_draws(v_id[id, parameter]) %>% 
    filter(parameter < 4) %>% 
    group_by(parameter, id) %>% 
    summarize(id_offset = mean(v_id))
  
  trip_res = list()
  c = 1
  for (j in stan_data$trips){
    
    # extract trip level data
    start = stan_data$trip_start[j]
    end = stan_data$trip_end[j]
    
    # extract data
    state = stan_data$state[start:end]
    catch = stan_data$catch[start:end]
    cumulative_catch = stan_data$cumulative_catch[start:end]
    angling_time = stan_data$angling_time[start:end]
    time_since_event = stan_data$time_since_event[start:end]
    
    # create spot identifier
    spot = c(1,lag(cumsum(state) + 1)[-1])
    
    # extract parameter identifiers
    lake = stan_data$lakes[start:end][1]
    trip = j
    id = stan_data$ids[start:end][1]
    
    # preallocate results
    res = list()
    
    # loop through all spots
    # counter for trip length
    trip_length_counter = 2
    
    # loop through all angling spots    
    for (s in unique(spot)){
      
      # subset predictors
      catch_tmp = catch[spot == s]
      time_since_event_tmp = time_since_event[spot == s]
      angling_time_tmp = angling_time[spot == s]
      
      # extend time varying predictors
      catch_tmp = c(catch_tmp, rep(0, 1500 - length(catch_tmp)))
      time_since_event_tmp = c(time_since_event_tmp, seq(from = time_since_event_tmp[length(time_since_event_tmp)]+1, to = 1500, by = 1))[1:1500]
      cumulative_catch_tmp = cumsum(catch_tmp)
      
      # compute regression weights
      weights = sim_main_effects$beta + sim_lake_effects$lake_offset[sim_lake_effects$lake == lake] + sim_trip_effects$trip_offset[sim_trip_effects$trip == trip] + sim_id_effects$id_offset[sim_id_effects$id == id]
      
      # compute leaving probability for each time bin; simulate until leaving
      
      # angling time
      t=2;
      leave = 0;
      while (leave == 0){
        
        # compute predictors
        if (cumulative_catch_tmp[t] == 0){
          
          predictors = c(1,
                         0,
                         time_since_event_tmp[t])
          
          
          
        } else {
          
          predictors = c(1,
                         1,
                         time_since_event_tmp[t])
          
        }
        
        # multiply weights by predictors
        p = inv_logit(weights %*% predictors)
        
        # sample state
        leave = rbinom(1, 1, p)
        
        # angling time can be no longer than 3 hrs
        if (t > 300){
          leave = 1
        }
        
        # add one to trip length
        trip_length_counter = trip_length_counter + 1
        
        # angling time counter
        t = t + 1
        
      }  
      
      # save spot info
      res[[s]] = data.frame(iter = iter,
                            trip = j,
                            spot = s,
                            id = id,
                            lake = lake,
                            time_sim = t-1,
                            time_data = max(angling_time_tmp),
                            gut_sim = time_since_event_tmp[t-1],
                            gut_data = time_since_event_tmp[max(angling_time_tmp)],
                            catch = sum(catch_tmp),
                            catch_sim = sum(catch_tmp[1:(t-1)]),
                            iter = iter,
                            lake_effect_beta1 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][1],
                            lake_effect_beta2 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][2],
                            lake_effect_beta3 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][3],
                            
                            trip_effect_beta1 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][1],
                            trip_effect_beta2 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][2],
                            trip_effect_beta3 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][3],
                            
                            id_effect_beta1 = sim_id_effects$id_offset[sim_id_effects$id == id][1],
                            id_effect_beta2 = sim_id_effects$id_offset[sim_id_effects$id == id][2],
                            id_effect_beta3 = sim_id_effects$id_offset[sim_id_effects$id == id][3])
      
    }
    
    # trip level results
    res = do.call(rbind, res)
    
    # save results for all trips
    trip_res[[c]] = res
    c=c+1
  }
  
  trip_res = do.call(rbind, trip_res)
  
  return(trip_res)
  
}
################################################################################
sim_time_wo_catch_patch_leaving_model <- function(iter, 
                                                  stan_data, 
                                                  model){
  
  # extract parameters
  sim_main_effects = model %>% 
    spread_draws(beta[parameter]) %>% 
    filter(parameter < 3) %>% 
    group_by(parameter) %>% 
    summarize(beta = mean(beta))
  
  sim_lake_effects = model %>% 
    spread_draws(v_lake[lake, parameter]) %>% 
    filter(parameter < 3) %>% 
    group_by(parameter, lake) %>% 
    summarize(lake_offset = mean(v_lake))
  
  sim_trip_effects = model %>% 
    spread_draws(v_trip[trip, parameter]) %>% 
    filter(parameter < 3) %>% 
    group_by(parameter, trip) %>% 
    summarize(trip_offset = mean(v_trip))
  
  sim_id_effects = model %>% 
    spread_draws(v_id[id, parameter]) %>% 
    filter(parameter < 3) %>% 
    group_by(parameter, id) %>% 
    summarize(id_offset = mean(v_id))
  
  trip_res = list()
  c = 1
  for (j in stan_data$trips){
    
    # extract trip level data
    start = stan_data$trip_start[j]
    end = stan_data$trip_end[j]
    
    # extract data
    state = stan_data$state[start:end]
    catch = stan_data$catch[start:end]
    cumulative_catch = stan_data$cumulative_catch[start:end]
    angling_time = stan_data$angling_time[start:end]
    time_since_event = stan_data$time_since_event[start:end]
    
    # create spot identifier
    spot = c(1,lag(cumsum(state) + 1)[-1])
    
    # extract parameter identifiers
    lake = stan_data$lakes[start:end][1]
    trip = j
    id = stan_data$ids[start:end][1]
    
    # preallocate results
    res = list()
    
    # loop through all spots
    # counter for trip length
    trip_length_counter = 2
    
    # loop through all angling spots    
    for (s in unique(spot)){
      
      # subset predictors
      catch_tmp = catch[spot == s]
      time_since_event_tmp = time_since_event[spot == s]
      angling_time_tmp = angling_time[spot == s]
      
      # extend time varying predictors
      catch_tmp = c(catch_tmp, rep(0, 1500 - length(catch_tmp)))
      time_since_event_tmp = c(time_since_event_tmp, seq(from = time_since_event_tmp[length(time_since_event_tmp)]+1, to = 1500, by = 1))[1:1500]
      cumulative_catch_tmp = cumsum(catch_tmp)
      
      # compute regression weights
      weights = sim_main_effects$beta + sim_lake_effects$lake_offset[sim_lake_effects$lake == lake] + sim_trip_effects$trip_offset[sim_trip_effects$trip == trip] + sim_id_effects$id_offset[sim_id_effects$id == id]
      
      # compute leaving probability for each time bin; simulate until leaving
      
      # angling time
      t=2;
      leave = 0;
      while (leave == 0){
        
        # compute predictors
        if (cumulative_catch_tmp[t] == 0){
          
          predictors = c(1,
                         time_since_event_tmp[t])
          
          
          
        } else {
          
          predictors = c(1,
                         time_since_event_tmp[t])
          
        }
        
        # multiply weights by predictors
        p = inv_logit(weights %*% predictors)
        
        # sample state
        leave = rbinom(1, 1, p)
        
        # angling time can be no longer than 3 hrs
        if (t > 300){
          leave = 1
        }
        
        # add one to trip length
        trip_length_counter = trip_length_counter + 1
        
        # angling time counter
        t = t + 1
        
      }  
      
      # save spot info
      res[[s]] = data.frame(iter = iter,
                            trip = j,
                            spot = s,
                            id = id,
                            lake = lake,
                            time_sim = t-1,
                            time_data = max(angling_time_tmp),
                            gut_sim = time_since_event_tmp[t-1],
                            gut_data = time_since_event_tmp[max(angling_time_tmp)],
                            catch = sum(catch_tmp),
                            catch_sim = sum(catch_tmp[1:(t-1)]),
                            iter = iter,
                            lake_effect_beta1 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][1],
                            lake_effect_beta2 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][2],
                            
                            trip_effect_beta1 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][1],
                            trip_effect_beta2 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][2],
                            
                            id_effect_beta1 = sim_id_effects$id_offset[sim_id_effects$id == id][1],
                            id_effect_beta2 = sim_id_effects$id_offset[sim_id_effects$id == id][2])
      
    }
    
    # trip level results
    res = do.call(rbind, res)
    
    # save results for all trips
    trip_res[[c]] = res
    c=c+1
  }
  
  trip_res = do.call(rbind, trip_res)
  
  return(trip_res)
  
}
################################################################################
sim_baseline_patch_leaving_model <- function(iter, 
                                             stan_data, 
                                             model){
  
  # extract parameters
  sim_main_effects = model %>% 
    spread_draws(beta[parameter]) %>% 
    filter(parameter < 2) %>% 
    group_by(parameter) %>% 
    summarize(beta = mean(beta))
  
  sim_lake_effects = model %>% 
    spread_draws(v_lake[lake, parameter]) %>% 
    filter(parameter < 2) %>% 
    group_by(parameter, lake) %>% 
    summarize(lake_offset = mean(v_lake))
  
  sim_trip_effects = model %>% 
    spread_draws(v_trip[trip, parameter]) %>% 
    filter(parameter < 2) %>% 
    group_by(parameter, trip) %>% 
    summarize(trip_offset = mean(v_trip))
  
  sim_id_effects = model %>% 
    spread_draws(v_id[id, parameter]) %>% 
    filter(parameter < 2) %>% 
    group_by(parameter, id) %>% 
    summarize(id_offset = mean(v_id))
  
  trip_res = list()
  c = 1
  for (j in stan_data$trips){
    
    # extract trip level data
    start = stan_data$trip_start[j]
    end = stan_data$trip_end[j]
    
    # extract data
    state = stan_data$state[start:end]
    catch = stan_data$catch[start:end]
    cumulative_catch = stan_data$cumulative_catch[start:end]
    angling_time = stan_data$angling_time[start:end]
    time_since_event = stan_data$time_since_event[start:end]
    
    # create spot identifier
    spot = c(1,lag(cumsum(state) + 1)[-1])
    
    # extract parameter identifiers
    lake = stan_data$lakes[start:end][1]
    trip = j
    id = stan_data$ids[start:end][1]
    
    # preallocate results
    res = list()
    
    # loop through all spots
    # counter for trip length
    trip_length_counter = 2
    
    # loop through all angling spots    
    for (s in unique(spot)){
      
      # subset predictors
      catch_tmp = catch[spot == s]
      time_since_event_tmp = time_since_event[spot == s]
      angling_time_tmp = angling_time[spot == s]
      
      # extend time varying predictors
      catch_tmp = c(catch_tmp, rep(0, 1500 - length(catch_tmp)))
      time_since_event_tmp = c(time_since_event_tmp, seq(from = time_since_event_tmp[length(time_since_event_tmp)]+1, to = 1500, by = 1))[1:1500]
      cumulative_catch_tmp = cumsum(catch_tmp)
      
      # compute regression weights
      weights = sim_main_effects$beta + sim_lake_effects$lake_offset[sim_lake_effects$lake == lake] + sim_trip_effects$trip_offset[sim_trip_effects$trip == trip] + sim_id_effects$id_offset[sim_id_effects$id == id]
      
      # compute leaving probability for each time bin; simulate until leaving
      
      # angling time
      t=2;
      leave = 0;
      while (leave == 0){
        
        # compute predictors
        if (cumulative_catch_tmp[t] == 0){
          
          predictors = c(1)
          
          
          
        } else {
          
          predictors = c(1)
          
        }
        
        # multiply weights by predictors
        p = inv_logit(weights %*% predictors)
        
        # sample state
        leave = rbinom(1, 1, p)
        
        # angling time can be no longer than 3 hrs
        if (t > 300){
          leave = 1
        }
        
        # add one to trip length
        trip_length_counter = trip_length_counter + 1
        
        # angling time counter
        t = t + 1
        
      }  
      
      # save spot info
      res[[s]] = data.frame(iter = iter,
                            trip = j,
                            spot = s,
                            id = id,
                            lake = lake,
                            time_sim = t-1,
                            time_data = max(angling_time_tmp),
                            gut_sim = time_since_event_tmp[t-1],
                            gut_data = time_since_event_tmp[max(angling_time_tmp)],
                            catch = sum(catch_tmp),
                            catch_sim = sum(catch_tmp[1:(t-1)]),
                            iter = iter,
                            lake_effect_beta1 = sim_lake_effects$lake_offset[sim_lake_effects$lake == lake][1],
                            
                            trip_effect_beta1 = sim_trip_effects$trip_offset[sim_trip_effects$trip == trip][1],
                            
                            id_effect_beta1 = sim_id_effects$id_offset[sim_id_effects$id == id][1])
      
    }
    
    # trip level results
    res = do.call(rbind, res)
    
    # save results for all trips
    trip_res[[c]] = res
    c=c+1
  }
  
  trip_res = do.call(rbind, trip_res)
  
  return(trip_res)
  
}
################################################################################
################################################################################
# run simulations
################################################################################
# baseline model
n_iter = 50
sim_res_baseline = list()
start_time = Sys.time()
for (i in 1:n_iter){
  sim_res_baseline[[i]] <- sim_baseline_patch_leaving_model(i,
                                                            stan_data = stan_data,
                                                            model = a_baseline)
  print(i)
}
end_time = Sys.time()
end_time-start_time

# save results
saveRDS(sim_res_baseline, file = "utils/data/processed_data/patch_leaving_simulations/a_sim_baseline.rds")

# time wo catch model
sim_res_time_wo_catch = list()
start_time = Sys.time()
for (i in 1:n_iter){
  sim_res_time_wo_catch[[i]] <- sim_time_wo_catch_patch_leaving_model(i,
                                                                      stan_data = stan_data,
                                                                      model = b_time_wo_catch)
  print(i)
}
end_time = Sys.time()
end_time-start_time

# save results
saveRDS(sim_res_time_wo_catch, file = "utils/data/processed_data/patch_leaving_simulations/b_sim_time_wo_catch.rds")

# patch discovery model
sim_res_patch_discovery = list()
start_time = Sys.time()
for (i in 1:n_iter){
  sim_res_patch_discovery[[i]] <- sim_patch_discovery_patch_leaving_model(i,
                                                                          stan_data = stan_data,
                                                                          model = c_patch_discovery)
  print(i)
}
end_time = Sys.time()
end_time-start_time

# save results
saveRDS(sim_res_patch_discovery, file = "utils/data/processed_data/patch_leaving_simulations/c_sim_patch_discovery.rds")

# global local updating model
sim_res_global_local = list()
start_time = Sys.time()
for (i in 1:n_iter){
  sim_res_global_local[[i]] <- sim_globallocal_patch_leaving_model(i,
                                                                   stan_data = stan_data,
                                                                   model = d_globallocal)
  print(i)
}
end_time = Sys.time()
end_time-start_time

# save results
saveRDS(sim_res_global_local, file = "utils/data/processed_data/patch_leaving_simulations/d_sim_global_local.rds")

# spatial features
sim_res_spatial_features = list()
start_time = Sys.time()
for (i in 1:n_iter){
  sim_res_spatial_features[[i]] <- sim_spatial_patch_leaving_model(i,
                                                                   stan_data = stan_data,
                                                                   model = e_spatial_features)
  print(i)
}
end_time = Sys.time()
end_time-start_time

# save results
saveRDS(sim_res_spatial_features, file = "utils/data/processed_data/patch_leaving_simulations/e_sim_spatial_features.rds")

}
################################################################################
# END
################################################################################