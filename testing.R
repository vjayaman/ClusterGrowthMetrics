#! /usr/bin/env Rscript
# input_args = commandArgs(trailingOnly = TRUE)

msg <- file("outputs/log_testing.txt", open="wt")
sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("functions/base_functions.R")))
source("functions/processing_functions.R")

time1_raw <- "data/timepoint1_data.csv" %>% readBaseData(., 1)
time2_raw <- "data/timepoint2_data.csv" %>% readBaseData(., 2)

# input_args = 125 2 , or 0 1
# cluster_x <- 1#input_args[2]
message("Part A")
# the TP1 cluster assignments
t1_melted <- time1_raw %>% 
  melt(id = "isolate") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% 
  mutate(across(variable, as.integer)) %>% 
  createID(., "tp1", "variable", "value") %>% 
  set_colnames(c("isolate", "tp1_h", "tp1_cl", "tp1_id"))

b1 <- t1_melted %>% 
  group_by(tp1_id) %>% 
  summarise(tp1_cl_size = n(), .groups = "drop") %>% 
  left_join(t1_melted, ., by = "tp1_id")

# the TP2 cluster assignments
t2_melted <- time2_raw %>% 
  melt(id = "isolate") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% 
  mutate(across(variable, as.integer)) %>% 
  createID(., "tp2", "variable", "value") %>% 
  set_colnames(c("isolate", "tp2_h", "tp2_cl", "tp2_id"))

b2 <- t2_melted %>% 
  group_by(tp2_id) %>% 
  summarise(tp2_cl_size = n(), .groups = "drop") %>% 
  left_join(t2_melted, ., by = "tp2_id")

heights <- b1$tp1_h %>% unique()
message("\nPart B")

for (height_x in heights) {
  message(paste0("TP1_h", height_x))
  # message("-Part 1")
  
  hx <- b1 %>% filter(tp1_h == height_x)
  b1clusters <- hx$tp1_cl %>% unique()
  
  sampled_clusters <- sample(b1clusters, 100, replace = FALSE)
  
  collected_data <- readRDS(paste0("outputs/height_data/h", height_x, ".Rds")) %>% 
    rename(tp1_id = id) %>% 
    createID(., "tp2", "tp2_h", "tp2_cl") %>% 
    rename(tp2_id = id) %>% 
    select(tp1_h, tp1_cl, tp1_cl_size, tp2_h, tp2_cl, tp2_cl_size) %>% 
    filter(tp1_cl %in% sampled_clusters) %>% 
    arrange(tp1_cl, tp2_h, tp2_cl) %>% 
    mutate(across(tp2_cl_size, as.integer))

  # message("-Part 2")
  
  actual_data <- lapply(1:length(sampled_clusters), function(i) {
    c1 <- hx %>% filter(tp1_cl == sampled_clusters[i])
    basecase <- b2 %>% filter(isolate %in% c1$isolate)
    
    sizes <- basecase %>% 
      group_by(tp2_id) %>% 
      summarise(num_c1 = n(), .groups = "drop") %>% 
      left_join(basecase, ., by = "tp2_id") %>% 
      filter(num_c1 >= nrow(c1)) %>% 
      select(-isolate, -num_c1) %>% 
      unique()
    
    c1 %>% 
      select(-isolate) %>% unique() %>% 
      bind_cols(., sizes) %>% 
      arrange(tp1_cl, tp2_h, tp2_cl) %>% return()
  }) %>% bind_rows() %>% 
    select(colnames(collected_data)) %>% 
    arrange(tp1_cl, tp2_h, tp2_cl) %>% 
    mutate(across(tp2_cl_size, as.integer))
  
  # message("-Part 3")
  
  if (!identical(collected_data, actual_data)) {
    message(paste0("Test case TP1_h", height_x, " failed."))
    saveRDS(actual_data, paste0("outputs/testing/h", height_x, "_actual_data.Rds"))
    saveRDS(collected_data, paste0("outputs/testing/h", height_x, "_collected_data.Rds"))
  }
}

# actual_data <- readRDS("outputs/testing/h1_actual_data.Rds")
# collected_data <- readRDS("outputs/testing/h1_collected_data.Rds")

# missed <- setdiff(actual_data, collected_data)
# saveRDS(missed, "tp1_h1_c105_missing.Rds")

# collected_data <- a1 %>% filter(tp1_cl == cluster_x) %>% arrange(tp1_cl, tp2_h, tp2_cl)
# actual_data <- hxclusters %>% filter(tp1_cl == cluster_x) %>% arrange(tp1_cl, tp2_h, tp2_cl)
# saveRDS(actual_data, "testseth0c1.Rds")



# idx <- time1_raw %>% select(isolate, '0') %>% slice(1) %>% pull(isolate)
# time2_raw %>% select(isolate, '0') %>% filter(isolate == idx)




# time1_raw %>% select('isolate', '1') %>% rename(h1 = '1') %>% filter(h1 == 105)
# idx <- time1_raw %>% select('isolate', '1') %>% 
#   rename(h1 = '1') %>% filter(h1 == 105) %>% pull(isolate)




# # identifiers for all h0 clusters
# tp1clustersx <- b1 %>% filter(tp1_h == x) %>% select(tp1_id) %>% unique() %>% pull()
# 
# # isos <- list()
# # starting at 4:42
# matching_x <- lapply(1:length(tp1clustersx), function(i) {
#   print(paste0(i, " / ", length(tp1clustersx)))
#   
#   cx <- tp1clustersx[i]
#   c2 <- b1 %>% filter(tp1_id == cx)
#   
#   part1 <- c2 %>% select(-isolate) %>% unique()
#   c2size <- c2$tp1_cl_size %>% unique()
# 
#   c2match <- b2 %>% 
#     filter(isolate %in% c2$isolate) %>% 
#     filter(tp2_cl_size >= c2size) %>% 
#     arrange(tp2_h, tp2_cl) %>% 
#     slice(1) %>% 
#     pull(tp2_id)
#   
#   part2 <- b2 %>% filter(tp2_id == c2match) %>% 
#     select(-isolate) %>% unique()
#   
#   bind_cols(part1, part2) %>% return()
# }) %>% bind_rows()
# # Sys.time()
# # saveRDS(matching_first_case, "outputs/testing/h0.Rds")
# # saveRDS(matching_x, "outputs/testing/h125.Rds")
# 
# # tp1clusters1 <- b1 %>% filter(tp1_h == 1) %>% select(tp1_id) %>% unique() %>% pull()
# # cx <- tp1clusters1[2]
# # c2 <- b1 %>% filter(tp1_id == cx)
# # b1 %>% filter(tp1_h == 0)
# # matching_second_case <- 
# 




# b2 <- time2_raw %>% 
#   melt(id = "isolate") %>% 
#   as_tibble() %>% 
#   mutate(across(variable, as.character)) %>% 
#   mutate(across(variable, as.integer)) %>% 
#   createID(., "tp2", "variable", "value") %>% 
#   set_colnames(c("isolate", "tp2_h", "tp2_cl", "tp2_id"))
#   
# b2 <- table(b2$tp2_id) %>% as.data.frame() %>% 
#   as_tibble() %>% set_colnames(c("tp2_id", "tp2_cl_size")) %>% 
#   mutate(across(tp2_id, as.character)) %>% 
#   left_join(b2, ., by = "tp2_id")

