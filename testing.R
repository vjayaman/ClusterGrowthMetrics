#! /usr/bin/env Rscript
# input_args = commandArgs(trailingOnly = TRUE)

msg <- file("outputs/log_testing.txt", open="wt")
sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("funcs.R")))

outputDetails(paste0("\n||---------------- Testing results ----------------||\n", 
                     "Started process at: ", Sys.time()))

time1_raw <- "data/timepoint1_data.csv" %>% readBaseData(., 1)
time2_raw <- "data/timepoint2_data.csv" %>% readBaseData(., 2)

message("Part A")

# the TP1 cluster assignments
t1_melted <- time1_raw %>% 
  melt(id = "isolate") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% 
  mutate(across(variable, as.integer)) %>% 
  createID(., "tp1", "variable", "value") %>% 
  set_colnames(c("isolate", "tp1_h", "tp1_cl", "tp1_id"))

# adding a column with the TP1 cluster sizes
b1 <- t1_melted %>% 
  group_by(tp1_id) %>% 
  summarise(tp1_cl_size = n(), .groups = "drop") %>% 
  left_join(t1_melted, ., by = "tp1_id")

message("Part B")
# the TP2 cluster assignments
t2_melted <- time2_raw %>% 
  melt(id = "isolate") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% 
  mutate(across(variable, as.integer)) %>% 
  createID(., "tp2", "variable", "value") %>% 
  set_colnames(c("isolate", "tp2_h", "tp2_cl", "tp2_id"))

# adding a column with the TP2 cluster sizes
b2 <- t2_melted %>% 
  group_by(tp2_id) %>% 
  summarise(tp2_cl_size = n(), .groups = "drop") %>% 
  left_join(t2_melted, ., by = "tp2_id")

message("Part C")
# going to run the test for every TP1 threshold
heights <- b1$tp1_h %>% unique()

# how many clusters is 5% of those at each threshold?
# we will sample and test this many clusters
tocheck <- lapply(heights, function(h) {
  b1clusters <- b1 %>% filter(tp1_h == h) %>% pull(tp1_cl) %>% unique()
  checked <- 0.05*length(b1clusters) %>% ceiling()
  tibble(height = h, total_num_clusters = length(b1clusters), num_clusters_tested = checked)
}) %>% bind_rows()

message("\nPart D")
for (height_x in heights) {
  # Part 1: For this TP1 height, we select 100 clusters at random for testing
  hx <- b1 %>% filter(tp1_h == height_x)
  b1clusters <- hx$tp1_cl %>% unique()
  
  sampled_clusters <- tocheck %>% 
    filter(height == height_x) %>% 
    pull(num_clusters_tested) %>% 
    sample(b1clusters, ., replace = FALSE)
  
  checked_size_range <- hx %>% filter(tp1_cl %in% sampled_clusters) %>% pull(tp1_cl_size) %>% unique()
  actual_size_range <- hx %>% pull(tp1_cl_size) %>% unique()
  
  message(paste0("TP1_h", height_x, " - testing ", length(sampled_clusters), " of ", length(b1clusters), 
                 " clusters, with sizes in range (", min(checked_size_range), ", ", max(checked_size_range), 
                 "). Actual size range (", min(actual_size_range), ", ", max(actual_size_range), ")."))
  
  # Part 2: Extracting the collected tracking and growth data for the given height, then filtering 
  #         and sorting so we will only compare for 100 randomly selected clusters at this height
  collected_data <- readRDS(paste0("outputs/height_data/h", height_x, ".Rds")) %>% 
    rename(tp1_id = id) %>% 
    createID(., "tp2", "tp2_h", "tp2_cl") %>% 
    rename(tp2_id = id) %>% 
    select(tp1_h, tp1_cl, tp1_cl_size, tp2_h, tp2_cl, tp2_cl_size) %>% 
    filter(tp1_cl %in% sampled_clusters) %>% 
    arrange(tp1_cl, tp2_h, tp2_cl) %>% 
    mutate(across(tp2_cl_size, as.integer))

  # Part 3: For each of the sampled clusters, we manually filter the TP2 set to find 
  #         all clusters that contain the isolates from that TP1 cluster, then count 
  #         to see which of these TP2 clusters contain (at least) all of the TP1 
  #         isolates from that cluster.
  #         We then bind the resulting TP2 cases to the TP1 cluster, sizes and all
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
  
  # Part 4: If our collected data is the same as the actual matching (for these 100 
  #         clusters), then we move onto the next height.
  #         Otherwise, we make a note in the log and save the results for further analysis.
  if (!identical(collected_data, actual_data)) {
    message(paste0("Test case TP1_h", height_x, " failed."))
    saveRDS(actual_data, paste0("outputs/testing/h", height_x, "_actual_data.Rds"))
    saveRDS(collected_data, paste0("outputs/testing/h", height_x, "_collected_data.Rds"))
  }
}

outputDetails(paste0("\n||---------------- Testing results ----------------||\n", 
                     "Finished process at: ", Sys.time()))

# missed <- setdiff(actual_data, collected_data)
# saveRDS(missed, "tp1_h1_c105_missing.Rds")

# actual_data <- readRDS("outputs/testing/h1_actual_data.Rds")
# collected_data <- readRDS("outputs/testing/h1_collected_data.Rds")

# time1_raw %>% select('isolate', '1') %>% rename(h1 = '1') %>% filter(h1 == 105)
# idx <- time1_raw %>% select('isolate', '1') %>% 
#   rename(h1 = '1') %>% filter(h1 == 105) %>% pull(isolate)



