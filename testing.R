#! /usr/bin/env Rscript
input_args = commandArgs(trailingOnly = TRUE)

outdir <- input_args[1] # "outputs/testing/"

msg <- file(paste0(outdir, "logfile_testing.txt"), open="wt")
sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("funcs.R")))

stopwatch <- rep(0,2) %>% set_names(c("start_time", "end_time"))
stopwatch[1] <- Sys.time()

outputDetails(paste0("\n||---------------- Testing results ----------------||\n", 
                     "Started process at: ", Sys.time()))

time1_raw <- input_args[3] %>% readBaseData(., 1) # "data/timepoint1_data.csv"
time2_raw <- input_args[4] %>% readBaseData(., 2) # "data/timepoint2_data.csv"
inpdir <- input_args[2]                           # "outputs/height_data/"

message("Part A")

# the TP1 cluster assignments
t1_melted <- time1_raw %>% melt(id = "isolate") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% mutate(across(variable, as.integer)) %>% 
  createID(., "tp1", "variable", "value") %>% 
  set_colnames(c("isolate", "tp1_h", "tp1_cl", "tp1_id"))

# adding a column with the TP1 cluster sizes
b1 <- t1_melted %>% group_by(tp1_id) %>% 
  summarise(tp1_cl_size = n(), .groups = "drop") %>% 
  left_join(t1_melted, ., by = "tp1_id")

message("Part B")
# the TP2 cluster assignments
t2_melted <- time2_raw %>% melt(id = "isolate") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% mutate(across(variable, as.integer)) %>% 
  createID(., "tp2", "variable", "value") %>% 
  set_colnames(c("isolate", "tp2_h", "tp2_cl", "tp2_id"))

# adding a column with the TP2 cluster sizes
b2 <- t2_melted %>% group_by(tp2_id) %>% 
  summarise(tp2_cl_size = n(), .groups = "drop") %>% 
  left_join(t2_melted, ., by = "tp2_id")

message("Part C")
# going to run the test for every TP1 threshold
heights <- b1$tp1_h %>% unique()

# how many clusters is 5% of those at each threshold?
# we will sample and test this many clusters
tocheck <- lapply(heights, function(h) {
  b1clusters <- b1 %>% filter(tp1_h == h) %>% pull(tp1_cl) %>% unique()
  checked <- ceiling(0.025*length(b1clusters))
  tibble(height = h, total_num_clusters = length(b1clusters), num_clusters_tested = checked)
}) %>% bind_rows()

saveRDS(stopwatch, paste0(outdir, "start_time.Rds"))

message("\nPart D")
for (x in heights) {
  # Part 1: For this TP1 height, we select 100 clusters at random for testing
  height_x <- b1 %>% filter(tp1_h == x)
  b1clusters <- height_x$tp1_cl %>% unique()
  
  sampled_clusters <- tocheck %>% 
    filter(height == x) %>% 
    pull(num_clusters_tested) %>% 
    sample(b1clusters, ., replace = FALSE)
  
  checked_sizes <- height_x %>% filter(tp1_cl %in% sampled_clusters) %>% pull(tp1_cl_size) %>% unique()
  actual_sizes <- height_x %>% pull(tp1_cl_size) %>% unique()

  # Part 2: Extracting the collected tracking and growth data for the given height, then filtering 
  #         and sorting so we will only compare for 100 randomly selected clusters at this height
  filename <- paste0(inpdir, "h", x, ".Rds")
  
  if (file.exists(filename)) {
    collected_data <- readRDS(filename) %>% 
      select(tp1_h, tp1_cl, tp1_cl_size, tp2_h, tp2_cl, tp2_cl_size) %>% 
      filter(tp1_cl %in% sampled_clusters) %>% 
      arrange(tp1_cl, tp2_h, tp2_cl) %>% 
      mutate(across(tp2_cl_size, as.integer))  
  }else {
    stop(paste0("Base data file not found for height h_", x))
  }

  # Part 3: For each of the sampled clusters, we manually filter the TP2 set to find all clusters 
  #         that contain the isolates from that TP1 cluster, then count to see which of these TP2 
  #         clusters contain (at least) all of the TP1 isolates from that cluster.
  #         We then bind the resulting TP2 cases to the TP1 cluster, sizes and all.
  actual_data <- lapply(1:length(sampled_clusters), function(i) {
    # TP1 cluster and the isolates in it
    c1 <- height_x %>% filter(tp1_cl == sampled_clusters[i])
    # All TP2 clusters with the isolates from the sampled TP1 cluster
    basecase <- b2 %>% filter(isolate %in% c1$isolate)
    
    sizes <- basecase %>% group_by(tp2_id) %>% summarise(cl_size = n(), .groups = "drop") %>% 
      left_join(basecase, ., by = "tp2_id") %>% filter(cl_size >= nrow(c1)) %>% 
      select(-isolate, -cl_size) %>% unique()
    
    c1 %>% select(-isolate) %>% unique() %>% bind_cols(., sizes) %>% 
      arrange(tp1_cl, tp2_h, tp2_cl) %>% return()
    
  }) %>% bind_rows() %>% 
    select(colnames(collected_data)) %>% 
    arrange(tp1_cl, tp2_h, tp2_cl) %>% 
    mutate(across(tp2_cl_size, as.integer))
  
  if (!exists(actual_data)) {stop(paste0("Actual data not correctly filtered."))}
  
  # Part 4: If our collected data is the same as the actual matching (for these 100 
  #         clusters), then we move onto the next height.
  #         Otherwise, we make a note in the log and save the results for further analysis.
  if (identical(collected_data, actual_data)) {
    message(paste0("TP1_h", x, " - sampled ", length(sampled_clusters), " of ", length(b1clusters), 
                 " clusters, with sizes in range (", min(checked_sizes), ", ", max(checked_sizes), 
                 "). Actual size range (", min(actual_sizes), ", ", max(actual_sizes), "): success."))
  }else {
    message(paste0("Test case at TP1_h", x, " failed. Check clusters: \n", toString(sampled_clusters)))
    saveRDS(actual_data, paste0(outdir, "h", x, "_actual_data.Rds"))
    saveRDS(collected_data, paste0(outdir, "h", x, "_collected_data.Rds"))
  }
}
# h70 - checked
outputDetails(paste0("\n||---------------- Testing results ----------------||\n", 
                     "Finished process at: ", Sys.time()))
stopwatch[2] <- Sys.time()
saveRDS(stopwatch, paste0(outdir, "stopwatch.Rds"))

message(timeTaken("testing", stopwatch))

# time1_raw %>% select('isolate', '1') %>% rename(h1 = '1') %>% filter(h1 == 105)
# idx <- time1_raw %>% select('isolate', '1') %>% 
#   rename(h1 = '1') %>% filter(h1 == 105) %>% pull(isolate)
