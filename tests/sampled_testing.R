#! /usr/bin/env Rscript
input_args = commandArgs(trailingOnly = TRUE)

msg <- file(paste0("outputs/", "logfile_testing.txt"), open="wt")
sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("functions/tracking_functions.R")))

stopwatch <- rep(0,2) %>% set_names(c("start_time", "end_time"))
stopwatch[1] <- Sys.time()

outputDetails(paste0("\n||---------------- Testing results ----------------||\n", 
                     "Started process at: ", Sys.time()), newcat = TRUE)

time1_raw <- input_args[1] %>% readBaseData(., 1)
time2_raw <- input_args[2] %>% readBaseData(., 2)
percent_clusters <- as.double(input_args[3])
# time1_raw <- readBaseData("data/timepoint1.csv", 1)
# time2_raw <- readBaseData("data/timepoint2.csv", 2); percent_clusters <- 0.25
outputDetails("Part 1/3 - data prep ...", newcat = TRUE)

# the TP1 cluster assignments
t1_melted <- time1_raw %>% melt(id = "isolate") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% mutate(across(variable, as.integer)) %>% 
  createID(., "tp1", "variable", "value") %>% 
  set_colnames(c("isolate", "tp1_h", "tp1_cl", "tp1_id"))

# adding a column with the TP1 cluster sizes
b1 <- t1_melted %>% group_by(tp1_id) %>% 
  summarise(tp1_cl_size = n(), .groups = "drop") %>% 
  left_join(t1_melted, ., by = "tp1_id")

# the TP2 cluster assignments
t2_melted <- time2_raw %>% melt(id = "isolate") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% mutate(across(variable, as.integer)) %>% 
  createID(., "tp2", "variable", "value") %>% 
  set_colnames(c("isolate", "tp2_h", "tp2_cl", "tp2_id"))

# adding a column with the TP2 cluster sizes
b2 <- t2_melted %>% group_by(tp2_id) %>% 
  summarise(tp2_cl_size = n(), .groups = "drop") %>% 
  left_join(t2_melted, ., by = "tp2_id")

novels <- setdiff(b2$isolate, b1$isolate)

# going to run the test for every TP1 threshold
collected_data <- read.csv(file = "outputs/summary/TP1_cluster_results.csv", 
                           stringsAsFactors = FALSE, numerals = "no.loss") %>% 
  as_tibble() %>% select(2, 3, 4, 8, 9, 10, 11, 12, 13) %>% 
  set_colnames(c("tp1_h", "tp1_cl", "tp1_cl_size", "tp2_h", "tp2_cl", 
                 "tp2_cl_size", "add_TP1", "num_nov", "size_change")) %>% 
  mutate(across(tp2_cl_size, as.integer)) %>% 
  arrange(tp1_cl, tp2_h, tp2_cl)

collected_data$tp1_h %<>% charToInt(., "h")
collected_data$tp1_cl %<>% charToInt(., "c")
collected_data$tp2_h %<>% charToInt(., "h")
collected_data$tp2_cl %<>% charToInt(., "c")

heights <- collected_data$tp1_h %>% unique()

# how many clusters is 5% of those at each threshold?
# we will sample and test this many clusters
outputDetails(paste0("Part 2/3 - cluster sampling - checking sizes and tracking for ", 
                     scales::percent(percent_clusters, suffix = " %"), 
                     " of the clusters at each height."), newcat = TRUE)
tocheck <- lapply(heights, function(h) {
  b1clusters <- b1 %>% filter(tp1_h == h) %>% pull(tp1_cl) %>% unique()
  checked <- ceiling(percent_clusters*length(b1clusters))
  tibble(height = h, total_num_clusters = length(b1clusters), num_clusters_tested = checked)
}) %>% bind_rows()

outputDetails(paste0("Part 3/3 - beginning check for thresholds ", toString(heights)), newcat = TRUE)
for (x in heights) {
  outputDetails(paste0("\n   Checking ", scales::percent(percent_clusters, suffix = " %"), 
                       " of clusters at threshold ", x, " in TP1"), newcat = TRUE)
  # Part 1: For this TP1 height, we select 100 clusters at random for testing
  heightx <- b1 %>% filter(tp1_h == x)
  b1clusters <- heightx$tp1_cl %>% unique()
  
  sampled_clusters <- tocheck %>% 
    filter(height == x) %>% 
    pull(num_clusters_tested) %>% 
    sample(b1clusters, ., replace = FALSE)
  
  checked_sizes <- heightx %>% filter(tp1_cl %in% sampled_clusters) %>% pull(tp1_cl_size) %>% unique()
  actual_sizes <- heightx %>% pull(tp1_cl_size) %>% unique()

  # Part 3: For each of the sampled clusters, we manually filter the TP2 set to find all clusters 
  #         that contain the isolates from that TP1 cluster, then count to see which of these TP2 
  #         clusters contain (at least) all of the TP1 isolates from that cluster.
  #         We then bind the resulting TP2 cases to the TP1 cluster, sizes and all.
  actual_data <- lapply(1:length(sampled_clusters), function(i) {
    
    c1 <- heightx %>% filter(tp1_cl == sampled_clusters[i])
    basecase <- b2 %>% filter(isolate %in% c1$isolate)
    
    sizes <- basecase %>% group_by(tp2_id) %>% 
      summarise(num_from_c1 = n(), .groups = "drop") %>% 
      left_join(basecase, ., by = "tp2_id")
    
    # want TP2 clusters that *at least* contain all the c1 isolates
    compared <- sizes %>% 
      filter(num_from_c1 >= unique(c1$tp1_cl_size)) %>% 
      select(-isolate, -num_from_c1) %>% unique() %>% 
      arrange(tp2_h, tp2_cl) %>% slice(1)
    
    c1a <- c1 %>% select(-isolate) %>% unique() %>% bind_cols(., compared)
    
    other_than <- b1 %>% filter(tp1_id == c1a$tp1_id) %>% pull(isolate) %>% c(., novels)
    in_tp2 <- b2 %>% filter(tp2_id == c1a$tp2_id) %>% pull(isolate)
    
    c1a %>% add_column(add_TP1 = setdiff(in_tp2, other_than) %>% length(), 
                       num_nov = intersect(in_tp2, novels) %>% length()) %>% return()
  }) %>% bind_rows()
  ad <- actual_data %>% 
    mutate(across(tp2_cl_size, as.integer)) %>% 
    select(tp1_h, tp1_cl, tp1_cl_size, tp2_h, tp2_cl, tp2_cl_size, add_TP1, num_nov) %>% 
    add_column(size_change = actual_data$tp2_cl_size - actual_data$tp1_cl_size) %>% 
    arrange(tp1_cl, tp2_h, tp2_cl)  

  # Part 2: Extracting the collected tracking and growth data for the given height, then filtering 
  #         and sorting so we will only compare for 100 randomly selected clusters at this height
  cd <- collected_data %>% filter(tp1_h == x) %>% filter(tp1_cl %in% sampled_clusters)
  
  # Part 4: If our collected data is the same as the actual matching (for these 100 
  #         clusters), then we move onto the next height.
  #         Otherwise, we make a note in the log and save the results for further analysis.
  if (identical(cd, ad)) {
    outputDetails(paste0("      TP1_h", x, " - sampled ", length(sampled_clusters), " of ", 
                         length(b1clusters), " clusters, with sizes in range (", min(checked_sizes), 
                         ", ", max(checked_sizes), "). ", 
                         "\n      Actual size range (", min(actual_sizes), ", ", max(actual_sizes), ").", 
                         "\n      Testing generated data == actual data: Success"), newcat = TRUE)
  }else {
    outputDetails(paste0("      Test case at TP1_h", x, " failed. \n      Check clusters: \n", 
                         toString(sampled_clusters)), newcat = TRUE)
  }
}

outputDetails(paste0("\n||---------------- Testing results ----------------||\n", 
                     "Finished process at: ", Sys.time()))
stopwatch[2] <- Sys.time()

outputDetails(timeTaken("testing", stopwatch))

