#! /usr/bin/env Rscript
# input_args = commandArgs(trailingOnly = TRUE)

suppressWarnings(suppressPackageStartupMessages(source("functions/tracking_functions.R")))

# time1_raw <- readBaseData(input_args[1], 1) # "data/timepoint1_data.csv"
# time2_raw <- readBaseData(input_args[2], 2) # "data/timepoint2_data.csv"

time1_raw <- readBaseData("data/timepoint1.csv", 1)
time2_raw <- readBaseData("data/timepoint2.csv", 2)

collected_data <- read.csv(file = "outputs/summary/TP1_cluster_results.csv", 
                           stringsAsFactors = FALSE, numerals = "no.loss") %>% 
  as_tibble() %>% select(2, 3, 4, 8, 9, 10, 11, 12, 13) %>% 
  set_colnames(c("tp1_h", "tp1_cl", "tp1_cl_size", "tp2_h", "tp2_cl", 
                 "tp2_cl_size", "add_TP1", "num_nov", "size_change")) %>% 
  arrange(tp1_cl, tp2_h, tp2_cl)

collected_data$tp1_h %<>% charToInt(., "h")
collected_data$tp1_cl %<>% charToInt(., "c")
collected_data$tp2_h %<>% charToInt(., "h")
collected_data$tp2_cl %<>% charToInt(., "c")

heights <- collected_data$tp1_h %>% unique()

b1 <- meltedSizing(time1_raw, 1)
b2 <- meltedSizing(time2_raw, 2)
novels <- setdiff(b2$isolate, b1$isolate)

for (j in 1:length(heights)) {
  cd <- collected_data %>% filter(tp1_h == heights[j])
  
  hx <- b1 %>% filter(tp1_h == heights[j])
  b1clusters <- hx$tp1_cl %>% unique()
  
  pb <- txtProgressBar(min = 0, max = length(b1clusters), initial = 0, style = 3)
  actual_data <- lapply(1:length(b1clusters), function(i) {
    setTxtProgressBar(pb, i)
    c1 <- hx %>% filter(tp1_cl == b1clusters[i])
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
  }) %>% bind_rows() %>% arrange(tp1_cl, tp2_h, tp2_cl) %>% 
    mutate(across(tp2_cl_size, as.integer))
  close(pb)
  
  ad <- actual_data %>% select(tp1_h, tp1_cl, tp1_cl_size, tp2_h, 
                               tp2_cl, tp2_cl_size, add_TP1, num_nov) %>% 
    add_column(size_change = actual_data$tp2_cl_size - actual_data$tp1_cl_size) %>% 
    arrange(tp1_cl, tp2_h, tp2_cl)
  
  print(paste0("Threshold ", heights[j], " (", j, " of ", length(heights), "): ", identical(cd, ad)))
}


