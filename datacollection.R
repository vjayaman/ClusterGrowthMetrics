#! /usr/bin/env Rscript
input_args = commandArgs(trailingOnly = TRUE)

# SUMMARY: this is the script for generating the cluster metrics, given input files for two time points
# 	- each time point dataset should be a file that can be read with the following statement 
# 	      read.csv(file = tp1_filename, stringsAsFactors = FALSE, numerals = "no.loss")
# 	  and should be composed of a single column with name "isolate" followed by a series of columns with 
# 	  cluster assignments at each of the thresholds
# 	- it may take some time if either time point dataset is very very large (>> 6000 genomes)

# sending errors and warnings to the log file
# msg <- file("logfile_datacollection.txt", open="wt")
# sink(msg, type="message")

stopwatch <- rep(0,2)
stopwatch[1] <- Sys.time()

cat(paste0("\n||-------------------------------- Cluster metric generation --------------------------------||\n"))

suppressWarnings(suppressPackageStartupMessages(source("functions/base_functions.R")))
source("functions/processing_functions.R")

cat(paste0("\nIf at any point the process cuts off abruptly with no success message, please see the log file.\n"))

cat(paste0("\nPreparing data for processing (see log file for details)\n"))

# FOR USER: replace the filename variables with quoted file paths if you don't want to input them each time
# the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||
message("\nLoading datafiles...")
# the coded datasets (where the isolates are now replaced by positive integers)
# t1_coded <- "data/timepoint1_data.csv" %>% readData(., 1) %>% codeIsolates(., "tp1")
# t2_coded <- "data/timepoint2_data.csv" %>% readData(., 2) %>% codeIsolates(., "tp2")
t1_coded <- input_args[1] %>% readData(., 1) %>% codeIsolates(., "tp1")
t2_coded <- input_args[2] %>% readData(., 2) %>% codeIsolates(., "tp2")

message("Successfully read in datafiles")

# Note: assumption made that first column is labelled "isolate"
message("Processing cluster composition for easier data handling - time point 1 (TP1)")
t1_comps <- t1_coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>% compsSet(., 1)

message("Collecting and counting novel isolates in each TP2 cluster...")
novels <- setdiff(t2_coded$isolate, t1_coded$isolate)
counting_novels <- t2_coded %>%
  filter(isolate %in% novels) %>%
  group_by(id) %>%
  summarise(num_novs = n(), .groups = "drop")

message("Processing cluster composition for easier data handling - time point 2 (TP2)")
t2_comps <- t2_coded %>% 
  rename(tp_h = tp2_h, tp_cl = tp2_cl) %>% 
  compsSet(., 2) %>% 
  left_join(., counting_novels, by = "id")
t2_comps$num_novs[is.na(t2_comps$num_novs)] <- 0

ids <- list()

### BASE CASE:
message("Collecting height data for base case, height 0:")
# this should be '0', the first column is the isolates

h_before <- unique(t1_coded$tp1_h)[1]
hdata <- t1_comps %>% filter(tp1_h == h_before) %>% arrange(tp1_h, tp1_cl)
message("\nTracking clusters")
cc <- trackClusters(hdata, t2_comps)

# note: the way these are flagged, we can see how many TP2 clusters were formed from the 
# selected TP1 cluster before we get to the "key growth" where a cluster went from its state
# at TP1 to its state at TP2 with a remarkable amount of growth
message("Flagging clusters")
clx <- unique(cc$tp1_cl)
fcb <- txtProgressBar(min = 0, max = length(clx), initial = 0, style = 3)
hbdata <- lapply(1:length(clx), function(i) {
  setTxtProgressBar(fcb, i)
  cc %>% filter(tp1_cl == clx[i]) %>% 
    arrange(tp2_h, tp2_cl) %>% 
    mutate(flag = 1:nrow(.)) %>% return()
}) %>% bind_rows() %>% 
  createID(., "tp1_h", "tp1_cl") %>% 
  add_column(flagged_heights = 0)
close(fcb)
saveData(dtype=3, tmp=hbdata, h=0)

bef_comps <- t1_comps %>% filter(tp1_h == h_before) %>%
  set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))

stopwatch <- rep(0,2) %>% set_names(c("start_time", "end_time"))
stopwatch[1] <- Sys.time()

heights <- unique(t1_comps$tp1_h)
cat(paste0("\nCollecting data for the rest (", length(heights) - 1, ") of the heights. ",  
           "This may take some time. \n(Note that for a more detailed look at progress, ", 
           "you can keep an eye on the newly created outputs directory, \nwhere updated ", 
           "data will be saved after each height.)\n"))

# ALL OTHER HEIGHTS
for (h_after in heights[-1]) {
  cat(paste0("\nCollecting data for height ", h_after, "/", length(heights), ":"))
  
  aft_comps <- t1_comps %>% filter(tp1_h == h_after) %>% 
    set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft"))
  
  # Part 1
  ss <- noChange(aft_comps, bef_comps, hbdata)
  
  # Part 2
  hdata <- aft_comps %>% 
    filter(!(id_aft %in% ss$id)) %>% 
    set_colnames(colnames(t1_comps))
  
  message("\nTracking clusters")
  cc <- trackClusters(hdata, t2_comps) %>% add_column(flag = NA)
  
  # Part 3
  message("Flagging clusters...")
  if (nrow(cc) > 0) {
    
    clx <- unique(cc$tp1_cl)

    fcb <- txtProgressBar(min = 0, max = length(clx), initial = 0, style = 3)
    fb <- lapply(1:length(clx), function(i) {
      setTxtProgressBar(fcb, i)
      cc %>% filter(tp1_cl == clx[i]) %>% 
        arrange(tp2_h, tp2_cl) %>% 
        mutate(flag = 1:nrow(.)) %>% return()
    }) %>% bind_rows()
    close(fcb)
    
    hnew <- fb %>% createID(., "tp1_h", "tp1_cl") %>% 
      add_column(flagged_heights = 0) %>% 
      bind_rows(., ss) %>% 
      arrange(tp1_h, tp1_cl)
    
  }else {
    hnew <- ss %>% arrange(tp1_h, tp1_cl)
  }
  saveData(dtype=3, tmp=hnew, h=h_after)
  
  # Part 4 - prepping data for next iteration
  hbdata <- hnew
  h_before <- h_after
  bef_comps <- aft_comps %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
  stopwatch[2] <- Sys.time()
}
# restarted process at 4:18
# restarted process at 8:22, starting height 5 at 8:42, height 11 at 8:47, height 428 at 9:18, height 630 at 9:29, 
# height 893 at 9:50, finished at 9:53
# AN HOUR AND A HALF!!!

# saveData(dtype = 2, sw = stopwatch)
# message("\nMerging height data into a single results file - note, just tracked cluster info.\n")
# saveData(dtype = 3)
# 
# timeTaken(pt = "data collection", stopwatch)
# message("Closing all connections")
# closeAllConnections()
# 
# cat(paste0("\nSuccessfully collected data for all heights. (Next step: run preparingmetrics.R)"))
# 