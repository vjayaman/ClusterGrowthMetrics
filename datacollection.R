#! /usr/bin/env Rscript
input_args = commandArgs(trailingOnly = TRUE)
# Note: a line of stars (********) will be placed after verified complicated operations, 
# involved statements that do what they are expected to do, and this has been checked for sure

# SUMMARY: this is the script for generating the cluster metrics, given input files for two time points
# 	- each time point dataset should be a file that can be read with the following statement 
# 	      read.csv(file = tp1_filename, stringsAsFactors = FALSE, numerals = "no.loss")
# 	  and should be composed of a single column with name "isolate" followed by a series of columns with 
# 	  cluster assignments at each of the thresholds

# sending errors and warnings to the log file
msg <- file("outputs/logfile_datacollection.txt", open="wt")
sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("funcs.R")))

stopwatch <- rep(0,2) %>% set_names(c("start_time", "end_time"))
stopwatch[1] <- Sys.time()

outputDetails(paste0("\n||", paste0(rep("-", 32), collapse = ""), " Cluster metric generation ", 
                     paste0(rep("-", 32), collapse = ""), "||\n", 
                     "Started process at: ", Sys.time()))

cat(paste0("\nIf at any point the process cuts off abruptly with no success message, ", 
           "please see the log file.\nPreparing data for processing.\n"))

# FOR USER: replace the filename variables with quoted file paths if you don't want to input them each time
# the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||

outputDetails(paste0("\nPART 1 OF 3: Data processing ", paste0(rep(".", 66), collapse = ""), "\n", 
                     "  Preparing data for processing\n  Loading and formatting datafiles..."), newcat = TRUE)

# the coded datasets (where the isolates are now replaced by positive integers)
# time1_raw <- "data/timepoint1_data.csv" %>% readBaseData(., 1)
# time2_raw <- "data/timepoint2_data.csv" %>% readBaseData(., 2)
time1_raw <- input_args[1] %>% readBaseData(., 1)
time2_raw <- input_args[2] %>% readBaseData(., 2)

all_isolates <- c(time1_raw$isolate, time2_raw$isolate) %>% unique() %>% 
  as_tibble() %>% set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

t1_coded <- time1_raw %>% codeIsolates(., "tp1", all_isolates)
t2_coded <- time2_raw %>% codeIsolates(., "tp2", all_isolates)

t2_colnames <- t2_coded$tp2_h %>% unique() %>% sort()
message("  Successfully read in datafiles")

# Note: assumption made that first column is labelled "isolate"
outputDetails("  Processing time point 1 (TP1) clusters for easier data handling", newcat = TRUE)
t1_comps <- t1_coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>%
  compsSet(., "TP1", indicate_progress = TRUE)

outputDetails("  Collecting and counting novel isolates in each TP2 cluster...", newcat = TRUE)
novels <- setdiff(t2_coded$isolate, t1_coded$isolate)
counting_novels <- t2_coded %>% filter(isolate %in% novels) %>%
  group_by(id) %>% summarise(num_novs = n(), .groups = "drop")

outputDetails("  Processing time point 2 (TP2) clusters for easier data handling...", newcat = TRUE)
t2_comps <- t2_coded %>% rename(tp_h = tp2_h, tp_cl = tp2_cl) %>%
  compsSet(., "TP2", indicate_progress = TRUE) %>% left_join(., counting_novels, by = "id")
t2_comps$num_novs[is.na(t2_comps$num_novs)] <- 0

### BASE CASE:
outputDetails(paste0("\nPART 2 OF 3: Tracking and flagging clusters for base case ", 
                     paste0(rep(".", 41), collapse = "")), newcat = TRUE)
outputDetails("  Collecting height data for base case, height 0...", newcat = TRUE)
# this should be '0', the first column is the isolates

outputDetails("  Tracking clusters", newcat = TRUE)
h_before <- unique(t1_coded$tp1_h)[1]
hdata <- t1_comps %>% filter(tp1_h == h_before) %>% arrange(tp1_h, tp1_cl)

cc <- trackClusters(hdata, t2_comps, t2_colnames, t1_coded, t2_coded, indicate_progress = TRUE) %>%
  createID(., "tp1", "tp1_h", "tp1_cl")

# note: the way these are flagged, we can see how many TP2 clusters were formed from the
# selected TP1 cluster before we get to the "key growth" where a cluster went from its state
# at TP1 to its state at TP2 with a remarkable amount of growth

outputDetails("  Flagging clusters", newcat = TRUE)
hbdata <- cc %>% add_column(flag = cc$id)

saveData(dtype=3, tmp=hbdata, h=0)

bef_comps <- t1_comps %>% filter(tp1_h == h_before) %>%
  set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))

heights <- unique(t1_comps$tp1_h)

outputDetails(paste0("\nPART 3 OF 3: Tracking and flagging clusters for the rest of the heights (", 
                     length(heights) - 1, " of them) ........."), newcat = TRUE)

outputDetails(paste0("  This may take some time. \n  ", 
                     "Note that for a more detailed look at progress, you can keep an eye on the outputs \n  ", 
                     "directory, where updated data will be saved after each height.\n"))

# ALL OTHER HEIGHTS
outputDetails("  Collecting data for other heights: ", newcat = TRUE)

fcb <- txtProgressBar(min = 0, max = length(heights[-1]), initial = 0, style = 3)
for (j in 1:length(heights[-1])) {
  setTxtProgressBar(fcb, j)
  
  h_after <- heights[-1][j]
  message(paste0("  Height ", j, " / ", length(heights)))
  
  aft_comps <- t1_comps %>% filter(tp1_h == h_after) %>%
    set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft"))

  # Part 1: identifying clusters that have not changed from the previous height
  ss <- noChange(aft_comps, bef_comps, hbdata)
  
  # Part 2: tracking clusters that changed
  hdata <- aft_comps %>%
    filter(!(id_aft %in% ss$id)) %>%
    set_colnames(colnames(t1_comps))
  
  # indicate_progress = FALSE
  cc <- trackClusters(hdata, t2_comps, t2_colnames, 
                      t1_coded, t2_coded, indicate_progress = FALSE) %>% 
    createID(., "tp1", "tp1_h", "tp1_cl")
  
  # Part 3: flagging clusters to indicate when they were first seen
  if (nrow(cc) > 0) {
    hnew <- cc %>% add_column(flag = cc$id) %>% bind_rows(., ss) %>% arrange(tp1_h, tp1_cl)
  }else {
    hnew <- ss %>% arrange(tp1_h, tp1_cl)
  }
  saveData(dtype=3, tmp=hnew, h=h_after)

  # Part 4 - prepping data for next iteration
  hbdata <- hnew
  h_before <- h_after
  bef_comps <- aft_comps %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
}
close(fcb)

stopwatch[2] <- Sys.time()
saveData(dtype = 7, sw = stopwatch)
outputDetails("  Merging height data into a single results file - note, just tracked clusters.", newcat = TRUE)

timeTaken(pt = "data collection", stopwatch) %>% outputDetails(., newcat = TRUE)
outputDetails(paste0("  Successfully collected data for all heights."), newcat = TRUE)
outputDetails(paste0("  Next step is to run the following: Rscript preparingmetrics.R tp1data.csv tp2data.csv"), newcat = TRUE)

outputDetails("  Closing all connections...", newcat = TRUE)
closeAllConnections()
