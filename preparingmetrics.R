#! /usr/bin/env Rscript
input_args = commandArgs(trailingOnly = TRUE)

# SUMMARY: this is the script for generating the cluster metrics, given input files for two time points
# 	- each time point dataset should be a file that can be read with the following statement 
# 	      read.csv(file = tp1_filename, stringsAsFactors = FALSE, numerals = "no.loss")
# 	  and should be composed of a single column with name "isolate" followed by a series of columns with 
# 	  cluster assignments at each of the thresholds
# 	- it may take some time if either time point dataset is very very large (>> 6000 genomes)

# sending errors and warnings to the log file
msg <- file("logfile_prepmetrics.txt", open="wt")
sink(msg, type="message")


stopwatch <- rep(0,2)
stopwatch[1] <- Sys.time()

message("\n--------- Result file generation - compilation of cluster metrics for all heights ---------\n")

source("functions/base_functions.R")
source("functions/processing_functions.R")

cat(paste0("\nIf at any point the process cuts off abruptly with no ", 
           "success message, please see the log file.\n"))

cat(paste0("\nPreparing data for processing (see log file for details)\n"))

pb <- txtProgressBar(min = 0, max = 20, initial = 0, style = 3)

## FOR USER: replace the filename variables with quoted file paths if you don't want to input them each time
message("Loading datafiles")
# the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||
# the raw datasets, no filtering or other changes made
time1_raw <- input_args[1] %>% readData(., 1)
time2_raw <- input_args[2] %>% readData(., 2)


message("Successfully read in datafiles")
setTxtProgressBar(pb, 1)

testing <- TRUE
stopwatch <- rep(0,2) %>% set_names(c("start_time", "end_time"))
stopwatch[1] <- Sys.time()

# USER: make sure the first column is the isolate labeling
# then replace with "isolate" for easier manipulation later on
colnames(time1_raw)[1] <- colnames(time2_raw)[1] <- "isolate"

# isolates found at TP2 (both original and novel isolates) - found at TP1 (the original isolates):
# isolates introduced at TP2 (novel isolates):
novels <- setdiff(time2_raw$isolate, time1_raw$isolate)

message("Reading in results of datacollection.R")

# tracked_clusters is the list of original genomes and their cluster assignments at all heights as well 
# as their originating clusters, with clusters flagged if their originating cluster has been found before
# the clusters that actually changed from TP1 to TP2 are flagged with an NA

# checking the clusters of the results from datacollection.R
original_tracking <- readRDS("outputs/results.Rds") %>% createID(., "tp2_h", "tp2_cl")

message("Checking the clusters tracked so far.")

every_case <- clusterIDS(time2_raw)
message("Testing...")
setTxtProgressBar(pb, 2)

if (testing) {
  check1a <- every_case %>% filter(!(isolate %in% novels)) %>% arrange(isolate, tp2_h, tp2_cl)
  
  message("Checking: that tracked clusters doesn't overlook anything.")
  
  check1b <- original_tracking %>% select(isolate, tp2_h, tp2_cl) %>%
    arrange(isolate, tp2_h, tp2_cl) %>% createID(., "tp2_h", "tp2_cl")
  
  # TRUE --> in the results file, we have all the TP2 cluster assignments for the original isolates as well 
  # as the metrics of their originating clusters and flags for when they first emerged from a TP1 cluster
  if (identical(check1a, check1b)) {
    message("\nChecked: that the output of datacollection.R covers all clusters that contain original isolates.")
  }else {
    message("Problem: There are some clusters with original isolates that were not tracked.")
  }
  setTxtProgressBar(pb, 3)
}
message("Identifying clusters with novels and original isolates.")

# WHAT WE NEED TO ADD:
# information regarding novels
#   a) clusters that contain both originals and novels (those that absorbed novels at TP2)
#     - we should have the origin information for these from the results file, we just need
#     to merge this data with that of the novels absorbed into their respective clusters
#   b) clusters that are entirely composed of novels at TP2 (i.e. did not exist at TP1)

# PART 1: novel only TP2 clusters
nov_cl <- every_case %>% filter(isolate %in% novels)
setTxtProgressBar(pb, 4)
orig_cl <- every_case %>% filter(!(isolate %in% novels))
both <- every_case %>% filter(id %in% intersect(nov_cl$id, orig_cl$id))

message("Testing...")
setTxtProgressBar(pb, 5)

if (testing) {
  if (all(both$id %in% original_tracking$id)) {
    message("Checked: that all the clusters that contain both novels and originals have been tracked.")
  }else {
    message("Problem: There are clusters with original and novel isolates that were not tracked.")
  }
}

# removing the originals, since we already have that information in original_tracking
both_n <- both %>% filter(isolate %in% novels)
setTxtProgressBar(pb, 6)

oan <- original_tracking %>% select(-isolate) %>% filter(id %in% both$id) %>% 
  unique() %>% left_join(both_n, ., by = c("tp2_h", "tp2_cl", "id"))
setTxtProgressBar(pb, 7)

metrics <- original_tracking %>% bind_rows(., oan)

message("Merging results for clusters with both novels and original isolates.")
setTxtProgressBar(pb, 8)

# ----------------------------------------------------------------------------------------------------
nov_only <- every_case %>% filter(id %in% nov_cl$id) %>% filter(!(id %in% orig_cl$id))

message("Checking the clusters with only novels.")
setTxtProgressBar(pb, 9)

nov_only <- nov_only %>% 
  select(-isolate, -id) %>% countCases("tp2_cl_size") %>% 
  left_join(nov_only, ., by = c("tp2_h", "tp2_cl")) %>% arrange(tp2_h, tp2_cl) %>% 
  add_column(tp1_h = NA, tp1_cl = NA, tp1_cl_size = 0, flagged = NA) %>% 
  select(colnames(metrics))

message("Testing...")
setTxtProgressBar(pb, 10)

if (testing) {
  if (!any(nov_only$id %in% metrics$id)) {
    message("Checked: that the novel-only clusters have not been dealt with yet.")
  }else {
    message("Problem: Some of the novel-only clusters have been tracked (somehow) so far.")
  }
  if (any(nov_only$isolate %in% time1_raw$isolate)) {
    message("Problem: Some of the isolates were mislabelled as novel.")
  }else {
    message("Checked: that the list of novels we have at this point is accurate.")
  }
}

metrics <- metrics %>% bind_rows(., nov_only)

if (testing) {
  a <- metrics$id %>% sort()
  b <- every_case$id %>% sort()
  if (identical(a, b)) {
    message("Checked: that all TP2 clusters have now been accounted for.")
  }else {
    message("Problem: There are TP2 clusters that have been overlooked, and are now not account for.")
  }
}

message("Adding the proportional increase column to the metrics set.")
setTxtProgressBar(pb, 11)

# TWO CASES FOR PROPORTIONAL INCREASE HANDLING
#   a) clusters did not exist at TP1 (denominator of zero)
#     --> prop_inc = NA
#   b) clusters did exist at TP1
#     --> prop_inc of (TP2 size - TP1 size) / (TP1 size)

metrics <- metrics %>% add_column(prop_inc = NA)

case2 <- which(metrics$tp1_cl_size != 0)
metrics$prop_inc[case2] <- metrics$tp2_cl_size[case2] %>% 
  '-'(metrics$tp1_cl_size[case2]) %>% '/'(metrics$tp1_cl_size[case2])

message("Going through the adaptive threshold development step.")
setTxtProgressBar(pb, 12)

# ADAPTIVE THRESHOLD DEVELOPMENT

# We can plot what a couple of different threshold models can look like.
# This part will be heavily altered depending on what the actual data looks like.
# Currently being used to model synthetic data, and what we might expect significant
# size change to look like.

message("Running local regression steps for adaptive threshold development")

x <- c(1, 25, 50, 100, 150)
y <- c(300, 150, 75, 30, 15)

# http://r-statistics.co/Loess-Regression-With-R.html
loessMod1 <- loess(y ~ x, span = 1)
smoothed1 <- predict(loessMod1)

loessMod5 <- loess(y ~ x, span = 5)
smoothed5 <- predict(loessMod5)

lm_df <- tibble(x, y, smoothed1, smoothed5) %>%
  set_colnames(c("x", "Preset values", "Local regression (span 1)", "Local regression (span 5)")) %>%
  meltData(., "x") %>% set_colnames(c("Cluster size", "Function", "Growth"))

message("Comparing models in the adaptive threshold step.")
setTxtProgressBar(pb, 13)

# The model selection step, and using it to determine what the adaptive threshold
# function values would be for each input of initial cluster size. Anything that
# needs to be extrapolated is set to a pre-determined plateau value of percent increase.

model_used <- loessMod1
# predict.lm(model_used, data.frame(x = 30)) # note there are different types of prediction methods
predicted_y <- predict(model_used, newdata = metrics$tp1_cl_size)
na_predicted <- metrics$tp1_cl_size[which(is.na(predicted_y))]

if (all(na_predicted > max(lm_df$`Cluster size`))) {
  predicted_y[is.na(predicted_y)] <- 15
}
message("Prep for ranking by fold change.")
setTxtProgressBar(pb, 14)

# FOLD CHANGE
# Ranking the data by the actual proportional change over the adaptive threshold requirement,
# to see by how much each cluster exceeds the growth prediction. The data is then saved to
# the outputs folder in the current working directory.
message("Ranking the data by the actual proportional change over the adaptive threshold requirement.")
transit <- predicted_y %>% round() %>% bind_cols(metrics, predicted = .)
transit$predicted <- transit$predicted*0.01

## Fold change columns, convert decimals to percent
transit$fold_change <- (transit$prop_inc / transit$predicted) %>% round(., digits = 3)

transit <- transit %>% arrange(., -fold_change)
message("Adding a column to show the number of novels in each TP2 cluster.")
setTxtProgressBar(pb, 15)

# not_na <- which(!is.na(transit$prop_inc))
# transit$prop_inc[not_na] <- transit$prop_inc[not_na] %>% scales::percent() 
# transit$predicted[not_na] <- transit$predicted[not_na] %>% scales::percent()

# We still need a column with the number of novels in each TP2 cluster (note that this may not correspond 
# exactly with the change in a cluster's size from TP1 to TP2)
nn <- every_case %>% filter(isolate %in% novels) %>% 
  select(tp2_h, tp2_cl) %>% countCases("number_novels") %>% 
  createID(., "tp2_h", "tp2_cl")
message("Formatting the metrics dataset for saving.")
setTxtProgressBar(pb, 16)

results <- transit %>% left_join(nn, by = c("tp2_h", "tp2_cl", "id")) %>% 
  select(isolate, tp1_h, tp1_cl, tp2_h, tp2_cl, tp1_cl_size, tp2_cl_size, 
         flagged, number_novels, prop_inc, predicted, fold_change) %>% 
  set_colnames(c("Isolate", "Height (TP1)", "Cluster (TP1)", "Height (TP2)", "Cluster (TP2)", 
                 "TP1 cluster size", "TP2 cluster size", "Flag (originator seen before)", 
                 "Number of novels", "Proportional growth", "Adaptive threshold", 
                 "Fold change (Growth / Threshold)"))
message("Saving data...")
setTxtProgressBar(pb, 17)

saveData(dtype = 4, m = results)

flagged_wi <- filter(metrics, !is.na(flagged))
saveData(dtype = 6, flagged = flagged_wi)
setTxtProgressBar(pb, 18)

flagged_jc <- filter(metrics, !is.na(flagged)) %>% select(-isolate) %>% unique()
saveData(dtype = 7, flagged = flagged_jc)
setTxtProgressBar(pb, 19)

stopwatch[2] <- Sys.time()
timeTaken(pt = "preparing metrics", stopwatch)
setTxtProgressBar(pb, 20)


