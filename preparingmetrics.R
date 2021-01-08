#! /usr/bin/env Rscript
input_args = commandArgs(trailingOnly = TRUE)

# sending errors and warnings to the log file
msg <- file("outputs/logfile_prepmetrics.txt", open="wt")
sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("functions/base_functions.R")))
source("functions/processing_functions.R")
source("model_pred.R")

stopwatch <- rep(0,6) %>% set_names(c("start_time", "t2", "t3", "t4", "end_time"))
stopwatch[1] <- Sys.time()

sb <- list(start = paste0(rep("=", 13), collapse = ""), end = paste0(rep("=", 34), collapse = ""))

outputDetails(paste0("\n||", sb$start, " Summary file of results - compiling cluster metrics for all heights ", 
                     sb$start, "||", "\nStarted process at: ", Sys.time(), collapse = ""))

cat(paste0("\nFor process details, please see the log file at: outputs/logfile_prepmetrics.txt.\n", 
           "If at any point the process cuts off abruptly with no success message, see the log.\n"))

outputDetails(paste0("\nPART 1 OF 4: Data processing ", paste0(rep(".", 70), collapse = ""), sep = ""))
message("Preparing data for processing; loading datafiles")

## FOR USER: replace the filename variables with quoted file paths if you don't want to input them each time
# the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||
# the raw datasets, no filtering or other changes made
time1_raw <- input_args[1] %>% readBaseData(., 1)
time2_raw <- input_args[2] %>% readBaseData(., 2)
message("Successfully read in datafiles")

allheights <- colnames(time1_raw)[-1]

# this part takes 15 minutes for the 910 heights
outputDetails(paste0("\nPART 2 OF 4: Adding growth columns ", 
                     paste0(rep(".", 64), collapse = ""), sep = ""), newcat = TRUE)

stopwatch[2] <- Sys.time()
message("Adding growth columns to the height data\nCheck the 'outputs/with_growth_rate/' directory for saved data")
# actual_growth_rate = (change in cluster size) / (original cluster size)
# b_ov_growth = (number of novel isolates) / (growth rate)
m <- length(allheights)
p2 <- txtProgressBar(min = 0, max = m, initial = 0, style = 3)
a <- lapply(1:m, function(i) {
  setTxtProgressBar(p2, i)
  df <- accData(allheights[i])
  saveData(dtype = 4, oh = df, h = allheights[i])
})
close(p2)

# takes ten minutes to run the following: 
outputDetails(paste0("\nPART 3 OF 4: Identifying most growth for each cluster ", 
                     paste0(rep(".", 45), collapse = ""), sep = ""), newcat = TRUE)

stopwatch[3] <- Sys.time()
message(paste0("Identifying the TP2 cluster with the 'most' growth for each TP1 cluster, saved for ", 
               "each threshold\nCheck the 'outputs/best_growth_per_clust/' directory for saved data"))

p3 <- txtProgressBar(min = 0, max = length(allheights), initial = 0, style = 3)
for (i in 1:length(allheights)) {
  h <- allheights[i]
  setTxtProgressBar(p3, i)
  df1 <- readData(dtype = 4, h) %>%
    arrange(-b_ov_growth, tp2_h, tp2_cl) %>%
    group_by(id) %>% slice(1) %>% ungroup() %>%
    arrange(tp1_h, tp1_cl)
  saveData(dtype = 5, df = df1, h = h)
}
close(p3)

outputDetails(paste0("\nPART 4 OF 4: Comparing change in cluster size to the predicted model ", 
                     paste0(rep(".", 30), collapse = ""), sep = ""), newcat = TRUE)

stopwatch[4] <- Sys.time()

p4 <- txtProgressBar(min = 0, max = length(allheights), initial = 0, style = 3)
saving_results <- lapply(1:length(allheights), function(i) {
  h <- allheights[i]
  setTxtProgressBar(p4, i)
  a1 <- readData(dtype = 5, h = h)
  predicted_y <- predict(model_used, newdata = a1$tp1_cl_size)
  na_predicted <- a1$tp1_cl_size[which(is.na(predicted_y))]

  if (all(na_predicted > max(lm_df$`Cluster size`))) {
    predicted_y[is.na(predicted_y)] <- 15
  }

  # FOLD CHANGE
  # Ranking the data by the actual proportional change over the adaptive threshold
  # requirement, to see by how much each cluster exceeds the growth prediction.
  # The data is then saved to the outputs folder in the current working directory.
  # message("Ranking the data by the actual proportional change over the adaptive threshold requirement.")
  transit <- predicted_y %>% round() %>% bind_cols(a1, predicted = .)
  transit$predicted <- transit$predicted*0.01

  ## Fold change columns, convert decimals to percent
  transit$fold_change <- (transit$actual_growth_rate / transit$predicted) %>% round(., digits = 3)
  transit %>% arrange(., -fold_change) %>% 
    set_colnames(c("Height (TP1)", "Cluster (TP1)", "TP1 cluster size", 
                   "Height (TP2)", "Cluster (TP2)", "TP2 cluster size", "Number of novels", 
                   "Number of identical clusters before this one (at first height)", "ID", 
                   "Number of heights through which this cluster has existed", "Actual growth rate", 
                   "b_Over_Growth", "Adaptive threshold", "Fold change (Growth / Threshold)")) %>% return()
}) %>% bind_rows()
close(p4)

saveData(dtype = 6, res = saving_results)
stopwatch[5] <- Sys.time()
saveData(dtype = 7, sw = stopwatch)

outputDetails(paste0("\nFinished process at: ", Sys.time(), "\n", 
                     "||", sb$end, " End of collecting results ", sb$end, "||\n", sep = ""), newcat = TRUE)

