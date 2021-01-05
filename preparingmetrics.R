#! /usr/bin/env Rscript
input_args = commandArgs(trailingOnly = TRUE)

# sending errors and warnings to the log file
# msg <- file("logfile_prepmetrics.txt", open="wt")
# sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("functions/base_functions.R")))
source("functions/processing_functions.R")

stopwatch <- rep(0,2) %>% set_names(c("start_time", "end_time"))
stopwatch[1] <- Sys.time()

cat(paste0("\n--------- Result file generation - compilation of cluster metrics ", 
           "for all heights ---------\n"))
cat(paste0("\nIf at any point the process cuts off abruptly with no ", 
           "success message, please see the log file.\n"))
cat(paste0("\nPreparing data for processing (see log file for details)\n"))

## FOR USER: replace the filename variables with quoted file paths if you don't want to input them each time
message("Loading datafiles")
# the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||
# the raw datasets, no filtering or other changes made
# time1_raw <- "data/timepoint1_data.csv" %>% readBaseData(., 1)
# time2_raw <- "data/timepoint2_data.csv" %>% readBaseData(., 1)
time1_raw <- input_args[1] %>% readBaseData(., 1)
time2_raw <- input_args[2] %>% readBaseData(., 2)
message("Successfully read in datafiles")


# dir.create("outputs/accdata")
# dir.create("outputs/bestacc")
accData <- function(h) {
  oneheight <- readRDS(paste0("outputs/height_data/h", h, ".Rds"))
  # growth = (change in cluster size) / (original cluster size)
  oneheight$actual_growth_rate <- (oneheight$tp2_cl_size - oneheight$tp1_cl_size) / oneheight$tp1_cl_size
  # growth acceleration = (number of novel isolates) / (growth rate)
  oneheight$acc <- oneheight$num_novs / oneheight$actual_growth_rate
  oneheight$acc[is.na(oneheight$acc)] <- 0
  saveRDS(oneheight, paste0("outputs/accdata/h", h, ".Rds"))
}

allheights <- colnames(time1_raw)[-1]

# this part takes 15 minutes for the 910 heights
m <- length(allheights)
pb <- txtProgressBar(min = 0, max = length(allheights), initial = 0, style = 3)
a <- lapply(1:m, function(i) {
  setTxtProgressBar(pb, i)
  accData(allheights[i])
})
close(pb)
Sys.time()


# takes ten minutes to run the following: 
for (i in 1:length(allheights)) {
  h <- allheights[i]
  print(h)
  readRDS(paste0("outputs/accdata/h", h, ".Rds")) %>% 
    arrange(-acc, tp2_h, tp2_cl) %>% 
    group_by(id) %>% slice(1) %>% ungroup() %>% 
    arrange(tp1_h, tp1_cl) %>% 
    saveRDS(., paste0("outputs/bestacc/h", h, ".Rds"))
}
Sys.time()


best_of_h <- lapply(1:m, function(i) {
  h <- allheights[i]
  a1 <- readRDS(paste0("outputs/bestacc/h", h, ".Rds")) %>% arrange(tp1_h, tp1_cl)
  a1[which.max(a1$acc),] %>% return()
}) %>% bind_rows()
# saveRDS(best_of_h, "bestgrowth.Rds")


source("model_pred.R")

saving_results <- lapply(allheights, function(h) {
  a1 <- readRDS(paste0("outputs/bestacc/h", h, ".Rds"))
  predicted_y <- predict(model_used, newdata = a1$tp1_cl_size)
  na_predicted <- a1$tp1_cl_size[which(is.na(predicted_y))]

  if (all(na_predicted > max(lm_df$`Cluster size`))) {
    predicted_y[is.na(predicted_y)] <- 15
  }
  # message("Prep for ranking by fold change.")

  # FOLD CHANGE
  # Ranking the data by the actual proportional change over the adaptive threshold 
  # requirement, to see by how much each cluster exceeds the growth prediction. 
  # The data is then saved to the outputs folder in the current working directory.
  # message("Ranking the data by the actual proportional change over the adaptive threshold requirement.")
  transit <- predicted_y %>% round() %>% bind_cols(a1, predicted = .)
  transit$predicted <- transit$predicted*0.01

  ## Fold change columns, convert decimals to percent
  transit$fold_change <- (transit$actual_growth_rate / transit$predicted) %>% round(., digits = 3)
  
  transit <- transit %>% arrange(., -fold_change)
  # message("Adding a column to show the number of novels in each TP2 cluster.")
  
  transit %>% set_colnames(c("Height (TP1)", "Cluster (TP1)", "TP1 cluster size", 
    "Height (TP2)", "Cluster (TP2)", "TP2 cluster size", "Number of novels", 
    "Number of identical clusters before this one (at first height)", "ID", 
    "Number of heights through which this cluster has existed", "Actual growth rate", 
    "Acc", "Adaptive threshold", "Fold change (Growth / Threshold)")) %>% return()
}) %>% bind_rows()

saveRDS(saving_results, "outputs/formatted_results.Rds")

