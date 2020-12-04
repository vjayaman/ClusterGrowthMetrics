
source("flagging_functions.R")
testing <- TRUE

pb <- progress_bar$new(total = 22)
pb$tick()

# the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||
# the raw datasets, no filtering or other changes made
time1_raw <- "data/timepoint1_data.csv" %>% readData(., X)
time2_raw <- "data/timepoint2_data.csv" %>% readData(., Y)

# USER: make sure the first column is the isolate labeling
# then replace with "isolate" for easier manipulation later on
colnames(time1_raw)[1] <- colnames(time2_raw)[1] <- "isolate"

# isolates found at TP2 (both original and novel isolates) - found at TP1 (the original isolates):
# isolates introduced at TP2 (novel isolates):
novels <- setdiff(time2_raw$isolate, time1_raw$isolate)
pb$tick()

# tracked_clusters is the list of original genomes and their cluster assignments at all heights as well 
# as their originating clusters, with clusters flagged if their originating cluster has been found before
# the clusters that actually changed from TP1 to TP2 are flagged with an NA

# checking the clusters of the results from datacollection.R
original_tracking <- readRDS("outputs/results.Rds") %>% createID(., "tp2_h", "tp2_cl")
pb$tick()

every_case <- time2_raw %>% 
  melt(., id = "isolate") %>% as_tibble() %>% 
  set_colnames(c("isolate","tp2_h", "tp2_cl")) %>% 
  arrange(isolate, tp2_h, tp2_cl) %>% 
  mutate(across(tp2_h, as.character)) %>% 
  mutate(across(tp2_h, as.integer)) %>% 
  createID(., "tp2_h", "tp2_cl")
pb$tick()
pb$tick()

if (testing) {
  check1a <- every_case %>%
    filter(!(isolate %in% novels))
  
  check1b <- original_tracking %>%
    select(isolate, tp2_h, tp2_cl) %>%
    arrange(isolate, tp2_h, tp2_cl) %>%
    createID(., "tp2_h", "tp2_cl")
  
  # TRUE --> in the results file, we have all the TP2 cluster assignments for the original isolates as well 
  # as the metrics of their originating clusters and flags for when they first emerged from a TP1 cluster
  if (identical(check1a, check1b)) {
    message("Confirming that the output of datacollection.R covers all clusters that contain original isolates.")
  }else {
    message("Problem: There are some clusters with original isolates that were not tracked.")
  }
}
pb$tick()

# WHAT WE NEED TO ADD:
# information regarding novels
#   a) clusters that contain both originals and novels (those that absorbed novels at TP2)
#     - we should have the origin information for these from the results file, we just need
#     to merge this data with that of the novels absorbed into their respective clusters
#   b) clusters that are entirely composed of novels at TP2 (i.e. did not exist at TP1)

# PART 1: novel only TP2 clusters
nov_cl <- every_case %>% filter(isolate %in% novels)
orig_cl <- every_case %>% filter(!(isolate %in% novels))
both <- every_case %>% filter(id %in% intersect(nov_cl$id, orig_cl$id))
pb$tick()

if (testing) {
  if (all(both$id %in% original_tracking$id)) {
    message("Confirming that all the clusters that contain both novels and originals have been tracked.")
  }else {
    message("Problem: There are clusters with original and novel isolates that were not tracked.")
  }
}
pb$tick()

# removing the originals, since we already have that information in original_tracking
both_n <- both %>% filter(isolate %in% novels)
oan <- original_tracking %>% select(-isolate) %>% 
  filter(id %in% both$id) %>% unique() %>% 
  left_join(both_n, ., by = c("tp2_h", "tp2_cl", "id"))

metrics <- original_tracking %>% bind_rows(., oan)
pb$tick()
pb$tick()
# ----------------------------------------------------------------------------------------------------
nov_only <- every_case %>% 
  filter(id %in% nov_cl$id) %>% 
  filter(!(id %in% orig_cl$id))
pb$tick()

if (testing) {
  # checking that these are clusters we haven't dealt with yet
  if (!any(nov_only$id %in% metrics$id)) {
    message("Confirming that the novel-only clusters have not been dealt with yet.")
  }else {
    message("Problem: Some of the novel-only clusters have been tracked (somehow) so far.")
  }
  # checking that all the isolates are not present in TP1
  if (any(nov_only$isolate %in% time1_raw$isolate)) {
    message("Problem: Some of the isolates were mislabelled as novel.")
  }else {
    message("Confirming that the list of novels we have at this point is accurate.")
  }
}
pb$tick()

nov_only <- nov_only %>% group_by(tp2_h, tp2_cl) %>% 
  summarise(tp2_cl_size = n(), .groups = "drop") %>% 
  left_join(nov_only, ., by = c("tp2_h", "tp2_cl")) %>% 
  arrange(tp2_h, tp2_cl)

nov_only <- nov_only %>% 
  add_column(tp1_h = NA, tp1_cl = NA, tp1_cl_size = 0, flagged = NA) %>% 
  select(colnames(metrics))

if (testing) {
  # checking that these clusters have not been accounted for
  if (any(nov_only$id %in% metrics$id)) {
    message("Problem: Some of the novel-only clusters have been tracked already.")
  }else {
    message("Confirming that the novel-only clusters have not been tracked yet.")
  }
}

metrics <- metrics %>% bind_rows(., nov_only)

if (testing) {
  # checking that all TP2 clusters have been now accounted before, and that there 
  # are no repeats or other anomalies
  a <- metrics$id %>% sort()
  b <- every_case$id %>% sort()
  if (identical(a, b)) {
    message("Confirming that all TP2 clusters have now been accounted for.")
  }else {
    message("Problem: There are TP2 clusters that have been overlooked, and are now not account for.")
  }
}

# TWO CASES FOR PROPORTIONAL INCREASE HANDLING
#   a) clusters did not exist at TP1 (denominator of zero)
#     --> prop_inc = NA
#   b) clusters did exist at TP1
#     --> prop_inc of (TP2 size - TP1 size) / (TP1 size)

metrics <- metrics %>% add_column(prop_inc = NA)

case2 <- which(metrics$tp1_cl_size != 0)
metrics$prop_inc[case2] <- metrics$tp2_cl_size[case2] %>% 
  '-'(metrics$tp1_cl_size[case2]) %>% '/'(metrics$tp1_cl_size[case2])

# ADAPTIVE THRESHOLD DEVELOPMENT

# We can plot what a couple of different threshold models can look like.
# This part will be heavily altered depending on what the actual data looks like.
# Currently being used to model synthetic data, and what we might expect significant
# size change to look like.

# message("Running local regression steps for adaptive threshold development")
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

# FOLD CHANGE
# # Ranking the data by the actual proportional change over the adaptive threshold requirement,
# # to see by how much each cluster exceeds the growth prediction. The data is then saved to
# # the current working directory.

final_df <- predicted_y %>% round() %>% bind_cols(metrics, predicted = .)
final_df$predicted <- final_df$predicted*0.01

## Fold change columns, convert decimals to percent
final_df$fold_change <- final_df$prop_inc %>% 
  '/'(final_df$predicted) %>% 
  round(., digits = 3)

final_df <- final_df %>% arrange(., -fold_change)

not_na <- which(!is.na(final_df$prop_inc))
final_df$prop_inc[not_na] <- final_df$prop_inc[not_na] %>% scales::percent() 
final_df$predicted[not_na] <- final_df$predicted[not_na] %>% scales::percent()

# We still need a column with the number of novels in each TP2 cluster (note that this may not correspond 
# exactly with the change in a cluster's size from TP1 to TP2)
nn <- every_case %>% 
  filter(isolate %in% novels) %>% 
  select(tp2_h, tp2_cl) %>% 
  group_by_all() %>% 
  summarise(number_novels = n(), .groups = "drop") %>% 
  createID(., "tp2_h", "tp2_cl")

results <- final_df %>% left_join(nn, by = c("tp2_h", "tp2_cl", "id")) %>% 
  select(isolate, tp1_h, tp1_cl, tp2_h, tp2_cl, tp1_cl_size, tp2_cl_size, 
         flagged, number_novels, prop_inc, predicted, fold_change) %>% 
  set_colnames(c("Isolate", "Height (TP1)", "Cluster (TP1)", "Height (TP2)", "Cluster (TP2)", 
                 "TP1 cluster size", "TP2 cluster size", "Flag (originator seen before)", 
                 "Number of novels", "Proportional growth", "Adaptive threshold", 
                 "Fold change (Growth / Threshold)"))

saveData(dtype = 4, m = results)
