
source("flagging_functions.R")
pb <- progress_bar$new(total = 22)
pb$tick()
# the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||
# the raw datasets, no filtering or other changes made
time1_raw <- "data/timepoint1_data.csv" %>% readData(., X)
time2_raw <- "data/timepoint2_data.csv" %>% readData(., Y)

# USER: make sure the first column is the isolate labeling
# then replace with "isolate" for easier manipulation later on
colnames(time1_raw)[1] <- colnames(time2_raw)[1] <- "isolate"
pb$tick()
# isolates found at TP2 (both original and novel isolates) - found at TP1 (the original isolates):
# isolates introduced at TP2 (novel isolates):
novels <- setdiff(time2_raw$isolate, time1_raw$isolate)
pb$tick()
# tracked_clusters is the list of original genomes and their cluster assignments at all heights as well 
# as their originating clusters, with clusters flagged if their originating cluster has been found before
# the clusters that actually changed from TP1 to TP2 are flagged with an NA
original_tracking <- readRDS("outputs2/results.Rds") %>% 
  createID(., "tp2_h", "tp2_cl")
pb$tick()
# in this dataset we have, for all ORIGINAL isolates, the cluster assignments at TP1 and TP2, 
# the change in cluster size, and whether or not the originating cluster was found before

# the next step:
#   - for each novel isolate
#     - get the cluster assignments at all heights (for TP2)
#       - then match with the originating cluster at TP1 (if one existed)
#       - use the originating cluster set (so the flags are included)
#     - get the cluster sizes (TP1 and TP2)
# then:
#   - proportional increase, adaptive threshold, fold change

# input: time2_raw, novels, result data (original_tracking)
# This part is just to check that all the TP2 clusters are represented
checkForMissing(time2_raw, original_tracking, chktype = 1, novels)
pb$tick(len = 3)

tracked <- time2_raw %>% addNovelsToResults(., novels) %>% add_column(prop_inc = NA)
# CHECKING: the following should be true, to show that we didn't lose any clusters on the way so far
# identical(sort(unique(all_sizes$id)), sort(unique(tracked$id)))

# HANDLE THESE CASES SEPARATELY: 
zero_denoms <- tracked %>% filter(tp1_cl_size == 0)
pb$tick()
# we now are only going to look at clusters that existed at TP1 first
# we'll add the zero_denominator cases back in later
tracked <- tracked %>% filter(tp1_cl_size > 0)
pb$tick()
tracked$prop_inc <- (tracked$tp2_cl_size - tracked$tp1_cl_size) / tracked$tp1_cl_size
pb$tick()
## Adaptive threshold development

# We can plot what a couple of different threshold models can look like.
# This part will be heavily altered depending on what the actual data looks like.
# Currently being used to model synthetic data, and what we might expect significant
# size change to look like.

# message("Running local regression steps for adaptive threshold development")
x <- c(1, 25, 50, 100, 150)
y <- c(300, 150, 75, 30, 15)
final_df <- addPredictedToResults(x, y, tracked)
pb$tick()
## Fold change columns, convert decimals to percent
final_df$fold_change <- (final_df$prop_inc / final_df$predicted) %>% round(., digits = 3)
pb$tick()
final_df <- final_df %>% 
  arrange(., -fold_change) %>% 
  mutate(across(c(predicted, prop_inc), scales::percent))
pb$tick()
# at this point, the columns we have data for are:
#   - isolate
#   - tp1_h, tp1_cl, tp1_cl_size, 
#   - tp2_h, tp2_cl, tp2_cl_size, 
#   - flagged, id, prop_inc, predicted, fold_change
# We still need a column with the number of novels in each TP2 cluster (note that this may not correspond 
# exactly with the change in a cluster's size from TP1 to TP2)
pb$tick()
nov_sizes <- time2_raw %>% 
  filter(isolate %in% novels) %>% 
  clusterIDS() %>% 
  select(-isolate, -id) %>%
  group_by_all() %>% 
  summarise(num_novels = n(), .groups = "drop")
pb$tick(len = 2)
zero_denoms$fold_change <- zero_denoms$predicted <- NA

final_df <- final_df %>% 
  left_join(., nov_sizes, by = c("tp2_h", "tp2_cl")) %>% 
  bind_rows(., zero_denoms)
pb$tick()
checkForMissing(time2_raw, final_df, chktype = 2)
pb$tick(len = 3)
metrics <- final_df %>% select(isolate, tp1_h, tp1_cl, tp2_h, tp2_cl, tp1_cl_size, tp2_cl_size, 
                               flagged, num_novels, prop_inc, predicted, fold_change) %>% 
  set_colnames(c("Isolate", "Height (TP1)", "Cluster (TP1)", "Height (TP2)", "Cluster (TP2)", 
                 "TP1 cluster size", "TP2 cluster size", "Flag (originator seen before)", 
                 "Number of novels", "Proportional growth", "Adaptive threshold", 
                 "Fold change (Growth / Threshold)"))
pb$tick()
saveData(dtype = 4, m = metrics)
pb$tick()