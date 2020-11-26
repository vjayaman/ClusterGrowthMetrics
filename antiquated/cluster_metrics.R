# # SUMMARY: this is the script for generating the cluster metrics, given input files for two time points
# # 	- each time point dataset should be a file that can be read with the following statement 
# # 	      read.csv(file = tp1_filename, stringsAsFactors = FALSE, numerals = "no.loss")
# # 	  and should be composed of a single column with name "isolate" followed by a series of columns with 
# # 	  cluster assignments at each of the thresholds
# # 	- it may take some time if either time point dataset is very very large (>> 6000 genomes)
#     - note that we only really use the cluster assignments of the TP2 dataset
#     - the TP1 data is used to identify which isolotes are the "original isolates" and not novel to TP2
# 

source("global.R")
source("tabulating_functions.R")

stopwatch <- rep(0,2)
stopwatch[1] <- Sys.time()

# ## FOR USER: replace the filename variables with quoted file paths if you don't want to input them each time
# message("Loading datafiles")
tp1_filename <- "data/timepoint1_data.csv" #input_args[1]      # "european-t1_clusters.csv"
tp2_filename <- "data/timepoint2_data.csv" #input_args[2]      # "european-t2_clusters.csv"

## Cluster metric generation
if (is.na(tp1_filename)) {
  stop(paste0("Time point 1 dataset not found."))
}else {
  time1 <- read.csv(file = tp1_filename, stringsAsFactors = FALSE, numerals = "no.loss", check.names = FALSE) %>% as_tibble()
}
X <- 1

if (is.na(tp2_filename)) {
  stop(paste0("Time point 2 dataset not found."))
}else {
  time2 <- read.csv(file = tp2_filename, stringsAsFactors = FALSE, numerals = "no.loss", check.names = FALSE) %>% as_tibble()
}
Y <- 2

# isolates found at TP1
isolates_tp1 <- time1$isolate

# isolates found at TP2
isolates_tp2 <- time2$isolate

# isolates introduced at TP2 (novel isolates)
new_isolates <- setdiff(isolates_tp2, isolates_tp1)

metrics_setup <- melt(time2, id = "isolate") %>% as_tibble() %>% 
  set_colnames(c("isolate", "height", "cluster"))
metrics_setup$height %<>% as.character()

# total cluster sizes for the TP2 dataset
sizes_for_tpY <- clusterSizes(time2) %>%
  set_colnames(c("height", "cluster", "TP2_cl_size", "ID"))
sizes_for_tpY$cluster %<>% as.character()

# general metrics for TP2 so far - sizes
# ordered by height and then column
tpY_data <- merge(metrics_setup, sizes_for_tpY) %>% as_tibble() %>% 
  dplyr::select(isolate, height, cluster, TP2_cl_size, ID)
tpY_data$height %<>% as.integer()
tpY_data <- tpY_data %>% dplyr::group_by(height, cluster)
# saveRDS(tpY_data, "tmpdata.Rds")

# multi-strain clusters at TP2 that contain novel isolates
novel_ms_clusters <- tpY_data %>% 
  dplyr::filter(TP2_cl_size > 1) %>% 
  dplyr::filter(isolate %in% new_isolates)

# multi-strain clusters at TP2 that contain original isolates
original_ms_clusters <- tpY_data %>% 
  dplyr::filter(TP2_cl_size > 1) %>% 
  dplyr::filter(!(isolate %in% new_isolates))

# NOTE: not looking at singletons at this point, will bring them back in later

# multi-strain clusters at TP2 that contain both novels and originals
# i.e. clusters that existed before the novels were added
# consider these KEY CLUSTERS
nmc_ids <- novel_ms_clusters %>% dplyr::select(-isolate) %>% unique()
omc_ids <- original_ms_clusters %>% dplyr::select(-isolate) %>% unique()
both_ms_clusters <- merge(nmc_ids, omc_ids) %>% as_tibble()

# counting the number of novel isolates in each key cluster
novel_sizes <- left_join(both_ms_clusters, novel_ms_clusters) %>% 
  dplyr::select(isolate, height, cluster, TP2_cl_size, ID) %>% 
  dplyr::group_by(height, cluster) %>% 
  dplyr::summarise(n = n()) %>% 
  set_colnames(c("height", "cluster", "num_novels"))

# counting the number of original isolates in each key cluster
original_sizes <- left_join(both_ms_clusters, original_ms_clusters) %>% 
  dplyr::select(isolate, height, cluster, TP2_cl_size, ID) %>% 
  dplyr::group_by(height, cluster) %>% 
  dplyr::summarise(n = n()) %>% 
  set_colnames(c("height", "cluster", "num_originals"))

# going to use actual TP2 size and the numbers of novels and original isolates in each cluster 
# to identify the change in theoretical cluster size
size_difs <- merge(novel_sizes, original_sizes, by = c("height", "cluster")) %>% as_tibble()
size_difs$height %<>% as.character()
size_difs$cluster %<>% as.character()
# || height | cluster | num_novels | num_originals | TP2_cl_size
# novel absorption (for all existing clusters that absorb at least one novel genome)
size_difs <- sizes_for_tpY %>% dplyr::select(-ID) %>% 
  merge(size_difs, ., by = c("height", "cluster")) %>% 
  as_tibble() %>% 
  dplyr::select(height, cluster, TP2_cl_size, num_originals, num_novels) %>% 
  set_colnames(c("TP2_h", "TP2_cl", "TP2_cl_size", "number_original", "number_novels"))
size_difs$ID <- paste0(size_difs$TP2_h, "-", size_difs$TP2_cl)

# the cluster assignments for each isolate
all_clusters <- tpY_data %>% set_colnames(c("isolate", "TP2_h", "TP2_cl", "TP2_cl_size", "ID")) %>% ungroup()
all_clusters$TP2_h %<>% as.character()
all_clusters$TP2_cl %<>% as.character()

# OUR METRICS FILE: will be adding to this one further
a1 <- left_join(all_clusters, size_difs)

# adding back in: singletons made of novel isolates
inds1 <- which(a1$TP2_cl_size == 1 & a1$isolate %in% new_isolates)
a1$number_original[inds1] <- 0
a1$number_novels[inds1] <- 1

# adding back in: singletons made of original isolates
inds2 <- which(a1$TP2_cl_size == 1 & !(a1$isolate %in% new_isolates))
a1$number_original[inds2] <- 1
a1$number_novels[inds2] <- 0

# adding back in: clusters that don't absorb any novels or are composed of only novels
inds3 <- a1[is.na(a1$number_original),] %>% pull(ID) %>% unique() %>% sort()
just_one_ids <- setdiff(unique(c(nmc_ids$ID, omc_ids$ID)), both_ms_clusters$ID) %>% sort()
test_for_leftovers <- identical(inds3, just_one_ids) # if TRUE, then we've accounted for all cases
# composed of only novels:
inds3a <- which(is.na(a1$number_original) & (a1$isolate %in% new_isolates))
a1$number_original[inds3a] <- 0
a1$number_novels[inds3a] <- a1$TP2_cl_size[inds3a]
# don't absorb any novels:
inds3b <- which(is.na(a1$number_original) & !(a1$isolate %in% new_isolates))
a1$number_original[inds3b] <- a1$TP2_cl_size[inds3b]
a1$number_novels[inds3b] <- 0

# size difference and proportional increase
a1$difference <- a1$number_novels
inds <- a1$number_original != 0
a1$prop_inc <- NA
a1$prop_inc[inds] <- a1$number_novels[inds] / a1$number_original[inds]
# so the rows with NA in the proportional increase column are those clusters that only contain novels
# WE REMOVE THESE FOR THE REST OF THIS PARTICULAR ANALYSIS - CAN SUB THEM BACK IN LATER
#   WANT TO AVOID DIVIDING BY ZERO
# a1 <- a1[!is.na(a1$prop_inc),] # any(is.na(a1$prop_inc) & a1$number_original != 0) --> FALSE
# ----------------------------------------------------------------------------------------------------------------------
## Adaptive threshold development

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
  melt(id = "x") %>% as_tibble() %>% set_colnames(c("Cluster size", "Function", "Growth"))

# The model selection step, and using it to determine what the adaptive threshold
# function values would be for each input of initial cluster size. Anything that
# needs to be extrapolated is set to a pre-determined plateau value of percent increase.

model_used <- loessMod1
# predict.lm(model_used, data.frame(x = 30)) # note there are different types of prediction methods
predicted_y <- predict(model_used, newdata = a1$number_original)
na_predicted <- a1$number_original[which(is.na(predicted_y))]

if (all(na_predicted > max(lm_df$`Cluster size`))) {
  predicted_y[is.na(predicted_y)] <- 15
}

sizes_predicted <- predicted_y %>% round() %>% bind_cols(a1, predicted = .)
sizes_predicted$predicted <- sizes_predicted$predicted*0.01

# ## Fold change
# # Ranking the data by the actual proportional change over the adaptive threshold requirement,
# # to see by how much each cluster exceeds the growth prediction. The data is then saved to
# # the current working directory.

sizes_predicted$fold_change <- sizes_predicted$prop_inc / sizes_predicted$predicted
sizes_predicted$fold_change <- sizes_predicted$fold_change %>% round(., digits = 3)
sizes_predicted <- sizes_predicted %>% arrange(., -fold_change)

sizes_predicted$predicted <- sizes_predicted$predicted %>% scales::percent()
sizes_predicted$prop_inc <- sizes_predicted$prop_inc %>% scales::percent()

metrics <- sizes_predicted %>% 
  dplyr::select(isolate, TP2_h, TP2_cl, number_original, TP2_cl_size, 
                number_novels, prop_inc, predicted, fold_change) %>% 
  set_colnames(c("Isolate", "Height (TP2)", "Cluster (TP2)", "TP1 cluster size", 
                 "TP2 cluster size", "Number of novels", "Proportional growth", 
                 "Adaptive threshold", "Fold change (Growth/Threshold)"))

write.csv(metrics, "outputs/all_clusters_table.csv", row.names = FALSE)
write.csv(metrics[1:10,], "outputs/cluster_metrics_first_ten.csv", row.names = FALSE)
stopwatch[2] <- Sys.time()