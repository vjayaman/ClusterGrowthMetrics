
source("processing_functions.R")

# the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||33333
# the raw datasets, no filtering or other changes made
time1_raw <- "data/timepoint1_data.csv" %>% readData(., X)
time2_raw <- "data/timepoint2_data.csv" %>% readData(., Y)

# USER: make sure the first column is the isolate labeling
# then replace with "isolate" for easier manipulation later on
colnames(time1_raw)[1] <- colnames(time2_raw)[1] <- "isolate"

# the coded TP2 datasets (where the isolates are now replaced by positive integers)
time2_coded <- time2_raw
time2_coded$isolate <- 1:nrow(time2_coded)

# isolates found at TP2 (both original and novel isolates) - found at TP1 (the original isolates):
# isolates introduced at TP2 (novel isolates):
novels <- setdiff(time2_raw$isolate, time1_raw$isolate)

# melted (long) format of the time point 1 dataset (note, contains assignments for all heights and isolates)
meltedTP1 <- time1_raw %>% 
  meltData(., "isolate") %>% 
  set_colnames(c("isolate", "tp1_h", "tp1_cl")) %>% 
  factorToInt(., "tp1_h")

ids <- list()

### BASE CASE:
# this should be '0', the first column is the isolates
base_case_h <- colnames(time2_raw)[2]
# finding the composition (in terms of coded isolates) of each of the clusters
# note that this includes both novel and original isolates
precc <- clustComp(time2_coded, base_case_h, "h_before")

single_height <- oneHeight(time2_raw, time2_coded, base_case_h, novels, meltedTP1, ids, precc$h_before, precc)
ids <- c(unlist(ids), single_height$flagged %>% unique()) %>% unique()
metrics <- addToMetrics(base_case_h, ids)

paste0("outputs/height_data/h_", base_case_h, ".Rds") %>% saveRDS(single_height, .)
saveData(dtype = 1, sh = single_height, h = "0")

stopwatch <- rep(0,2) %>% set_names(c("start_time", "end_time"))
stopwatch[1] <- Sys.time()

for (j in 1:length(colnames(time2_raw)[-1][-1])) {
  height2 <- colnames(time2_raw)[-1][-1][j]
  message(paste0("Threshold h_", height2, " - ", j, " / ", length(colnames(time2_raw)[-1][-1])))
  
  postcc <- clustComp(time2_coded, height2, "h_after")
  
  list('i'=ids, 'po'=postcc, 'pr'=precc, 'm'=metrics) %>% saveData(dtype = 5, h = height2, transit = .)
  
  # these are the clusters that changed from TP1 to TP2, and their new composition
  changed_comp <- changedClusters(precc, postcc)
  # these are the new clusters, for those that changed in composition
  ac <- changed_comp$h_after
  transit <- noChange(postcc, changed_comp, precc, single_height, ids)
  precc <- postcc %>% set_colnames(c("composition", "h_before"))
  
  # there is >= one cluster different in composition from previous height
  if (length(ac) > 0) {
    new_height <- oneHeight(time2_raw, time2_coded, height2, novels, meltedTP1, ids, ac, precc)
    
    # > 1 so "sb" doesn't trigger an error
    if (length(intersect(ids, new_height$flagged)) > 1) {
      stop(paste0("Some clusters at height ", height2, " are being flagged even though they've ", 
                  "been seen before!: ", paste0(intersect(ids, new_height$flagged), collapse = ",")))
    }
    ids <- c(ids, new_height$flagged) %>% unique()
    single_height <- new_height %>% bind_rows(transit, .)
    
  }else {
    message("No clusters changed from the previous height")
    single_height <- transit
  }
  
  metrics <- addToMetrics(height2, ids, metrics)
  
  saveData(dtype = 1, sh = single_height, h = height2)
  message("\n")
  stopwatch[2] <- Sys.time()
}
saveData(dtype = 2, sw = stopwatch)
saveData(dtype = 3)

timeTaken(pt = "data collection", stopwatch)
