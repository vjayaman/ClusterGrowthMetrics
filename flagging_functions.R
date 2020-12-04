x <- c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress")
lapply(x, require, character.only = TRUE)

# DATA HANDLING -----------------------------------------------------------------------------------------------------

get <- .Primitive("[[")

getUserInput <- function(msg) {
  cat(msg)
  op = readLines(con = "stdin", 1)
}

# given the defining filename in:
# C:\Users\vasen\Documents\Pre-MSc courses\Honours Research Project\SummerProject-2020\, read in the data
readData <- function(filename, file_number) {
  if (is.na(filename)) {
    stop(paste0("Time point ", file_number, " dataset not found."))
  }else {
    read.csv(file = filename, stringsAsFactors = FALSE, numerals = "no.loss", 
             check.names = FALSE) %>% as_tibble() %>% return()
  }
}

saveData <- function(dtype = 1, h = NULL, sh = NULL, sw = NULL, m = NULL) {
  if (dtype == 1) {
    paste0("outputs/height_data/h_", h, ".Rds") %>% saveRDS(sh, .)
  }else if (dtype == 2) {
    saveRDS(sw, "outputs/stopwatch.Rds")
    message("\nSaved stopwatch data.")
  }else if (dtype == 3) {
    mergeResults("outputs/height_data/") %>% saveRDS(., "outputs/results.Rds")
    message("\nSaved transitory results data - now need to run preparingmetrics.R")
  }else if (dtype == 4) {
    write.csv(m, "outputs/all_clusters_table.csv", row.names = FALSE)
    write.csv(m[1:10,], "outputs/first_ten_rows.csv", row.names = FALSE)
    message("\nSaved metrics to outputs folder.")
  }
}

meltData <- function(dataset, id_val) {
  melt(dataset, id = id_val) %>% as_tibble() %>% return()
}

factorToInt <- function(dataset, cname) {
  dataset %>% 
    mutate(across(all_of(cname), as.character)) %>% 
    mutate(across(all_of(cname), as.integer)) %>% return()
}

mergeResults <- function(data_dir) {
  hfiles <- list.files(data_dir)
  tracked_clusters <- paste0(data_dir, hfiles[1]) %>% readRDS()

  pb <- progress_bar$new(total = length(hfiles))
  for (h_file in hfiles[-1]) {
    pb$tick()
    next_file <- paste0(data_dir, h_file) %>% readRDS()
    tracked_clusters <- bind_rows(tracked_clusters, next_file)
  }
  tracked_clusters %>%
    mutate(across(tp2_h, as.integer)) %>% return()
}

# CLUSTER CHANGE PROCESSING -----------------------------------------------------------------------------------------
noChange <- function(postcomp, cc, precomp, single_height, ids) {
  # note that there are many clusters that exist at height1 but not at height2
  # (we are not interested in this representation) the focus is on clusters that now exist at height2
  indices_new <- cc %>% pull(composition)
  
  # the cluster assignments, before and after, for all the clusters that stayed the same 
  # *in composition* from the previous height
  stayed_the_same <- postcomp %>% 
    filter(!(composition %in% indices_new)) %>% 
    left_join(., precomp, by = "composition") %>% 
    select(-composition)
  
  # the result output for clusters that did not change in composition
  # since the clusters did not change at all, the "flags" should now become "sb"!
  transit <- inner_join(single_height, stayed_the_same, by = c("tp2_cl" = "h_before")) %>% 
    select(-tp2_cl) %>% rename(tp2_cl = h_after) %>% 
    select(colnames(single_height))
  transit$tp2_h <- height2
  
  # just checking:
  if (all(transit$flagged %in% ids)) {
    transit$flagged <- "sb"
  }else {
    stop("A cluster with a new originating cluster was misclassified. See noChange() function.")
  }

  return(transit)
}

# these are the clusters that change in composition from h0 to h1
# then for h1 at TP2, we only need to find the originating TP1 clusters for these ones, 
# all other clusters have the same originating cluster as the TP2 h0 clusters
changedClusters <- function(precomp, postcomp) {
  # cluster composition that is found at the new height that were different before
  indices_new <- setdiff(postcomp$composition, precomp$composition)
  # going to analyze the clusters that are actually different at height2
  changed_comp <- postcomp %>% filter(composition %in% indices_new)  
  # ac <- changed_comp$h_after
  return(changed_comp)
}
  
# we find the cluster assignments for all isolates at a particular height at TP2 
# that can be found in the given list of clusters (have not yet originating clusters for these)
clusterAssignments <- function(h, time2_data, key_clusters, novels) {
  # the cluster assignments for all clusters found in key_clusters, those that 
  # changed in composition from h1 to h2
  tmp <- time2_data %>% select(isolate, all_of(h)) %>% 
    set_colnames(c("isolates", "clusters")) %>% 
    dplyr::filter(clusters %in% key_clusters)
  
  # we group by cluster and then count the number of isolates found in each cluster
  # these are the actual cluster sizes
  tmp2 <- tmp %>% select(-isolates) %>% 
    countCases(., "tp2_cl_size") %>% 
    left_join(tmp, ., by = "clusters") %>% 
    set_colnames(c("isolate", "tp2_cl", "tp2_cl_size"))
  
  # We filter so it only includes the original isolates (but still has the actual TP2 cluster sizes)
  # then return a dataset with || isolate | height (at tp2) | cluster (at tp2) | cluster size (at tp2) ||
  # We will use this dataset to match with TP1 potential originating clusters (just_originals)
  tmp2 %>% filter(!(isolate %in% novels)) %>% add_column(tp2_h = h, .after = 1) %>% return()
}

createID <- function(df, c1, c2) {
  df %>% add_column(id = paste0(pull(df, c1), "-", pull(df, c2))) %>% return()
}

# THE PROBLEM: we have three cases
# (1) the originating cluster has not been seen before --> flag is given
# (2) the originating cluster has been seen before --> no flag is given 
#     --> should be clearer, give it a flag "seen before" or "sb"
# (3) the originating cluster does not exist --> no flag is given (this case was overlooked )
#     --> that way this remains NA
#     --> easier to understand the results of the entire analysis later on

checkID <- function(identifier, kc, x1, ids) {
  if (identifier %in% ids) {
    # if this originating cluster has been seen before, return an unflagged result
    list(new_ids = ids, 
         new_kc = kc %>% add_column(x1, .after = 1) %>% add_column(flagged = "sb")) %>% return()
  }else {
    # if this originating cluster has never been seen before, return a flagged result, with the ID
    list(new_ids = c(ids, identifier), 
         new_kc = kc %>% add_column(x1, .after = 1) %>% add_column(flagged = identifier)) %>% return()
  }
}

clusterIDS <- function(dataset, dtype = 2) {
  df <- dataset %>% 
    meltData(., "isolate") %>% 
    factorToInt("variable") %>% 
    set_colnames(c("isolate", "tp2_h", "tp2_cl")) %>% 
    createID(., "tp2_h", "tp2_cl")
  
  if (dtype == 1) {
    df %>% select(-isolate) %>% unique() %>% 
      arrange(tp2_h, tp2_cl) %>% return()
  }else {
    df %>% arrange(tp2_h, tp2_cl) %>% return()
  }
}

countCases <- function(dataset, last_col) {
  df <- dataset %>% 
    group_by_all() %>% 
    summarise(n = n(), .groups = "drop")
  colnames(df)[ncol(df)] <- last_col
  df %>% return()
}

# --------------------------------------------------------------------------------------------------------------
resultsProcess <- function(time2_data, all_clusters, df1, ids, jo) {
  # given the heights and cluster assignments for all TP1 isolates, we filter to keep only the genomes found 
  # in a particular cluster kc, then we identify all the clusters at TP1 that only contain these genomes
  # and then return the first height-and cluster pair where these isolates are found in a single cluster  
  
  # all the jo clusters should have at least one original isolate
  # note that the ac (all_clusters) list just indicates which clusters changed in composition
  #   - it is possible for one or more to contain only novels
  #     - these will be ignored for now, and dealt with in a separate step
  pb <- progress_bar$new(total = length(all_clusters))
  ids_copy <- ids
  a <- list()
  for (i in 1:length(all_clusters)) {
    pb$tick()
    # print(paste0("h_", h, "-", i, "/", length(all_clusters)))
    kc <- jo %>% dplyr::filter(tp2_cl == all_clusters[i])
    
    # this is true if the cluster actually has original isolates, and is not just made of novels
    # if it is just novels, we can skip it for now, it will be dealt with later
    if (nrow(kc) > 0) {
      # the first cluster in the TP1 dataset to contain only the originals found in the 
      # TP2 cluster, cluster_x, in form || tp1_h | tp1_cl | tp1_cl_size ||
      x1 <- df1 %>% 
        filter(isolate %in% kc$isolate) %>% 
        select(-isolate) %>%
        countCases(., "tp1_cl_size") %>% 
        filter(tp1_cl_size == nrow(kc)) %>% 
        slice(n = 1)
      
      # We add this "originating cluster" to the dataframe with this TP2 cluster.
      # If this "originating cluster" is already in the ID list, return with an unflagged notation, 
      # otherwise add to the list of IDs and flag the cluster.
      # NOTE: if we use an lapply, any IDs added to the list during this particular lapply 
      # are not compared to new IDs --> so if you add a new ID for a cluster, it can still be 
      # added for subsequent clusters for the same height.
      # 
      # This can only happen in cases where a cluster was larger at TP1, where these original isolates 
      # were actually clustered into two or more different clusters at TP1.
      # e.g. cluster 7-1310-TP1 (size 5), clustered into 1-1777-TP2 (size 3) and 1-1778-TP2 (size 2).
      # As this is just something that can happen with clustering, we will not make an exception in the 
      # flagging methodology --> will run as a for-loop so we can flag the secondary cases as "seen before", 
      # even if it was seen before at an earlier cluster for the same height
      # 
      # NOTE: need to update a copy of the ID list as well, so you can see when a cluster was flagged before, 
      # even if it was just at the same height
      verify_ID <- paste0(x1$tp1_h, "-", x1$tp1_cl) %>% checkID(., kc, x1, ids_copy)
      ids_copy <- verify_ID %>% '[['("new_ids")
      a[[i]] <- verify_ID %>% '[['("new_kc")  
    }
  }
  bind_rows(a) %>% return()
}

# these are all the clusters and their composition
clustComp <- function(df, height, dtype) {
  message(paste0("Preparing cluster composition IDs for height ", height))
  a1 <- df %>% 
    select(isolate, all_of(height)) %>% 
    set_colnames(c("isolate", "tp2_cl"))
  clusters <- a1$tp2_cl %>% unique()
  
  hx <- list("a1" = a1, "cl" = clusters)
  
  lapply(1:length(hx$cl), function(i) {
    filter(hx$a1, tp2_cl == hx$cl[i]) %>% 
      pull(isolate) %>% 
      paste0(collapse = ",") %>% return()
  }) %>% unlist() %>% 
    tibble(., hx$cl) %>% 
    set_colnames(c("composition", dtype)) %>% 
    return()
}

# cluster 7-1310 (is the TP1 originating cluster for both 1-1777 and 1-1778 at TP2)
# it is flagged for both, when it should only be flagged for 1-1777 --> regardless of the change in cluster 
# size, new clusters should only be flagged once per originating cluster

oneHeight <- function(df, df_coded, h, novel_isolates, meltedTP1, ids, ac, precomp) {
  # we find the cluster assignments for all isolates at a particular height at TP2 that can be 
  # found in the given list of clusters (have not yet found originating clusters for these)
  # dataset of form || isolate (originals) | tp2 height | tp2 cluster | tp2 cluster size ||
  just_originals <- clusterAssignments(h, df, ac, novel_isolates)
  # "just_originals" only includes the original isolates (but still has the actual TP2 cluster sizes) and 
  # is a dataset of form || isolate | height (at tp2) | cluster (at tp2) | cluster size (at tp2) ||
  pull(precomp, h_before) %>% collectionMsg(h, df, .)
  resultsProcess(df, ac, meltedTP1, ids, just_originals) %>% return()
}

# MESSAGES FOR USER -------------------------------------------------------------------------------------------------

collectionMsg <- function(h, time2_data, all_clusters) {
  if (match(h, colnames(time2_data)) == 2) { # if h == "0"
    message(paste0("Collecting data for the ", length(all_clusters), 
                   " clusters at height ", h, " (base case)"))
  }else {
    message(paste0("Collecting data for the ", length(all_clusters), " cluster(s) ", 
                   "that changed in composition from the previous height to the ", 
                   "current height, ", h))    
  }
}

# FUNCTIONS FOR TESTINGS --------------------------------------------------------------------------------------------

# This part is just to check that all the TP2 clusters are represented
checkForMissing <- function(df, results, chktype = 1, novel_isolates = NULL) {
  
  if (chktype == 1) {
    # definitely have originals and might have novels
    omn <- df %>% filter(!(isolate %in% novel_isolates)) %>% clusterIDS(., 1)
    omn_ids <- omn$id %>% unique() %>% sort()
    
    # definitely have novels and might have originals
    nmo <- df %>% filter(isolate %in% novel_isolates) %>% clusterIDS(., 1)
    nmo_ids <- nmo$id %>% unique() %>% sort()
    
    # novels only - present in the nmo set but not in the omn set
    novels_only <- setdiff(nmo_ids, omn_ids)
    
    # novels and originals - present in both the omn and nmo sets
    present_in_both <- intersect(omn_ids, nmo_ids)
    
    # originals only - present in the omn and not the nmo set
    originals_only <- setdiff(omn_ids, nmo_ids)
    
    in_tracked <- results %>%
      select(tp2_h, tp2_cl, id) %>%
      pull(id) %>% unique() %>% sort()
    x <- sum(length(present_in_both), length(originals_only)) == length(in_tracked)
  }else {
    actual_clusters <- df %>% clusterIDS(., 1)
    result_clusters <- results %>% 
      select(tp2_h, tp2_cl) %>% unique() %>% 
      arrange(tp2_h, tp2_cl) %>% 
      createID(., "tp2_h", "tp2_cl")
    
    # if the following is true, then all clusters have been accounted for
    x <- identical(actual_clusters$id, result_clusters$id)
  }
  
  # checking that the actual clusters match with the output of the tracking process
  if (isTRUE(x)) {
    paste0("\nNo missing clusters at stage ", chktype) %>% message() %>% return()
  }else {
    stop(message("\nSome of the clusters are missing - not found in the set of tracked clusters"))
  }
}


addNovelsToResults <- function(df, novel_isolates) {
  # get the TP2 cluster sizes
  tmp <- clusterIDS(df)
  
  all_sizes <- tmp %>% group_by(tp2_h, tp2_cl) %>% 
    summarise(tp2_cl_size = n(), .groups = "drop") %>% 
    left_join(tmp, ., by = c("tp2_h", "tp2_cl"))
  
  # novel cluster assignments
  cwn <- tmp %>% filter(isolate %in% novel_isolates)
  # non-novel cluster assignments
  cnn <- tmp %>% filter(!(isolate %in% novel_isolates))
  
  # the original_tracking dataset has only original isolates
  # this includes originals that are found in clusters with novels
  
  #   - and create a new set with only the novels
  nov_only <-  all_sizes %>% filter(id %in% setdiff(cwn$id, cnn$id)) %>% 
    add_column(tp1_h = NA, tp1_cl = NA, tp1_cl_size = 0, flagged = NA) %>% 
    select(isolate, tp1_h, tp1_cl, tp1_cl_size, tp2_h, tp2_cl, tp2_cl_size, flagged, id)
  
  original_tracking %>% bind_rows(., nov_only) %>% return()
}


addPredictedToResults <- function(x, y, tracked_cl) {
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
  predicted_y <- predict(model_used, newdata = tracked_cl$tp1_cl_size)
  na_predicted <- tracked_cl$tp1_cl_size[which(is.na(predicted_y))]
  
  if (all(na_predicted > max(lm_df$`Cluster size`))) {
    predicted_y[is.na(predicted_y)] <- 15
  }
  
  # ## Fold change
  # # Ranking the data by the actual proportional change over the adaptive threshold requirement,
  # # to see by how much each cluster exceeds the growth prediction. The data is then saved to
  # # the current working directory.
  
  final_df <- predicted_y %>% round() %>% bind_cols(tracked, predicted = .)
  final_df$predicted <- final_df$predicted*0.01
  return(final_df)
}


addToMetrics <- function(height, ids, metrics) {
  if (missing(metrics)) {
    # must be at height 0
    tibble(h = height, number_of_ids = length(ids)) %>% return()
  }else {
    metrics %>% add_row(tibble(h = height, number_of_ids = length(ids))) %>% return()
  }
}
