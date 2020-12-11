saveData <- function(dtype = 1, h = NULL, sh = NULL, sw = NULL, m = NULL, transit, flagged) {
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
  }else if (dtype == 5) {
    saveRDS(transit[['i']], paste0("outputs/transit/ids/h_", h, ".Rds"))
    saveRDS(transit[['po']], paste0("outputs/transit/postcc//h_", h, ".Rds"))
    saveRDS(transit[['pr']], paste0("outputs/transit/precc//h_", h, ".Rds"))
    saveRDS(transit[['m']], paste0("outputs/transit/metrics/h_", h, ".Rds"))
  }else if (dtype == 6) {
    write.csv(flagged, "outputs/flagged_cases_w_isolates.csv", row.names = FALSE)
  }else if (dtype == 7) {
    write.csv(flagged, "outputs/flagged_cases_just_clusters.csv", row.names = FALSE)
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
