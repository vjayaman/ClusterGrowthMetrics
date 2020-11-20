x <- c("tibble", "magrittr", "dplyr", "reshape2", "scales")
lapply(x, require, character.only = TRUE)

get <- .Primitive("[[")

getUserInput <- function(msg) {
  cat(msg)
  op = readLines(con = "stdin", 1)
}

# given the defining filename in:
# C:\Users\vasen\Documents\Pre-MSc courses\Honours Research Project\SummerProject-2020\, 
# read in the data
readData <- function(filename, file_number) {
  
  if (is.na(filename)) {
    stop(paste0("Time point ", file_number, " dataset not found."))
  }else {
    time_X <- read.csv(file = filename, stringsAsFactors = FALSE, numerals = "no.loss", check.names = FALSE) %>% as_tibble()
  }
  return(time_X)
}

# --------------------------------------------------------------------------------------------------------------
# we find the cluster assignments for all isolates at a particular height at TP2 
# that can be found in the given list of clusters (have not yet originating clusters for these)
clusterAssignments <- function(h, time2_data, key_clusters, novels) {
  tmp <- time2_data %>% dplyr::select(isolate, all_of(h)) %>% 
    set_colnames(c("isolates", "clusters")) %>% 
    dplyr::filter(clusters %in% key_clusters)
  
  # TEST THIS PART FURTHER
  # we group by cluster and then count the number of isolates found in each cluster, 
  tmp2 <- tmp %>% select(-isolates) %>% 
    countCases(., "tp2_cl_size") %>% 
    left_join(tmp, ., by = "clusters") %>% 
    set_colnames(c("isolate", "tp2_cl", "tp2_cl_size"))
  
  # We filter so it only includes the original isolates (but still has the actual TP2 cluster sizes)
  # then return a dataset with || isolate | height (at tp2) | cluster (at tp2) | cluster size (at tp2) ||
  # We will use this dataset to match with TP1 potential originating clusters (just_originals)
  tmp2 %>% dplyr::filter(!(isolate %in% novels)) %>% 
    add_column(tp2_h = h, .after = 1) %>% return()
}

# --------------------------------------------------------------------------------------------------------------

checkID <- function(identifier, kc, x1, ids) {
  if (identifier %in% ids) {
    # if this originating cluster has been seen before, return an unflagged result
    kc %>% add_column(x1, .after = 1) %>% add_column(flagged = NA) %>% return()
  }else {
    # if this originating cluster has never been seen before, return a flagged result, with the ID
    kc %>% add_column(x1, .after = 1) %>% add_column(flagged = identifier) %>% return()
  }  
}
# --------------------------------------------------------------------------------------------------------------

countCases <- function(dataset, last_col) {
  df <- dataset %>% 
    group_by_all() %>% 
    summarise(n = n(), .groups = "drop")
  colnames(df)[ncol(df)] <- last_col
  df %>% return()
}

meltData <- function(dataset, id_val) {
  melt(dataset, id = id_val) %>% as_tibble() %>% return()
}

factorToInt <- function(dataset, cname) {
  dataset %>% 
    mutate(across(all_of(cname), as.character)) %>% 
    mutate(across(all_of(cname), as.integer)) %>% return()
}

createID <- function(df, c1, c2) {
  df %>% add_column(id = paste0(pull(df, c1), "-", pull(df, c2))) %>% return()
}
# --------------------------------------------------------------------------------------------------------------
resultsProcess <- function(h, time2_data, all_clusters, df1, ids, jo) {
  # given the heights and cluster assignments for all TP1 isolates, we filter to keep only the genomes found 
  # in a particular cluster kc, then we identify all the clusters at TP1 that only contain these genomes
  # and then return the first height-and cluster pair where these isolates are found in a single cluster  
    
  lapply(1:length(all_clusters), function(i) {
    print(paste0("h_", h, "-", i, "/", length(all_clusters)))
    kc <- jo %>% dplyr::filter(tp2_cl == all_clusters[i])
      
    # the first cluster in the TP1 dataset to contain only the originals found in the 
    # TP2 cluster, cluster_x, in form || tp1_h | tp1_cl | tp1_cl_size ||
    x1 <- df1 %>% 
      dplyr::filter(isolate %in% kc$isolate) %>% 
      dplyr::select(-isolate) %>%
      countCases(., "tp1_cl_size") %>% 
      dplyr::filter(tp1_cl_size == nrow(kc)) %>% 
      slice(tp1_cl_size = 1)
      
    # we add this "originating cluster" to the dataframe with this TP2 cluster
    # if this "originating cluster" is already in the ID list, return with an unflagged notation, 
    # otherwise add to the list of IDs and flag the cluster
    paste0(x1$tp1_h, "-", x1$tp1_cl) %>% 
      checkID(., kc, x1, ids) %>% 
      return()
  }) %>% bind_rows() %>% return()
}

# --------------------------------------------------------------------------------------------------------------

addToMetrics <- function(height, ids, metrics) {
  if (missing(metrics)) {
    # must be at height 0
    tibble(h = height, number_of_ids = length(ids)) %>% return()
  }else {
    metrics %>% add_row(tibble(h = height, number_of_ids = length(ids))) %>% return()
  }
}
