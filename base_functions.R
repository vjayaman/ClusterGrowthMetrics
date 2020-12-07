x <- c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress", "tcltk")
lapply(x, require, character.only = TRUE)

# DATA HANDLING -----------------------------------------------------------------------------------------------------

get <- .Primitive("[[")

getUserInput <- function(msg) {
  cat(msg)
  op = readLines(con = "stdin", 1)
}

# given the defining filename, read in the data (need the full path from your working directory)
readData <- function(filename, file_number) {
  if (is.na(filename)) {
    stop(paste0("Time point ", file_number, " dataset not found."))
  }else {
    read.csv(file = filename, stringsAsFactors = FALSE, numerals = "no.loss", 
             check.names = FALSE) %>% as_tibble() %>% return()
  }
}

# Used after datacollection.R is done, merges the separate height files into one and saves them as a 
# results file of the tracked clusters, to be used in preparingmetrics.R to generate the metrics file
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

meltData <- function(dataset, id_val) {
  melt(dataset, id = id_val) %>% as_tibble() %>% return()
}

factorToInt <- function(dataset, cname) {
  dataset %>% 
    mutate(across(all_of(cname), as.character)) %>% 
    mutate(across(all_of(cname), as.integer)) %>% return()
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

countCases <- function(dataset, last_col) {
  df <- dataset %>% 
    group_by_all() %>% 
    summarise(n = n(), .groups = "drop")
  colnames(df)[ncol(df)] <- last_col
  df %>% return()
}

# these are all the clusters and their composition
clustComp <- function(df, height, dtype) {
  message(paste0("Preparing cluster composition IDs for height ", height))
  a1 <- df %>% select(isolate, all_of(height)) %>% 
    set_colnames(c("isolate", "tp2_cl"))
  clusters <- a1$tp2_cl %>% unique()
  
  hx <- list("a1" = a1, "cl" = clusters)
  
  lapply(1:length(hx$cl), function(i) {
    filter(hx$a1, tp2_cl == hx$cl[i]) %>% pull(isolate) %>% 
      paste0(collapse = ",") %>% return()
  }) %>% unlist() %>% 
    tibble(., hx$cl) %>% set_colnames(c("composition", dtype)) %>% return()
}

# MESSAGES FOR USER -------------------------------------------------------------------------------------------------

collectionMsg <- function(h, time2_data, all_clusters) {
  if (match(h, colnames(time2_data)) == 2) { # if h == "0"
    message(paste0("Collecting data for the ", length(all_clusters), " clusters at height ", h, " (base case)"))
  }else {
    message(paste0("Collecting data for the ", length(all_clusters), " cluster(s) that ", 
                   "changed in composition from the previous height to the current height, ", h))    
  }
}

addToMetrics <- function(height, ids, metrics) {
  if (missing(metrics)) { # at height 0
    tibble(h = height, number_of_ids = length(ids)) %>% return()
  }else {
    metrics %>% add_row(tibble(h = height, number_of_ids = length(ids))) %>% return()
  }
}

timeTaken <- function(pt, sw) {
  t1 <- (sw[['end_time']] - sw[['start_time']])/60/60
  t2 <- abs(t1 - trunc(t1))*60
  t3 <- abs(t2 - trunc(t2))*60
  message(paste0("The ", pt, " process took ", trunc(t1), " hours, ", 
                 trunc(t2), " minutes, and ", trunc(t3), " seconds.\n"))
}

externalProgressBar <- function(x, i, msg) {
  setTkProgressBar(x, i, title = paste0("Progress of metrics prep: ", round(i/total*100, 0), "% done"), label = msg)
}

