
source("functions/formattingdata.R")

# --------------------------------------------------------------------------------------------------------------
# Given a TP dataset, we output a set where the composition of each cluster is indicated with a string where 
# each character represents a single isolate (easier to compare cluster composition later on)
# --------------------------------------------------------------------------------------------------------------
compsSet <- function(tp_coded, tp, indicate_progress) {
  # sort the data and assign general labels
  cb <- tp_coded %>% mutate(num_iso = as.integer(gsub("-", "", isolate)), 
                            composition = isolate, tp_cl_size = isolate) %>% arrange(tp_h, tp_cl, num_iso)
  allheights <- cb$tp_h %>% unique()
  
  # This part returns a data frame of CL ID | CL COMP | CL SIZE, for all heights and clusters at the given TP
  if (indicate_progress) {pb <- txtProgressBar(min = 0, max = length(allheights), initial = 0, style = 3)}
  
  tmp <- lapply(1:length(allheights), function(i) {
    if (indicate_progress) {setTxtProgressBar(pb, i)}
    
    # for height i, we first collect the data from tp_coded, then arrange by isolate
    x <- cb %>% filter(tp_h == allheights[i])
    
    # if there are multiple clusters at a given height:
    if (length(unique(x$tp_id)) > 1) {
      # output a table where each cluster has its isolate composition to the right: iso 1, iso 2, iso 3, ...
      a2 <- aggregate(composition ~ tp_id, data = x, FUN = toString)
      a3 <- aggregate(tp_cl_size ~ tp_id, data = x, FUN = length)
      left_join(a2, a3, by = "tp_id") %>% as_tibble() %>% return()
    }else {
      # if at a height there is only one cluster, we handle it differently
      tibble(tp_id = unique(x$tp_id), composition = paste0(x$isolate, collapse=","), 
             tp_cl_size = length(x$isolate)) %>% return()
    }
  }) %>% bind_rows()
  
  if (indicate_progress) {close(pb)}
  
  # This part binds the height and cluster columns with the cluster composition table
  tpcomps <- tp_coded %>% select(-isolate) %>% unique() %>% 
    left_join(tmp, ., by = "tp_id") %>% select(tp_h, tp_cl, tp_id, composition, tp_cl_size) %>% 
    set_colnames(gsub("tp", tolower(tp), colnames(.))) %>% return()
}

# --------------------------------------------------------------------------------------------------------------
# Identifying TP1 clusters at height b that did not change from height a --> they inherit the metrics, 
# since the values did not change
# --------------------------------------------------------------------------------------------------------------
noChange <- function(ca, cb, hb) {
  # these are clusters whose composition did not change from cb to ca (height x to height x+1)
  in_both_heights <- inner_join(ca, cb, by = "comp") %>% select(grep("aft", colnames(.)), id_bef)
  
  stayed_the_same <- left_join(in_both_heights, hb, by = c("id_bef"="tp1_id")) %>% 
    select(-grep("tp1", colnames(.)), -id_bef) %>% 
    rename(tp1_h = h_aft, tp1_cl = cl_aft, tp1_cl_size = size_aft, tp1_id = id_aft)
  
  assert("This dataset only has clusters that did not change", !any(is.na(stayed_the_same$tp2_h)))
  return(stayed_the_same)
}

# --------------------------------------------------------------------------------------------------------------
# When using grep doesn't return the needed results (several possible reasons), we track the individual 
# strains in a cluster (to find out which TP2 clusters contain at least these strains)
# --------------------------------------------------------------------------------------------------------------
checkEachIsolate <- function(cluster_i, t2_coded, t2_comps) {
  # no matching TP2 cluster found - the order of the cluster composition is muddling this up
  isolates <- strsplit(cluster_i$composition, split = ",|, ") %>% unlist()

  ids <- t2_coded %>% filter(isolate %in% isolates) %>% 
    group_by(tp2_id) %>% summarise(size = n(), .groups = "drop") %>% 
    filter(size == length(isolates)) %>% pull(tp2_id)
  
  df2 <- t2_comps %>% filter(tp2_id %in% ids) %>% select(-composition)
  cluster_i %>% select(-composition) %>% bind_cols(., df2) %>% return()
}

# --------------------------------------------------------------------------------------------------------------
# For a set of TP1 singletons, we track them to the TP2 clusters they're found in
# --------------------------------------------------------------------------------------------------------------
trackSingletons <- function(singles, t1_coded, t2_coded, composition_set) {
  isos_in_tp1 <- t1_coded %>% filter(tp1_id %in% singles$tp1_id)
  
  a4 <- t2_coded %>% filter(isolate %in% isos_in_tp1$isolate) %>% 
    left_join(isos_in_tp1, ., by = "isolate") %>% select(-isolate) %>% 
    left_join(., singles, by = "tp1_id")
  
  left_join(a4, composition_set, by = "tp2_id") %>% 
    select(tp1_h, tp1_cl, tp1_id, tp1_cl_size, tp2_h, tp2_cl, tp2_id, tp2_cl_size, num_novs) %>% return()
}

# --------------------------------------------------------------------------------------------------------------
# Tracking TP1 clusters to their TP2 equivalents, using grep primarily, with other / backup methods for 
# TP1 singletons and multistrain cases where the first method fails or is not comprehensive enough
# --------------------------------------------------------------------------------------------------------------
# When using grep: if looking for 1, will return 1, 10, 214, etc. So we sandwich with hyphens: -1-,-10-,...
trackClusters <- function(hdata, t2_comps, t2names, t1_coded, t2_coded, indicate_progress) {
  
  t1set <- t2set <- tibble()
  singletons <- hdata %>% filter(tp1_cl_size == 1) %>% select(tp1_id, tp1_cl_size)
  # at least one singleton cluster
  if (nrow(singletons) > 0) {
    t1set <- t2_comps %>% select(tp2_id, tp2_cl_size, num_novs) %>% 
      trackSingletons(singletons, t1_coded, t2_coded, .)
  }
  
  multistrain <- hdata %>% filter(tp1_cl_size > 1)
  if (nrow(multistrain) > 0) { # at least one cluster with size larger than 1
    if (indicate_progress) {tc <- txtProgressBar(min = 0, max = nrow(multistrain), initial = 0, style = 3)}
    
    t2set <- lapply(1:nrow(multistrain), function(i) {
      if (indicate_progress) {setTxtProgressBar(tc, i)}
      results_i <- tibble()
      cluster_i <- multistrain[i,]
      
      if (nchar(cluster_i$composition) < 2000) {
        inds <- grep(cluster_i$composition, t2_comps$composition)
        results_i <- t2_comps[inds,] %>% select(-composition) %>% 
          bind_cols(cluster_i, .) %>% select(-composition)
        
        if (nchar(cluster_i$composition) >= 2000 | (nrow(results_i) == 0) | (!last(t2names) %in% results_i$tp1_h)) {
          # (1) too many isolates in the cluster for effective grep use
          # (2) no matching TP2 cluster found - the order of the cluster composition is muddling this up
          # (3) if the process is being cut off partway - not tracking to all heights (e.g. TP2 height 1634 not found)
          results_i <- checkEachIsolate(cluster_i, t2_coded, t2_comps)
        }
      }
      results_i %>% arrange(tp2_h, tp2_cl) %>% return()
    }) %>% bind_rows()
    if (indicate_progress) {close(tc)}
    assert("\nNo NA values\n", all(!is.na(t2set)))
  }
  bind_rows(t1set, t2set) %>% arrange(tp1_h, tp1_cl) %>% return()
}

# --------------------------------------------------------------------------------------------------------------
# We introduce the metric columns, once we have collected all the tracking information for the TP1 clusters
# Include: cluster size change, growth rate (using actual sizes), growth rate (using number of novels), 
# number of sneakers, number of novels
# --------------------------------------------------------------------------------------------------------------
oneHeight <- function(novels, q1, q2, matched) {
  compmatches <- matched[!duplicated(matched$tp1_id),]
  
  did_not_chg <- compmatches %>% filter(tp1_cl_size == tp2_cl_size) %>% add_column(add_TP1 = 0)
  chg <- compmatches %>% filter(tp1_cl_size != tp2_cl_size)
  
  # identifying the number of additional TP1 strains (sneakers) that show up in the TP2 cluster 
  # that each TP1 cluster was first tracked to
  sneakers <- lapply(1:nrow(chg), function(j) {
    kc <- chg %>% slice(j) # key cluster j
    tbl1 <- q1 %>% filter(tp1_id == kc$tp1_id) %>% mutate(status = "tp1_cl_size")
    tbl2 <- q2 %>% filter(tp2_id == kc$tp2_id)
    # number of novels in the TP2 cluster it kc was tracked to
    x2 <- tbl2 %>% filter(status == "novs") %>% nrow()
    # number of sneakers in the TP2 cluster it kc was tracked to
    x3 <- tbl2 %>% filter(!(isolate %in% tbl1$isolate) & is.na(status)) %>% 
      mutate(status = "additional_TP1") %>% nrow()
    
    c2_tally <- tibble(tp1_cl_size = nrow(tbl1), num_novs = x2, add_TP1 = x3, 
                       tp2_cl_size = sum(nrow(tbl1), x2, x3))
    
    assert(paste0("oneHeight(): Novels check failed for ", kc$tp1_id), kc$num_novs == c2_tally$num_novs)
    assert(paste0("Composition not calculated properly for ", kc$tp2_id, 
                  " when tracking ", kc$tp1_id), c2_tally$tp2_cl_size == kc$tp2_cl_size)
    
    left_join(kc, c2_tally, by = c("tp1_cl_size", "tp2_cl_size", "num_novs")) %>% return()
  }) %>% bind_rows()
  
  bind_rows(did_not_chg, sneakers) %>% 
    mutate(actual_size_change = tp2_cl_size - tp1_cl_size, 
           actual_growth_rate = ((tp2_cl_size - tp1_cl_size) / tp1_cl_size) %>% round(., digits = 3), 
           new_growth = (tp2_cl_size / (tp2_cl_size - num_novs)) %>% round(., digits = 3)) %>% 
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% return()
}

