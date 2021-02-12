
source("functions/formattingdata.R")

compsSet <- function(tp_coded, tp, indicate_progress) {
  # sort the data and assign general labels
  cb <- tp_coded %>% arrange(tp_h, tp_cl, isolate) %>% set_colnames(c("isolate", "tp_h", "tp_cl", "id"))
  allheights <- cb$tp_h %>% unique()
  
  # This part returns a data frame of CLUSTER ID | CLUSTER COMPOSITION | CLUSTER SIZE, 
  # for all heights and clusters at the given time point
  if (indicate_progress) {pb <- txtProgressBar(min = 0, max = length(allheights), initial = 0, style = 3)}
  tmp <- lapply(1:length(allheights), function(i) {
    if (indicate_progress) {setTxtProgressBar(pb, i)}
    
    # for height i, we first collect the data from tp_coded, then arrange by isolate
    x1 <- cb %>% filter(tp_h == allheights[i]) %>% arrange(isolate)
    x <- x1 %>%
      add_column(num_iso = x1$isolate %>% gsub("-", "", .) %>% as.integer()) %>%
      arrange(num_iso)
    
    # if at a height there is only one cluster, we handle it differently
    num_clusters <- x$id %>% unique() %>% length()
    # if there are multiple clusters at a given height:
    if (num_clusters > 1) {
      # output a table where each cluster has its isolate composition to the right: 
      #   isolate 1, isolate 2, isolate 3, ...
      a2 <- aggregate(isolate ~ id, data = x, FUN = toString) %>% as_tibble() %>% 
        set_colnames(c("id", "composition"))
      a3 <- aggregate(isolate ~ id, data = x, FUN = length) %>% as_tibble() %>% 
        set_colnames(c("id", "size"))
      left_join(a2, a3, by = "id") %>% return()
    }else {
      x$isolate %>% paste0(., collapse=",") %>% 
        tibble(id = unique(x$id), composition = ., size = length(x$isolate)) %>% return()
    }
  }) %>% bind_rows()
  
  if (indicate_progress) {close(pb)}
  
  # This part uses the cluster IDs to get the height and cluster as separate columns, 
  # then binds it with the cluster composition table
  tpcomps <- paste0(c(tp, "_h", "_"), collapse = "|") %>% gsub(., "", tmp$id) %>% 
    str_split_fixed(., "c", 2) %>% set_colnames(c("tp_h", "tp_cl")) %>% as_tibble() %>% 
    bind_cols(., tmp) %>% mutate(across(c(tp_h, tp_cl), as.integer))
  
  if (tp == "TP1") {
    tpcomps %>% rename(tp1_h = tp_h, tp1_cl = tp_cl, tp1_cl_size = size) %>% return()
  }else {
    tpcomps %>% rename(tp2_h = tp_h, tp2_cl = tp_cl, tp2_cl_size = size) %>% return()
  }
}

noChange <- function(ca, cb, hb) {
  in_both_heights <- ca %>% 
    inner_join(., cb, by = "comp") %>% 
    select(h_aft, cl_aft, size_aft, id_bef, id_aft)
  
  stayed_the_same <- left_join(in_both_heights, hb, by = c("id_bef" = "id")) %>% 
    select(-tp1_h, -tp1_cl, -tp1_cl_size, -id_bef) %>% 
    rename(tp1_h=h_aft, tp1_cl=cl_aft, tp1_cl_size=size_aft, id=id_aft) %>% 
    select(colnames(hb))
  # stayed_the_same$flagged_heights <- stayed_the_same$flagged_heights + 1  
  
  # not found in TP2 in some capacity
  inds <- which(is.na(stayed_the_same$tp2_h))
  if (length(inds) > 0) {stop("Did not stay the same")}
  # stayed_the_same$tp2_cl_size[inds] <- stayed_the_same$num_novs[inds] <- 0
  return(stayed_the_same)
}

checkEachIsolate <- function(cluster_i, t2_coded, t2_comps) {
  # no matching TP2 cluster found - the order of the cluster composition is muddling this up
  if (length(grep(" ", cluster_i$composition)) == 0) {
    isolates <- strsplit(cluster_i$composition, split = ",") %>% unlist()
  }else {
    isolates <- strsplit(cluster_i$composition, split = ", ") %>% unlist()  
  }
  ids <- t2_coded %>% filter(isolate %in% isolates) %>%
    pull(id) %>% table() %>% 
    as.data.frame() %>% as_tibble() %>% set_colnames(c("id", "size")) %>% 
    filter(size == length(isolates)) %>% pull(id)
  df2 <- t2_comps %>% filter(id %in% ids)
  
  cluster_i %>% rmIDComp() %>% bind_cols(., df2) %>% rmIDComp() %>% return()
}

# When using grep: if looking for 1, will return 1, 10, 214, etc. So we sandwich with hyphens: -1-,-10-,...
trackClusters <- function(hdata, t2_comps, t2names, t1_coded, t2_coded, indicate_progress) {
  
  singletons <- hdata %>% filter(tp1_cl_size == 1)
  t1set <- t2set <- tibble()
  if (nrow(singletons) > 0) { # at least one singleton cluster
    isos_in_tp1 <- t1_coded %>% filter(id %in% singletons$id) %>% rename(tp1_id = id)
    
    a3 <- t2_coded %>% filter(isolate %in% isos_in_tp1$isolate) %>% 
      rename(tp2_id = id) %>% left_join(isos_in_tp1, ., by = "isolate") %>% select(-isolate)
    
    a4 <- singletons %>% select(id, tp1_cl_size) %>% rename(tp1_id = id) %>% left_join(a3, ., by = "tp1_id")
    
    # basic TP1 data | basic TP2 data | number novels
    # i.e. the singletons are tracked to the TP2 clusters they're found in
    t1set <- t2_comps %>% select(id, tp2_cl_size, num_novs) %>% rename(tp2_id = id) %>% 
      left_join(a4, ., by = "tp2_id") %>% select(tp1_h, tp1_cl, tp1_cl_size, tp2_h, tp2_cl, tp2_cl_size, num_novs)
  }
  
  hx <- hdata %>% filter(tp1_cl_size > 1)
  if (nrow(hx) > 0) { # at least one cluster with size larger than 1
    if (indicate_progress) {
      tc <- txtProgressBar(min = 0, max = nrow(hx), initial = 0, style = 3)  
    }
    t2set <- lapply(1:nrow(hx), function(i) {
      if (indicate_progress) {setTxtProgressBar(tc, i)}
      cluster_i <- hx[i,]
      
      if (nchar(cluster_i$composition) > 2000) {
        results_i <- checkEachIsolate(cluster_i, t2_coded, t2_comps)
      }else {
        inds <- grep(cluster_i$composition, t2_comps$composition)
        results_i <- t2_comps[inds,] %>% rmIDComp() %>% bind_cols(cluster_i, .) %>% rmIDComp()
        
        if (nrow(results_i) == 0) {
          results_i <- checkEachIsolate(cluster_i, t2_coded, t2_comps)
        }
      }
      # if the process is being cut off partway
      # not tracking to all heights (e.g. TP2 height 1634 has not been found)
      if (!(last(t2names) %in% results_i$tp2_h)) {
        results_i <- checkEachIsolate(cluster_i, t2_coded, t2_comps)
      }
      results_i %>% arrange(tp2_h, tp2_cl) %>% return()
      
    }) %>% bind_rows()
    if (indicate_progress) {close(tc)}
    if (any(is.na(t2set))) {message("\nNA values!\n")}
  }
  
  bind_rows(t1set, t2set) %>% arrange(tp1_h, tp1_cl) %>% return()
}

# additional TP1s are the TP1 strains that were not found in the TP1 cluster, 
# but are found in the TP2 cluster (and are not novels)
additionalTP1 <- function(b1, b2, id1, id2, novs) {
  # id1 <- w$tp1_id[i]; id2 <- w$tp2_id[i]; novs <- novels
  tp1_isos <- b1 %>% filter(tp1_id == id1) %>% pull(isolate)
  tp2_isos <- b2 %>% filter(tp2_id == id2) %>% pull(isolate)
  m2 <- length(tp2_isos)
  m2a <- length(which(tp2_isos %in% novs))
  m2b <- length(which(tp2_isos %in% tp1_isos))
  return(tibble(tp1_id = id1, tp2_id = id2, add_tp1 = m2 - m2a - m2b))
}

oneHeight <- function(h_i, novels, t1_composition, t2_composition, oneh, b1, b2) {
  # check: (no cluster numbers skipped) - none should be skipped by definition, but just in case
  # identical(unique(oneh$tp1_cl), min(oneh$tp1_cl):max(oneh$tp1_cl))
  
  # first TP2 match for each TP1 cluster at this height
  d2 <- oneh[diff(c(0, oneh$tp1_cl)) != 0,] %>% rename(tp1_id = id) %>% 
    createID(., "tp2", "tp2_h", "tp2_cl") %>% rename(tp2_id = id)
  matched <- d2 %>% select(tp1_id, tp2_id)
  
  e2 <- t2_composition %>% select(id, composition) %>% rename(tp2_id = id, comp2 = composition)
  e1 <- t1_composition %>% select(id, composition) %>% rename(tp1_id = id, comp1 = composition)
  
  compmatches <- e1 %>% filter(tp1_id %in% matched$tp1_id) %>% 
    left_join(., matched, by = "tp1_id") %>% left_join(., e2, by = "tp2_id")
  
  did_not_change <- compmatches %>% filter(comp1 == comp2) %>% 
    select(tp1_id, tp2_id) %>% add_column(add_tp1 = 0)
  
  chg <- compmatches %>% filter(comp1 != comp2) %>% select(tp1_id, tp2_id)
  
  novs <- setdiff(b2$isolate, b1$isolate) %>% unique()
  changed_additional <- lapply(1:nrow(chg), function(i) {
    additionalTP1(b1, b2, chg$tp1_id[i], chg$tp2_id[i], novs)
  }) %>% bind_rows()
  
  c1 <- bind_rows(changed_additional, did_not_change) %>% 
    left_join(d2, ., by = c("tp1_id", "tp2_id"))
 
  c1 %>% add_column(
    actual_size_change = (c1$tp2_cl_size - c1$tp1_cl_size), 
    actual_growth_rate = ((c1$tp2_cl_size - c1$tp1_cl_size) / c1$tp1_cl_size) %>% round(., digits = 3),
    new_growth = (c1$tp2_cl_size / (c1$tp2_cl_size - c1$num_novs)) %>% round(., digits = 3)) %>% 
    rename(additional_tp1 = add_tp1) %>% 
    select(tp1_id, tp1_h, tp1_cl, tp1_cl_size, flag, last_flag, tp2_id, tp2_h, tp2_cl, tp2_cl_size, 
           additional_tp1, num_novs, actual_size_change, actual_growth_rate, new_growth) %>% 
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
    leadingZeros(., "tp1_cl", "c") %>% leadingZeros(., "tp2_cl", "c") %>%
    leadingZeros(., "tp1_h", "h", w = nchar(max(t1_composition$tp1_h))) %>%
    leadingZeros(., "tp2_h", "h", w = nchar(max(t2_composition$tp2_h))) %>% return()
}

