countingNovels <- function(time2, clusters) {
  tmp <- colnames(time2)[-1] %>% 
    lapply(., function(h_x) {
      h <- as.name(h_x)
      
      # at height h, the clusters that contain novel genomes
      h_cl <- pull(clusters, h) %>% unique()
      dfx <- time2 %>% dplyr::select(isolate, all_of(h)) %>%
        filter(!!h %in% h_cl)
      dfx$novel <- "original"
      dfx$novel[dfx$isolate %in% isolates_t2] <- "novel"
      
      # sizes of these clusters
      a <- dfx %>% group_by(!!h) %>% tally() %>% set_colnames(c("clusters","cl_size"))
      
      # number of novel genomes and original genomes in each cluster
      b <- dfx %>% group_by(!!h, novel) %>% tally()
      
      # clusters with novel and original genomes
      multistrain_cl <- b %>% filter(novel == "original") %>% pull(h)
      b %<>% filter(!!h %in% multistrain_cl) %>% set_colnames(c("clusters", "type", "count"))
      
      left_join(b, a, by = "clusters") %>%
        bind_cols(height = as.character(h), .) %>% return()
    }) %>% bind_rows()
  return(tmp)
}


absorbingClusters <- function(isolates_t2, time2, novel_counts) {
  notable_clusters <- isolates_t2 %>%
    lapply(., function(isolate_x) {
      isolate_x_df <- time2 %>% filter(isolate == isolate_x) %>%
        dplyr::select(-isolate) %>% t() %>% as.data.frame() %>%
        rownames_to_column() %>% as_tibble() %>%
        set_colnames(c("height","clusters"))
      
      absorbing_cluster <- left_join(isolate_x_df, novel_counts,
                                     by = c("height", "clusters")) %>%
        bind_cols(isolate = isolate_x)
      
      # the first cluster-height pair for isolate_x that includes non-novel genomes
      absorbing_cluster[complete.cases(absorbing_cluster),][1:2,] %>% return()
    }) %>% bind_rows()
  
  h_num <- notable_clusters$height %>% gsub("h_", "", .) %>% as.numeric()
  notable_clusters %<>% add_column(., h_num, .before = 1)
  
  return(notable_clusters)
}

identifyTP2Isolates <- function(key_clusters_T2, time2) {
  # for each of the absorbing clusters, we identify the isolates found in them at TP2
  T2_key_assign <- (1:nrow(key_clusters_T2)) %>%
    lapply(., function(i) {
      row_x <- key_clusters_T2[i,]
      time2 %>% filter(!!as.name(row_x$T2_h) == row_x$T2_cl) %>%
        dplyr::select(isolate, row_x$T2_h) %>% return()
    }) %>% set_names(key_clusters_T2$id)
  return(T2_key_assign)
}

clusterTP1 <- function(T2_key_assign, time1) {
  # identifies cluster at TP1 that absorbs the novel genomes, pairs with TP2 height-cluster id
  key_clusters_T1 <- (1:length(T2_key_assign)) %>%
    lapply(., function(i) {
      # timepoint 2 isolates in key_cluster i
      iso_in_t2 <- T2_key_assign[[i]] %>% pull(isolate)    # timepoint 2 isolates
      dfx <- time1 %>% filter(isolate %in% iso_in_t2) # T1 cluster assignments
      
      # first T1 height at which all T2 isolates of key cluster are in one T1 cluster
      T1_h <- dfx %>% melt(id = "isolate") %>% as_tibble() %>%
        group_by(variable, value) %>% tally() %>%
        filter(n == nrow(dfx)) %>% `[[`(1,1) %>% as.character()
      
      id <- T2_key_assign[i] %>% names() %>% as.character()
      T1_cl <- dfx %>% pull(T1_h) %>% unique()
      
      T1_cl_actual_size <- time1 %>% filter(!!as.name(T1_h) == T1_cl) %>% nrow()
      T1_cl_wo_novel_size <- dfx %>% pull(T1_h) %>% length()
      
      tibble(id, T1_h, T1_cl, T1_cl_wo_novel_size, T1_cl_actual_size) %>% return()
    }) %>% bind_rows()
  return(key_clusters_T1)
}

prepKeyClusters <- function(key_clusters_T1, key_clusters_T2) {
  key_clusters <- full_join(key_clusters_T2, key_clusters_T1, by = "id")
  
  key_clusters$m1 <- key_clusters$T2_cl_size - key_clusters$T1_cl_actual_size
  key_clusters$m2 <- key_clusters$T2_cl_size - key_clusters$T1_cl_wo_novel_size
  
  key_clusters$T1_cl %<>% as.character()
  key_clusters$T2_cl %<>% as.character()
  
  key_clusters$h1 <- key_clusters$T1_h %>% gsub("h_", "", .) %>% as.integer()
  key_clusters$h2 <- key_clusters$T2_h %>% gsub("h_", "", .) %>% as.integer()
  
  key_clusters %<>% 
    dplyr::select(c("id", "h1", "T1_h", "T1_cl", "T1_cl_wo_novel_size", 
                    "T1_cl_actual_size", "h2", "T2_h", "T2_cl", 
                    "T2_cl_size", "m1", "m2")) %>% rownames_to_column("x")
  key_clusters$x %<>% as.integer()
  return(key_clusters)
}

percForSize <- function(gth, key_clusters, size_x) {
  lapply(gth$perc, function(g) {
    key_clusters %>% filter(T1_cl_wo_novel_size >= size_x) %>%
      filter(prop_inc >= g) %>% nrow()
  }) %>% unlist() %>% set_names(gth$perc) %>%
    as.data.frame() %>% rownames_to_column() %>% as_tibble() %>%
    set_colnames(c("perc", paste0("num_cl_size_above_", size_x))) %>%
    return()
}

