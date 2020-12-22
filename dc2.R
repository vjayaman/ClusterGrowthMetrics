x <- c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress", "stringr")
lapply(x, require, character.only = TRUE)

# given the defining filename, read in the data (need the full path from your working directory)
readData <- function(filename, file_number) {
  if (is.na(filename)) {
    stop(paste0("Time point ", file_number, " dataset not found."))
  }else {
    read.csv(file = filename, stringsAsFactors = FALSE, numerals = "no.loss", 
             check.names = FALSE) %>% as_tibble() %>% return()
  }
}

factorToInt <- function(dataset, cname) {
  dataset %>% 
    mutate(across(all_of(cname), as.character)) %>% 
    mutate(across(all_of(cname), as.integer)) %>% return()
}

countRows <- function(df, newname) {
  df2 <- df %>% 
    group_by(tp2_h, tp2_cl) %>% 
    summarise(n = n(), .groups = "drop")
  colnames(df2)[which(colnames(df2) == "n")] <- newname
  df2 %>% return()
}

createID <- function(df, c1, c2) {
  df %>% add_column(id = paste0(pull(df, c1), "-", pull(df, c2))) %>% return()
}

# tp_coded <- t1_coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>% arrange(tp_h, tp_cl, isolate)
compsSet <- function(tp_coded, tp) {
  cb <- tp_coded %>% set_colnames(c("isolate", "tp_h", "tp_cl", "id"))
  a1 <- cb$tp_h %>% unique()
  
  pb <- txtProgressBar(min = 0, max = length(a1), initial = 0, style = 3)
  tmp <- lapply(1:length(a1), function(i) {
    setTxtProgressBar(pb, i)
    x <- cb %>% filter(tp_h == a1[i]) %>% arrange(isolate)
    if (length(unique(x$id)) > 1) {
      a2 <- aggregate(isolate ~ id, data = x, FUN = toString) %>% as_tibble() %>% 
        set_colnames(c("id", "composition"))
      a3 <- aggregate(isolate ~ id, data = x, FUN = length) %>% as_tibble() %>% 
        set_colnames(c("id", "size"))
      left_join(a2, a3, by = "id") %>% return()
    }else {
      x$isolate %>% paste0(., collapse=",") %>% 
        tibble(id = unique(x$id), composition = ., 
               size = length(x$isolate)) %>% return()
    }
  }) %>% bind_rows()
  close(pb)
  
  tpcomps <- str_split_fixed(tmp$id, "-", 2) %>% 
    set_colnames(c("tp_h", "tp_cl")) %>% 
    as_tibble() %>% 
    bind_cols(., tmp) %>% 
    mutate(across(c(tp_h, tp_cl), as.integer))
  
  if (tp == 1) {
    tpcomps %>% rename(tp1_h = tp_h, tp1_cl = tp_cl, tp1_cl_size = size) %>% return()
  }else {
    tpcomps %>% rename(tp2_h = tp_h, tp2_cl = tp_cl, tp2_cl_size = size) %>% return()
  }
}
# ca <- aft_comps
# cb <- bef_comps
noChange <- function(ca, cb, hb) {
  # cb <- t1_comps %>% filter(tp1_h == 0)
  # ca <- t1_comps %>% filter(tp1_h == 1)
  
  in_both_heights <- ca %>% 
    inner_join(., cb, by = "comp") %>% 
    select(h_aft, cl_aft, size_aft, id_bef, id_aft)
  
  # f1 <- left_join(in_both_heights, hb, by = c("id_bef" = "id")) %>% 
  #   select(-tp1_h, -tp1_cl, -tp1_cl_size, -id_bef)
  
  stayed_the_same <- left_join(in_both_heights, hb, by = c("id_bef" = "id")) %>% 
    select(-tp1_h, -tp1_cl, -tp1_cl_size, -id_bef) %>% 
    rename(tp1_h=h_aft, tp1_cl=cl_aft, tp1_cl_size=size_aft, id=id_aft) %>% 
    select(colnames(hb))
  stayed_the_same$flagged_heights <- stayed_the_same$flagged_heights + 1  
  
  # not found in TP2 in some capacity
  inds <- which(is.na(stayed_the_same$tp2_h))
  stayed_the_same$tp2_cl_size[inds] <- stayed_the_same$num_novs[inds] <- 0
  return(stayed_the_same)
}

checkEachIsolate <- function(allc, t2_comps) {
  if (length(allc) > 1) {
    ids1 <- t2_comps$id[grep(allc[1], t2_comps$composition)]
    
    for (i in 2:length(allc)) {
      if (length(ids1) == 0) {message("\nRan out of TP2 clusters - check composition!\n")}
      newset <- t2_comps %>% filter(id %in% ids1)
      checkclusters <- grep(allc[i], newset$composition)
      ids <- t2_comps %>% 
        filter(id %in% newset$id[checkclusters]) %>% 
        pull(id)
    }
    return(ids)  
  }else {
    message("Only one element in the cluster - different issue")
  }
}

rmIDComp <- function(df) {
  df %>% select(-id, -composition) %>% return()
}
#   x <- df$composition %>% strsplit(split = ",") %>% unlist() %>% gsub(" ", "", .)
#   lapply(x, function(xi) grep(xi, t2_comps$composition)) %>% 
#     Reduce(intersect, .) %>% return()
# }
# hdata <- hdata %>% arrange(tp1_h, tp1_cl)
trackClusters <- function(hdata, t2_comps) {
  # FOUND A PROBLEM WITH USING GREP!
  #   - if you're looking for composition 1, it considers 10, 107, 1247, ... as hits
  #   - need to include the commas in our actual checks
  
  singletons <- hdata %>% filter(tp1_cl_size == 1)
  if (nrow(singletons) > 0) { # at least one singleton cluster
    # TP1 clusters that contain these TP1 singletons for this height
    a1 <- t1_coded %>% filter(id %in% singletons$id) %>% rename(tp1_id = id)
    # the TP2 clusters that contain at least the TP1 isolates from the above TP1 clusters
    a2 <- t2_coded %>% filter(isolate %in% a1$isolate) %>% rename(tp2_id = id)
    # the TP1 singletons matched to TP2 clusters that contain at least those genomes
    a3 <- left_join(a1, a2, by = "isolate") %>% select(tp1_id, tp2_id)
      
    a4 <- t2_comps %>% filter(id %in% a3$tp2_id) %>% 
      left_join(., a3, by = c("id" = "tp2_id")) %>% 
      select(-composition, -id)
    t1set <- singletons %>% rename(tp1_id = id) %>% 
      select(-composition) %>% 
      right_join(., a4, by = "tp1_id") %>% select(-tp1_id)
  }
  
  hx <- hdata %>% filter(tp1_cl_size > 1)
  
  if (nrow(hx) > 0) { # at least one cluster with size larger than 1
    pb <- txtProgressBar(min = 0, max = nrow(hx)+1, initial = 0, style = 3)
    k <- 0
    t2set <- lapply(1:nrow(hx), function(i) {
      k <- k + i
      setTxtProgressBar(pb, k)
      i <- 1
      cluster_i <- hx[i,]
      results_i <- t2_comps[grep(cluster_i$composition, t2_comps$composition),] %>% 
        rmIDComp() %>% bind_cols(cluster_i, .) %>% rmIDComp()
      
      if (nrow(results_i) == 0) {
        # no matching TP2 cluster found - the order of the cluster composition is muddling this up
        isolates <- strsplit(cluster_i$composition, split = ", ") %>% unlist()
        ids <- t2_coded %>%
          filter(isolate %in% isolates) %>%
          pull(id) %>% table() %>%
          as.data.frame() %>% as_tibble() %>%
          set_colnames(c("id", "size")) %>%
          filter(size == length(isolates)) %>%
          pull(id)
        df2 <- t2_comps %>% filter(id %in% ids)
        
        results_i <- cluster_i %>% rmIDComp() %>% bind_cols(., df2) %>% rmIDComp()
      }
      return(results_i)
    }) %>% bind_rows()
    setTxtProgressBar(pb, k + 1)
    close(pb)
    
    if (any(is.na(t2set))) {
      message("\nNA values!\n")  
    }
  }
  
  if (nrow(t1set) > 0 & nrow(t2set) > 0) {
    bind_rows(t1set, t2set) %>% return()
  }else if (nrow(t1set) == 0 & nrow(t2set) > 0) {
    # message("Note: no singletons")
    t2set %>% return()
  }else if (nrow(t1set) > 0 & nrow(t2set) == 0) {
    # message("Note: only singletons")
    t1set %>% return()
  }else {
    # nrow(t1set) == 0 & nrow(t2set) == 0
    message("No results!")
  }
  
  # pb <- txtProgressBar(min = 0, max = nrow(hx)+1, initial = 0, style = 3)
  # k <- 0
  # tmp <- lapply(1:nrow(hx), function(i) {
  #   k <- k + i
  #   setTxtProgressBar(pb, k)
  #   # the composition (and cluster size) data for the cluster in row i of 
  #   # the TP1 composition data (for the height in hx)
  #   cxdata <- hx[i,]
  #   # a1 <- t2_comps[grep(cxdata$composition, t2_comps$composition),] %>% select(-composition, -id)
  #   # cases:
  #   # number of characters in the cluster's pattern of membership is > 2000
  #   # number of characters in the cluster's pattern of membership is < 2000
  #   #   the order of clusters composition is off
  #   #   e.g. in T1_5_2:     1795, 1796, 1797, 1806, 1807
  #   #        in T2_232_635: 1795, 1796, 1797, 1798, 1801, 1802, 1803, 1804, ...
  #   #   in these cases: need to check for each member of the clusters
  #   
  #   allc <- strsplit(cxdata$composition, split = ",") %>% unlist() %>% str_trim(., "both")
  #   ids <- NULL
  #   
  #   if (nchar(cxdata$composition) < 2000) {
  #     # indices of the TP2 composition dataset where clusters contain >= all the isolates from the cxdata cluster
  #     inds <- grep(cxdata$composition, t2_comps$composition)
  #     ids <- t2_comps$id[inds]
  #   }else {
  #     # checking each cluster member, using recursion
  #     ids <- checkEachIsolate(allc, t2_comps)
  #   }
  #   
  #   # if the composition of a cluster's members at TP2 is in a different order than at TP1
  #   if (length(ids) == 0) {
  #     ids <- checkEachIsolate(allc, t2_comps)
  #   }
  #   
  #   if (length(ids) == 0) {message("\nLength of indices is still zero! Some other case not considered.\n")}
  #   
  #   kc <- t2_comps %>% filter(id %in% ids) %>% select(-composition, -id)
  #   
  #   # merged the data for the single TP1 clusters and all corresponding TP2 clusters
  #   cxdata %>% select(-composition, -id) %>% bind_cols(., kc) %>% return()
  # }) %>% bind_rows()
  # setTxtProgressBar(pb, k + 1)
  # close(pb)
  # return(tmp)
}

# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
## FOR USER: replace the filename variables with quoted file paths if you don't want to input them each time
# # the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||33333
# # the raw datasets, no filtering or other changes made

t1_coded <- "data/timepoint1_data.csv" %>% readData(., 1) %>% 
  mutate(isolate = 1:nrow(.)) %>% 
  melt(id = "isolate") %>% as_tibble() %>% 
  set_colnames(c("isolate", "tp1_h", "tp1_cl")) %>% 
  factorToInt("tp1_h") %>% 
  createID(., "tp1_h", "tp1_cl")
t1_coded$isolate <- paste0("-", t1_coded$isolate, "-")

t2_coded <- "data/timepoint2_data.csv" %>% readData(., 2) %>% 
  mutate(isolate = 1:nrow(.)) %>% 
  melt(id = "isolate") %>% as_tibble() %>% 
  set_colnames(c("isolate", "tp2_h", "tp2_cl")) %>% 
  factorToInt("tp2_h") %>% 
  createID(., "tp2_h", "tp2_cl")
t2_coded$isolate <- paste0("-", t2_coded$isolate, "-")

# t1_comps <- t1_coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>%
#   arrange(tp_h, tp_cl, isolate) %>%
#   compsSet(., 1)
# # saveRDS(t1_comps, "t1_comps.Rds")
t1_comps <- readRDS("t1_comps.Rds")

# t2_comps <- t2_coded %>% rename(tp_h = tp2_h, tp_cl = tp2_cl) %>%
#   arrange(tp_h, tp_cl, isolate) %>%
#   compsSet(., 2)

novels <- setdiff(t2_coded$isolate, t1_coded$isolate)
counting_novels <- t2_coded %>%
  filter(isolate %in% novels) %>%
  group_by(id) %>%
  summarise(num_novs = n(), .groups = "drop")

# t2_comps <- t2_comps %>% left_join(., counting_novels, by = "id")
# t2_comps$num_novs[is.na(t2_comps$num_novs)] <- 0
# # saveRDS(t2_comps, "t2_comps.Rds")
t2_comps <- readRDS("t2_comps.Rds")
# --------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------
# BASE CASE
h_before <- 0
message(paste0("Collecting data for height ", h_before))
hdata <- t1_comps %>% filter(tp1_h == h_before) %>% arrange(tp1_h, tp1_cl)
cc <- trackClusters(hdata, t2_comps)

# now need to make metric calculations and extract the best result:
#     - the TP2 cluster that has the largest proportion of novels to growth rate
#     i.e. the "best" TP2 cluster where we look at clusters that absorbed lots of novels very quickly
# the actual growth of these clusters (TP2 size - TP1 size) / TP1 size
cc$actual_growth_rate <- (cc$tp2_cl_size - cc$tp1_cl_size) / cc$tp1_cl_size
# the number of novels in the TP2 cluster over the growth rate --> growth acceleration?
cc$acc <- cc$num_novs / cc$actual_growth_rate
cc$acc[is.na(cc$acc)] <- 0
cc$flag <- NA

# note: the way these are flagged, we can see how many TP2 clusters were formed from the 
# selected TP1 cluster before we get to the "key growth" where a cluster went from its state
# at TP1 to its state at TP2 with a remarkable amount of growth
clx <- unique(cc$tp1_cl)
k <- 0
pb <- txtProgressBar(min = 0, max = length(clx), initial = 0, style = 3)
hb <- lapply(1:length(clx), function(i) {
  k <- k + i
  setTxtProgressBar(pb, k)
  df <- cc %>% filter(tp1_cl == clx[i]) %>% arrange(tp2_h, tp2_cl)
  df$flag <- 1:nrow(df)
  df$acc[which(df$actual_growth_rate == 0)] <- 0
  df %>% arrange(-acc, tp2_h, tp2_cl) %>% slice(1) %>% return()
}) %>% bind_rows() %>% 
  createID(., "tp1_h", "tp1_cl") %>% 
  add_column(flagged_heights = 0)
close(pb)  

bef_comps <- t1_comps %>% filter(tp1_h == h_before) %>%
  set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))

alldata <- hb
# --------------------------------------------------------------------------------------------------------
# h4 is too slow - too many clusters? clusters too small? why is is so much slower than for h2 or h3?
# --------------------------------------------------------------------------------------------------------
# ALL OTHER HEIGHTS
heights <- unique(t1_comps$tp1_h)
for (h_after in heights[-1]) {
  message(paste0("Collecting data for height ", h_after))

  aft_comps <- t1_comps %>% filter(tp1_h == h_after) %>% 
    set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft"))

  print("part1")
  ss <- noChange(aft_comps, bef_comps, hb)
  
  print("part2")
  cc <- aft_comps %>% 
    filter(!(id_aft %in% ss$id)) %>% 
    set_colnames(colnames(t1_comps)) %>% 
    trackClusters(., t2_comps) %>% 
    add_column(actual_growth_rate = NA, acc = NA, flag = NA)
  
  print("part3")
  if (nrow(cc) > 0) {
    natp2h <- which(is.na(cc$tp2_h))
    cc$actual_growth_rate <- (cc$tp2_cl_size - cc$tp1_cl_size) / cc$tp1_cl_size
    cc$actual_growth_rate[natp2h] <- cc$acc[natp2h] <- 0
    
    cc$acc <- cc$num_novs / cc$actual_growth_rate
    cc$acc[which(cc$actual_growth_rate == 0)] <- 0
    
    clx <- unique(cc$tp1_cl)
    hnew <- lapply(1:length(clx), function(i) {
      df <- cc %>% filter(tp1_cl == clx[i]) %>% arrange(tp2_h, tp2_cl)
      df$flag <- 1:nrow(df)
      df$acc[which(df$actual_growth_rate == 0)] <- 0
      df %>% arrange(-acc, tp2_h, tp2_cl) %>% slice(1) %>% return()
    }) %>% bind_rows() %>% 
      createID(., "tp1_h", "tp1_cl") %>% 
      add_column(flagged_heights = 0) %>% 
      bind_rows(., ss) %>% 
      arrange(tp1_h, tp1_cl)
  }else {
    hnew <- ss %>% arrange(tp1_h, tp1_cl)
  }

  print("part4")
  alldata <- bind_rows(alldata, hnew)
  hb <- hnew
  h_before <- h_after
  bef_comps <- aft_comps %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
  
  saveRDS(alldata, "alldata.Rds")
  saveRDS(hb, "hb.Rds")
  saveRDS(h_before, "h_before.Rds")
  saveRDS(bef_comps, "bef_comps.Rds")
}

# NA TP2 height:
# tp1_h = 1
# tp1_cl = 2
# tp1_cl_size = 14
# tp2_h = tp2_cl = NA

# restarted process at around 3:45 (at height 59 by 4:02, at height 125 by 4:06)
# got to height 197 before it ran into a problem
# alldata <- readRDS("alldata.Rds") %>% filter(tp1_h < 197)
# saveRDS(hb, "hb.Rds")
#   saveRDS(h_before, "h_before.Rds")
#   saveRDS(bef_comps, "bef_comps.Rds")


# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
# # some of the heights still have clusters with length zero
# # e.g. h309, 310
# 
# # Figuring out why some heights are not tracked properly
# 
# # alldata <- readRDS("alldata.Rds")
# # alldata %>% arrange(tp1_h, tp1_cl) %>% filter(is.na(tp2_h))
# # first case is at height 5, cluster 18
# 
# # # going to go back to that height and go through each cluster
# # alldata <- readRDS("alldata.Rds") %>% filter(tp1_h < 5) %>% arrange(tp1_h, tp1_cl)
# # hb <- readRDS("alldata.Rds") %>% filter(tp1_h == 4) %>% arrange(tp1_h, tp1_cl)
# # h_before <- 4
# # bef_comps <-  t1_comps %>% filter(tp1_h == 4) %>%
# #   set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
# 
# # going to go back to that height and go through each cluster
# alldata <- readRDS("alldata.Rds") %>% filter(tp1_h < 309) %>% arrange(tp1_h, tp1_cl)
# hb <- readRDS("alldata.Rds") %>% filter(tp1_h == 308) %>% arrange(tp1_h, tp1_cl)
# h_before <- 308
# bef_comps <-  t1_comps %>% filter(tp1_h == 308) %>%
#   set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
# 
# # -------------------------------------------------------------------------------------------------
# 

h_after <- 1
message(paste0("Collecting data for height ", h_after))

aft_comps <- t1_comps %>% filter(tp1_h == h_after) %>% 
  set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft"))

print("part1")
ss <- noChange(aft_comps, bef_comps, hb)
  
print("part2")
cc <- aft_comps %>% 
  filter(!(id_aft %in% ss$id)) %>% 
  set_colnames(colnames(t1_comps)) %>% 
  trackClusters(., t2_comps) %>% 
  add_column(actual_growth_rate = NA, acc = NA, flag = NA)
  
print("part3")
if (nrow(cc) > 0) {
  natp2h <- which(is.na(cc$tp2_h))
  cc$actual_growth_rate <- (cc$tp2_cl_size - cc$tp1_cl_size) / cc$tp1_cl_size
  cc$actual_growth_rate[natp2h] <- cc$acc[natp2h] <- 0
    
  cc$acc <- cc$num_novs / cc$actual_growth_rate
  cc$acc[which(cc$actual_growth_rate == 0)] <- 0
    
  clx <- unique(cc$tp1_cl)
  hnew <- lapply(1:length(clx), function(i) {
    df <- cc %>% filter(tp1_cl == clx[i]) %>% arrange(tp2_h, tp2_cl)
    df$flag <- 1:nrow(df)
    df$acc[which(df$actual_growth_rate == 0)] <- 0
    df %>% arrange(-acc, tp2_h, tp2_cl) %>% slice(1) %>% return()
  }) %>% bind_rows() %>% 
    createID(., "tp1_h", "tp1_cl") %>% 
    add_column(flagged_heights = 0) %>% 
    bind_rows(., ss) %>% 
    arrange(tp1_h, tp1_cl)
}else {
  hnew <- ss %>% arrange(tp1_h, tp1_cl)
}

print("part4")
alldata <- bind_rows(alldata, hnew)
hb <- hnew
h_before <- h_after
bef_comps <- aft_comps %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
  