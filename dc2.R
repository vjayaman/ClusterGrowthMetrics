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

## FOR USER: replace the filename variables with quoted file paths if you don't want to input them each time
# # the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||33333
# # the raw datasets, no filtering or other changes made
time1_raw <- "data/timepoint1_data.csv" %>% readData(., 1) %>% 
  melt(id = "isolate") %>% as_tibble() %>% 
  set_colnames(c("isolate", "tp1_h", "tp1_cl")) %>% 
  factorToInt("tp1_h")

time2_raw <- "data/timepoint2_data.csv" %>% readData(., 2) %>% 
  melt(id = "isolate") %>% as_tibble() %>% 
  set_colnames(c("isolate", "tp2_h", "tp2_cl")) %>% 
  factorToInt("tp2_h")

novels <- setdiff(time2_raw$isolate, time1_raw$isolate)
withnovs <- time2_raw %>% 
  filter(isolate %in% novels) %>% 
  group_by(tp2_h, tp2_cl) %>% 
  summarise(num_novs = n(), .groups = "drop")

t1clusters <- time1_raw %>% select(-isolate) %>% unique()
t1_raw <- time1_raw
t1_coded <- "data/timepoint1_data.csv" %>% readData(., 1) %>% 
  mutate(isolate = 1:nrow(.)) %>% 
  melt(id = "isolate") %>% as_tibble() %>% 
  set_colnames(c("isolate", "tp1_h", "tp1_cl")) %>% 
  factorToInt("tp1_h") %>% 
  createID(., "tp1_h", "tp1_cl")

t2_raw <- time2_raw
t2_raw$id <- paste0(t2_raw$tp2_h, "-", t2_raw$tp2_cl)

t2_coded <- "data/timepoint2_data.csv" %>% readData(., 2) %>% 
  mutate(isolate = 1:nrow(.)) %>% 
  melt(id = "isolate") %>% as_tibble() %>% 
  set_colnames(c("isolate", "tp2_h", "tp2_cl")) %>% 
  factorToInt("tp2_h") %>% 
  createID(., "tp2_h", "tp2_cl")

# t1_comps <- t1_coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>% 
#   arrange(tp_h, tp_cl, isolate) %>% 
#   compsSet(., 1)
# # saveRDS(t1_comps, "t1_comps.Rds")
t1_comps <- readRDS("t1_comps.Rds")

# t2_comps <- t2_coded %>% rename(tp_h = tp2_h, tp_cl = tp2_cl) %>% 
#   arrange(tp_h, tp_cl, isolate) %>% 
#   compsSet(., 2)
# 
# novels <- setdiff(t2_coded$isolate, t1_coded$isolate)
# counting_novels <- t2_coded %>% 
#   filter(isolate %in% novels) %>% 
#   group_by(id) %>% 
#   summarise(num_novs = n(), .groups = "drop")
# t2_comps <- t2_comps %>% left_join(., counting_novels, by = "id")
# t2_comps$num_novs[is.na(t2_comps$num_novs)] <- 0
# # saveRDS(t2_comps, "t2_comps.Rds")
t2_comps <- readRDS("t2_comps.Rds")

# --------------------------------------------------------------------------------------------------------
# part2 <- collectGrowth(changed_clusters, h_before, t1_raw, t2_raw)
# hxallc <- changed_clusters
trackClusters <- function(hx, t2_comps) {
  pb <- txtProgressBar(min = 0, max = nrow(hx), initial = 0, style = 3)
  k <- 0
  tmp <- lapply(1:nrow(hx), function(i) {
    # print(paste0(i, "/", nrow(hxallc)))
    k <- k + i
    setTxtProgressBar(pb, k)
    # the composition (and cluster size) data for the cluster in row i of 
    # the TP1 composition data (for the height in hx)
    cxdata <- hx[i,]
    
    # the indices of the TP2 composition dataset where clusters contain at least all 
    # the isolates from the cxdata cluster
    inds <- grep(cxdata$composition, t2_comps$composition)
    kc <- t2_comps[inds,] %>% select(-composition, -id)
    
    # merged the data for the single TP1 clusters and all corresponding TP2 clusters
    cxdata %>% select(-composition, -id) %>% bind_cols(., kc) %>% return()
  }) %>% bind_rows()
  close(pb)
  return(tmp)
}

h_before <- 0
print(paste0("Collecting data for height ", h_before))
hx <- t1_comps %>% filter(tp1_h == h_before)
# tmp <- trackClusters(hx, t2_comps)
# saveRDS(tmp, "tmp.Rds")
tmp <- readRDS("tmp.Rds")

# now need to make metric calculations and extract the best result:
#     - the TP2 cluster that has the largest proportion of novels to growth rate
#     i.e. the "best" TP2 cluster where we look at clusters that absorbed lots of novels very quickly
# the actual growth of these clusters (TP2 size - TP1 size) / TP1 size
tmp$actual_growth_rate <- (tmp$tp2_cl_size - tmp$tp1_cl_size) / tmp$tp1_cl_size
# the number of novels in the TP2 cluster over the growth rate --> growth acceleration?
tmp$acc <- tmp$num_novs / tmp$actual_growth_rate
tmp$acc[is.na(tmp$acc)] <- 0
tmp$flag <- NA


# note: the way these are flagged, we can see how many TP2 clusters were formed from the 
# selected TP1 cluster before we get to the "key growth" where a cluster went from its state
# at TP1 to its state at TP2 with a remarkable amount of growth
pb <- txtProgressBar(min = 0, max = nrow(hx), initial = 0, style = 3)
clx <- unique(hx$tp1_cl)
k <- 0
h0 <- lapply(1:length(clx), function(i) {
  k <- k + i
  setTxtProgressBar(pb, k)
  df <- tmp %>% filter(tp1_cl == clx[i]) %>% arrange(tp2_h, tp2_cl)
  df$flag <- 1:nrow(df)
  df %>% arrange(-acc, tp2_h, tp2_cl) %>% slice(1) %>% return()
}) %>% bind_rows()
close(pb)

h0$acc[which(h0$actual_growth_rate == 0)] <- 0

tmp2 <- h0


# comps_before <- lapply(unique(cb$tp1_cl), function(j) {
#   cb %>% filter(tp1_cl == j) %>% pull(isolate) %>% sort() %>%
#     paste0(collapse = ",") %>% tibble(composition = .) %>%
#     bind_cols(tp1_h = h_before, tp1_cl = j)
# }) %>% bind_rows() %>%
#   set_colnames(c("composition", "h_before", "cl_before"))


for (h in unique(t1_coded$tp1_h)[-1]) {
  
  h_after <- h
  print(paste0("Collecting data for height ", h_after))
  
  ca <- t1_coded %>% filter(tp1_h == h_after)
  comps_after <- lapply(unique(ca$tp1_cl), function(j) {
    ca %>% filter(tp1_cl == j) %>% pull(isolate) %>% sort() %>% 
      paste0(collapse = ",") %>% tibble(composition = .) %>% 
      bind_cols(tp1_h = h_after, tp1_cl = j)
  }) %>% bind_rows() %>% 
    set_colnames(c("composition", "h_after", "cl_after"))

  # clusters that have not changed from the previous height
  did_not_change <- intersect(comps_after$composition, comps_before$composition)
  stayed_the_same <- comps_after %>% filter(composition %in% did_not_change) %>% 
    left_join(., comps_before, by = "composition")

  part1 <- stayed_the_same %>% 
    left_join(., tmp3, by = c("h_before" = "tp1_h", "cl_before" = "tp1_cl")) %>% 
    select(-h_before, -cl_before, -composition) %>% 
    rename(tp1_h = h_after, tp1_cl = cl_after)
  
  # clusters that did change from the previous height
  changed_clusters <- comps_after %>% 
    filter(!(composition %in% did_not_change)) %>% 
    select(-composition) %>% set_colnames(c("tp1_h", "tp1_cl"))

  part2 <- collectGrowth(changed_clusters, t1_raw, t2_raw)
  part2$flag <- NA
  
  for (j in 1:nrow(part2)) {
    if (part2$id[j] %in% ids) {
      part2$flag[j] <- "sb"
    }else {
      ids <- c(unlist(ids), part2$id[j])
      part2$flag[j] <- part2$id[j]
    }
  }
  
  h_before <- h_after
  comps_before <- comps_after %>% set_colnames(c("composition", "h_before", "cl_before"))
  
  tmp3 <- bind_rows(part1, part2) %>% bind_rows(tmp3, .)
  saveRDS(tmp3, "tmp3.Rds")
}



# saveRDS(ids, "ids.Rds")
# saveRDS(tmp3, "tmp4.Rds")








# these are the TP2 clusters that each have at least some isolates 
# from cluster x and also have size at least as large as cluster x  
# 
# these are the isolates found in the TP2 clusters that contain 
# at least all the isolates from cluster x and possibly more
# 
# these are the actual sizes of these TP2 clusters
# 
# a3 is the set of height-clusters and their actual TP2 cluster sizes that 
# contain all of the isolates from hx-cx (and possibly more isolates)


# ggplot(a5, aes(x = tp2_h, y = acc)) + geom_point()























