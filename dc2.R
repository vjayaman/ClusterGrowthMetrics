x <- c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress", "stringr")
lapply(x, require, character.only = TRUE)

# --------------------------------------------------------------------------------------------------------
## FOR USER: replace the filename variables with quoted file paths if you don't want to input them each time
# # the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||33333
# # the raw datasets, no filtering or other changes made

message("Reading timepoint 1 cluster assignment, processing...\n")
t1_coded <- "data/timepoint1_data.csv" %>% readData(., 1) %>% codeIsolates(., "tp1")
t1_comps <- t1_coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>% compsSet(., 1)

message("Reading timepoint 2 cluster assignment, processing...\n")
t2_coded <- "data/timepoint2_data.csv" %>% readData(., 2) %>% codeIsolates(., "tp2")

novels <- setdiff(t2_coded$isolate, t1_coded$isolate)
counting_novels <- t2_coded %>%
  filter(isolate %in% novels) %>%
  group_by(id) %>%
  summarise(num_novs = n(), .groups = "drop")

t2_comps <- t2_coded %>% 
  rename(tp_h = tp2_h, tp_cl = tp2_cl) %>% 
  compsSet(., 2) %>% 
  left_join(., counting_novels, by = "id")
t2_comps$num_novs[is.na(t2_comps$num_novs)] <- 0

# --------------------------------------------------------------------------------------------------------
# saveData(dtype=2, tp=1, tpdata=t1_comps)
# t1_comps <- readRDS("t1_comps.Rds")

# saveData(dtype=2, tp=2, tpdata=t2_comps)
# t2_comps <- readRDS("t2_comps.Rds")
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
hbdata <- lapply(1:length(clx), function(i) {
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

alldata <- hbdata

saveData(dtype=1, a=alldata, hbdata=hbdata, hb=h_before, cb=bef_comps)
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
  ss <- noChange(aft_comps, bef_comps, hbdata)
  
  print("part2")
  hdata <- aft_comps %>% 
    filter(!(id_aft %in% ss$id)) %>% 
    set_colnames(colnames(t1_comps))
  cc <- trackClusters(hdata, t2_comps) %>% 
    add_column(actual_growth_rate = NA, acc = NA, flag = NA)
  
  print("part3")
  if (nrow(cc) > 0) {
    cc$actual_growth_rate <- (cc$tp2_cl_size - cc$tp1_cl_size) / cc$tp1_cl_size
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
  hbdata <- hnew
  h_before <- h_after
  bef_comps <- aft_comps %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
  
  saveData(dtype=1, a=alldata, hbdata=hbdata, hb=h_before, cb=bef_comps)
}

# restarted process at 8:22, starting height 5 at 8:42, height 11 at 8:47, height 428 at 9:18, height 630 at 9:29, 
# height 893 at 9:50, finished at 9:53
# AN HOUR AND A HALF!!!
