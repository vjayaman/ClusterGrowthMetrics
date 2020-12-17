x <- c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress")
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

# source("functions/base_functions.R")
# source("functions/processing_functions.R")

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
  factorToInt("tp1_h")

t2_raw <- time2_raw
t2_raw$id <- paste0(t2_raw$tp2_h, "-", t2_raw$tp2_cl)

collectGrowth <- function(hxallc, height, t1_raw, t2_raw) {
  print(height)
  
  lapply(1:nrow(hxallc), function(i) {
    print(paste0(i, "/", nrow(hxallc)))
    hx <- hxallc$tp1_h[i]
    clx <- hxallc$tp1_cl[i]
    
    t1_hxcx <- t1_raw %>% 
      filter(tp1_h == hx & tp1_cl == clx) %>% 
      add_column(tp1_cl_size = nrow(.))
    
    cxdata <- t1_hxcx %>% select(tp1_h, tp1_cl, tp1_cl_size) %>% slice(1)
    
    a3 <- t2_raw %>% filter(isolate %in% t1_hxcx$isolate)
    
    kc <- table(a3$id) %>% `>=`(nrow(t1_hxcx)) %>% which() %>% names()
    
    a5 <- t2_raw %>% filter(id %in% kc) %>% 
      select(id) %>% table() %>% as.data.frame() %>% as_tibble() %>% 
      set_colnames(c("id", "tp2_cl_size")) %>% 
      mutate(across(id, as.character))
    
    a9 <- a3 %>% select(-isolate) %>% unique() %>% 
      right_join(., a5, by = "id") %>% 
      left_join(., withnovs, by = c("tp2_h", "tp2_cl")) %>% 
      mutate(num_novs = ifelse(is.na(num_novs), 0, num_novs)) %>% 
      bind_cols(cxdata, .)
    
    a9$growth <- (a9$tp2_cl_size - a9$tp1_cl_size) / a9$tp1_cl_size
    a9$acc <- a9$num_novs / a9$growth
    
    a9[which.max(a9$acc),] %>% return()
  }) %>% bind_rows() %>% return()
}

h_before <- 0
hxallc <- t1clusters %>% filter(tp1_h == h_before)
# tmp2 <- collectGrowth(hxallc, h_before, t1_raw, t2_raw)
tmp2 <- readRDS("tmp2.Rds")


h_before <- 0
for (h in unique(t1_coded$tp1_h)[-1]) {
  h_after <- h
  ca <- t1_coded %>% filter(tp1_h == h_after)
  comps_after <- lapply(unique(ca$tp1_cl), function(j) {
    ca %>% filter(tp1_cl == j) %>% pull(isolate) %>% sort() %>% 
      paste0(collapse = ",") %>% tibble(composition = .) %>% 
      bind_cols(tp1_h = h_after, tp1_cl = j)
  }) %>% bind_rows() %>% 
    set_colnames(c("composition", "h_after", "cl_after"))

  cb <- t1_coded %>% filter(tp1_h == h_before)
  comps_before <- lapply(unique(cb$tp1_cl), function(j) {
    cb %>% filter(tp1_cl == j) %>% pull(isolate) %>% sort() %>% 
      paste0(collapse = ",") %>% tibble(composition = .) %>% 
      bind_cols(tp1_h = h_before, tp1_cl = j)
  }) %>% bind_rows() %>% 
    set_colnames(c("composition", "h_before", "cl_before"))

  # clusters that have not changed from the previous height
  did_not_change <- intersect(comps_after$composition, comps_before$composition)
  stayed_the_same <- comps_after %>% filter(composition %in% did_not_change) %>% 
    left_join(., comps_before, by = "composition")

  part1 <- stayed_the_same %>% 
    left_join(., tmp2, by = c("h_before" = "tp1_h", "cl_before" = "tp1_cl")) %>% 
    select(-h_before, -cl_before, -composition) %>% 
    rename(tp1_h = h_after, tp1_cl = cl_after)
  
  
  # clusters that did change from the previous height
  changed_clusters <- comps_after %>% filter(!(composition %in% did_not_change)) %>% 
    select(-composition) %>% set_colnames(c("tp1_h", "tp1_cl"))

  h_before <- h_after
  part2 <- collectGrowth(changed_clusters, h_before, t1_raw, t2_raw)
  tmp2 <- bind_rows(tmp2, new_height)
  saveRDS(tmp2, "tmp3.Rds")
}

tmp3 <- readRDS("tmp3.Rds")













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









# notin <- setdiff(time1_raw$isolate, t1_h_cx$isolate)
# ids <- list()
# i <- 1
# hx <- a3$height[i]
# cx <- a3$cluster[i]
# x <- time2_raw %>% filter(height == hx & cluster == cx) %>% pull(isolate)
# x1 <- list(x)
# ids[paste0(hx, "-", cx)] <- x1
# tmp2 <- a3 %>% filter(height == hx & cluster == cx) %>% 
#   add_column(fo = length(intersect(x, t1_h_cx$isolate))) %>% 
#   add_column(novs = length(intersect(x, novels))) %>% 
#   add_column(o = length(intersect(x, notin))) %>% 
#   add_column(type = "new")
# 
# for (i in 2:nrow(a2)) {
#   print(i)
#   hx <- a3$height[i]
#   cx <- a3$cluster[i]
#   x <- time2_raw %>% filter(height == hx & cluster == cx) %>% pull(isolate)
#   
#   x1 <- list(x)
#   
#   if (x1 %in% ids) {
#     id <- names(ids)[which(x1 %in% ids)] %>% 
#       strsplit(., split = "-") %>% unlist() %>% as.integer()
#     tmp <- tmp2 %>% filter(height == id[1] & cluster == id[2]) %>% 
#       mutate(height = hx, cluster = cx, type = "old")
#   }else {
#     ids[paste0(hx, "-", cx)] <- x1
#     tmp <- a3 %>% filter(height == hx & cluster == cx) %>% 
#       add_column(fo = length(intersect(x, t1_h_cx$isolate))) %>% 
#       add_column(novs = length(intersect(x, novels))) %>% 
#       add_column(o = length(intersect(x, notin))) %>% 
#       add_column(type = "new")
#   }
#   tmp2 <- bind_rows(tmp2, tmp)
# }

# tmp2 <- readRDS("tmp2.Rds")
















