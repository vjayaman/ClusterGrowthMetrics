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
t2_raw <- time2_raw

hxallc <- t1clusters %>% filter(tp1_h == 0)



# h0 <- lapply(1:nrow(hxallc), function(i) {
for (i in 1:nrow(hxallc))  {
  print(paste0(i, "/", nrow(hxallc)))
  hx <- hxallc$tp1_h[i]
  clx <- hxallc$tp1_cl[i]
  
  t1_hxcx <- t1_raw %>% 
    filter(tp1_h == hx & tp1_cl == clx) %>% 
    add_column(tp1_cl_size = nrow(.))
  
  # these are the TP2 clusters that each have at least some isolates 
  # from cluster x and also have size at least as large as cluster x  
  a3 <- t2_raw %>% 
    filter(isolate %in% t1_hxcx$isolate) %>% 
    createID(., "tp2_h", "tp2_cl")
    
  a4 <- table(a3$id) %>% as.data.frame() %>% as_tibble() %>% 
    filter(Freq >= nrow(t1_hxcx)) %>% select(Var1)
  
  a3 %>% filter(id %in% a4$Var1) %>% 
  countRows(., "num_from_cx") %>% 
    filter(num_from_cx >= nrow(t1_hxcx)) %>% 
    # these are the isolates found in the TP2 clusters that contain 
    # at least all the isolates from cluster x and possibly more
    select(-num_from_cx) %>% 
    right_join(t2_raw, ., c("tp2_h", "tp2_cl")) %>% 
    # these are the actual sizes of these TP2 clusters
    countRows(., "tp2_cl_size")
  # a3 is the set of height-clusters and their actual TP2 cluster sizes that 
  # contain all of the isolates from hx-cx (and possibly more isolates)
  
  a4 <- left_join(a3, withnovs, by = c("tp2_h", "tp2_cl"))
  a4$num_novs[is.na(a4$num_novs)] <- 0
  
  a5 <- t1_hxcx %>% 
    select(-isolate) %>% 
    unique() %>% 
    bind_cols(., a4)
  
  a5$growth <- (a5$tp2_cl_size - a5$tp1_cl_size) / a5$tp1_cl_size
  a5$acc <- a5$num_novs / a5$growth
  
  if (i == 1) {
    tmp2 <- a5[which.max(a5$acc),]
  }else {
    tmp2 <- a5[which.max(a5$acc),] %>% bind_rows(tmp2, .)
  }
  
}
# })




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

tmp2 <- readRDS("tmp2.Rds")
















