#! /usr/bin/env Rscript

suppressWarnings(suppressPackageStartupMessages(source("functions/tracking_functions.R")))
source("class_definitions.R")

option_list <- list(
  make_option(c("-a", "--tp1"), metavar = "file", default = NULL, help = "Time point 1 file name"),
  make_option(c("-b", "--tp2"), metavar = "file", default = NULL, help = "Time point 2 file name"),
  # make_option(c("-d", "--delimiter"), metavar = "character", default = "\t"), 
  make_option(c("-x", "--heights"), metavar = "character", default = NULL,
              help = paste0("A string of comma-delimited numbers, e.g. '50,75,100' to ", 
                            "use as heights for which to generate cluster and strain tables")))

arg <- parse_args(OptionParser(option_list=option_list))

stopwatch <- list(as.character.POSIXt(Sys.time()))

# f1 <- readBaseData(arg$tp1, 1, "\t")#arg$delimiter)
# f2 <- readBaseData(arg$tp2, 2, "\t")#arg$delimiter)
# heights <- strsplit(arg$heights, split = ",") %>% unlist()

f1 <- readBaseData("t1_clusters_processed.csv", 1, "\t")
f2 <- readBaseData("t2_clusters_processed.csv", 2, "\t")

colnames(f1)[1] <- colnames(f2)[1] <- "isolate"

all_isolates <- unique(c(f1$isolate, f2$isolate)) %>% as_tibble() %>% 
  set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

tp1 <- Timedata$new("tp1", raw = f1, all_isolates)
tp2 <- Timedata$new("tp2", raw = f2, all_isolates)

tp1$coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>% tp1$set_comps(.)

t2_colnames <- tp2$coded$tp2_h %>% unique() %>% sort()

novels <- setdiff(tp2$coded$isolate, tp1$coded$isolate)
counting_novels <- tp2$coded %>% filter(isolate %in% novels) %>%
  group_by(id) %>% summarise(num_novs = n(), .groups = "drop")

tp2$coded %>% rename(tp_h = tp2_h, tp_cl = tp2_cl) %>% tp2$set_comps()
tp2$comps <- tp2$comps %>% left_join(., counting_novels, by = "id")
tp2$comps$num_novs[is.na(tp2$comps$num_novs)] <- 0


collected_data <- read.csv(file = "outputs/CGM_strain_results.txt", sep = "\t",
                           stringsAsFactors = FALSE, numerals = "no.loss") %>% as_tibble() %>% 
  select(TP1.ID, TP1.height, TP1.cluster, First.time.this.cluster.was.seen.in.TP2, TP2.height, 
         TP2.cluster, TP1.cluster.size, TP2.cluster.size, Number.of.additional.TP1.strains.in.the.TP2.match, 
         Number.of.novels.in.the.TP2.match, Actual.cluster.size.change..TP2.size...TP1.size.) %>% 
  set_colnames(c("tp1_id", "tp1_h", "tp1_cl", "tp2_id", "tp2_h", "tp2_cl", 
                 "tp1_cl_size", "tp2_cl_size", "add_TP1", "num_nov", "size_change")) %>% 
  unique() %>% 
  arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)

collected_data$tp1_h %<>% charToInt(., "h")
collected_data$tp1_cl %<>% charToInt(., "c")
collected_data$tp2_h %<>% charToInt(., "h")
collected_data$tp2_cl %<>% charToInt(., "c")

heights <- collected_data %>% filter(!is.na(tp1_h)) %>% pull(tp1_h) %>% unique()

b1 <- meltedSizing(tp1$raw, "tp1", padding_heights, padding_clusters)
b2 <- meltedSizing(tp2$raw, "tp2", padding_heights, padding_clusters)
novels <- setdiff(b2$isolate, b1$isolate)

for (j in 1:length(heights)) {
  cd <- collected_data %>% filter(tp1_h == heights[j]) %>% 
    select(tp1_id, tp1_h, tp1_cl, tp2_id, tp2_h, 
           tp2_cl, tp1_cl_size, tp2_cl_size, add_TP1, num_nov, size_change)
  
  hx <- b1 %>% filter(tp1_h == heights[j])
  b1clusters <- hx$tp1_cl %>% unique() %>% sort()
  
  pb <- txtProgressBar(min = 0, max = length(b1clusters), initial = 0, style = 3)
  ad <- lapply(1:length(b1clusters), function(i) {
    setTxtProgressBar(pb, i)
    # cluster assignments for cluster b1clusters[i] in TP1 at heights[j]
    c1 <- hx %>% filter(tp1_cl == b1clusters[i])
    # tp1_csize <- c1$tp1_cl_size %>% unique() # will be length 1
    
    # what cluster(s) do these isolates belong to in TP2
    basecase <- b2 %>% filter(isolate %in% c1$isolate)
    
    # how many of the original TP1 isolates from b1clusters[i] make up these TP2 clusters?
    sizes <- basecase %>% 
      group_by(tp2_id) %>% 
      summarise(num_from_c1 = n(), .groups = "drop") %>% 
      left_join(basecase, ., by = "tp2_id")
    
    # want TP2 clusters that *at least* contain all the b1clusters[i] isolates
    compared <- sizes %>% filter(num_from_c1 >= c1$tp1_cl_size[1]) %>% 
      select(-isolate) %>% unique() %>% 
      arrange(tp2_h, tp2_cl) %>% 
      slice(1)
    
    compared$num_nov <- b2 %>% 
      filter(tp2_id == compared$tp2_id) %>% 
      filter(isolate %in% novels) %>% 
      nrow()
    
    original_and_novel <- c(c1$isolate, novels)
    compared$add_TP1 <- b2 %>% 
      filter(tp2_id == compared$tp2_id) %>% 
      filter(!(isolate %in% original_and_novel)) %>% 
      nrow()

    c1 %>% select(tp1_id, tp1_h, tp1_cl, tp1_cl_size) %>% slice(1) %>% 
      bind_cols(., compared) %>% 
      add_column(size_change = compared$tp2_cl_size - c1$tp1_cl_size[1]) %>% 
      select(tp1_id, tp1_h, tp1_cl, tp2_id, tp2_h, tp2_cl, tp1_cl_size, 
             tp2_cl_size, add_TP1, num_nov, size_change) %>% return()
    
  }) %>% bind_rows() %>% arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
    mutate(across(tp2_cl_size, as.integer))
  close(pb)
  
  print(paste0("Threshold ", heights[j], " (", j, " of ", length(heights), "): ", identical(cd, ad)))
}

stopwatch <- append(stopwatch, values = as.character.POSIXt(Sys.time())) %>% 
  set_names(c("start_time", "end_time"))

