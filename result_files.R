
suppressWarnings(suppressPackageStartupMessages(source("funcs.R")))

option_list <- list(
  
  make_option(c("-a", "--tp1"), metavar = "file", default = NULL,
              help = "Time point 1 file name"),
  
  make_option(c("-b", "--tp2"), metavar = "file", default = NULL,
              help = "Time point 2 file name"),
  
  make_option(c("-t", "--type"), metavar = "character", default = NULL,
              help = paste0("Collection type: 'nawc' to use local minima ",
                            "heights of stability plot, otherwise leave empty ",
                            "and provide a particular height in the next argument")), 
  
  make_option(c("-x", "--height"), metavar = "numeric", default = NULL, 
              help = paste0("A number, e.g. 50, to use as a particular height for ", 
                            "which to generate cluster and strain tables"))
)

arg <- parse_args(OptionParser(option_list=option_list))

oneHeightV2 <- function(h_i, novels, t2_composition, t1_composition) {
  oneh <- readRDS(paste0("outputs/height_data/h", h_i, ".Rds")) %>%
    left_join(., dfx, by = c("flag" = "first_flag")) %>%
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)

  # check: (no cluster numbers skipped) - none should be skipped by definition, but just in case
  # identical(unique(oneh$tp1_cl), min(oneh$tp1_cl):max(oneh$tp1_cl))
  
  # first TP2 match for each TP1 cluster at this height
  d2 <- oneh[diff(c(0, oneh$tp1_cl)) != 0,] %>% rename(tp1_id = id) %>% 
    createID(., "tp2", "tp2_h", "tp2_cl") %>% rename(tp2_id = id)
  
  matched <- d2 %>% select(tp1_id, tp2_id)
  
  e2 <- t2_composition %>% select(id, composition) %>% rename(tp2_id = id, comp2 = composition)
  
  e1 <- t1_composition %>% select(id, composition) %>% rename(tp1_id = id, comp1 = composition) %>% 
    filter(tp1_id %in% matched$tp1_id) %>% 
    left_join(., matched, by = "tp1_id") %>% 
    left_join(., e2, by = "tp2_id")
  
  did_not_change <- e1 %>% filter(comp1 == comp2) %>% select(tp1_id, tp2_id) %>% add_column(add_tp1 = 0)
  chg <- e1 %>% filter(comp1 != comp2) %>% select(tp1_id, tp2_id)
  
  c1 <- lapply(1:nrow(chg), function(i) {
    additionalTP1(b1, b2, chg$tp1_id[i], chg$tp2_id[i], novels)
  }) %>% bind_rows(., did_not_change) %>% left_join(d2, ., by = c("tp1_id", "tp2_id"))
  
  c1 %>% 
    add_column(
      actual_size_change = (c1$tp2_cl_size - c1$tp1_cl_size), 
      actual_growth_rate = ((c1$tp2_cl_size - c1$tp1_cl_size) / c1$tp1_cl_size) %>% round(., digits = 3),
      new_growth = (c1$tp2_cl_size / (c1$tp2_cl_size - c1$num_novs)) %>% round(., digits = 3)) %>% 
    rename(additional_tp1 = add_tp1) %>% 
    select(tp1_id, tp1_h, tp1_cl, tp1_cl_size, flag, last_flag, tp2_id, tp2_h, tp2_cl, tp2_cl_size, 
           additional_tp1, num_novs, actual_size_change, actual_growth_rate, new_growth) %>% 
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
    leadingZeros(., "tp1_cl", "c") %>% leadingZeros(., "tp2_cl", "c") %>%
    leadingZeros(., "tp1_h", "h", w = max(nchar(colnames(time1_raw)[-1]))) %>%
    leadingZeros(., "tp2_h", "h", w = max(nchar(colnames(time2_raw)[-1]))) %>% return()
}

# time1_raw <- readBaseData("data/timepoint1_data.csv", 1)
# time2_raw <- readBaseData("data/timepoint2_data.csv", 2)
# data_type <- "nawc"

if (!is.null(arg$tp1)) {time1_raw <- readBaseData(arg$tp1, 1)}
if (!is.null(arg$tp2)) {time2_raw <- readBaseData(arg$tp2, 2)}
if (!is.null(arg$type)) {data_type <- arg$type}else {data_type <- NULL}
if (!is.null(arg$height)) {height_choice <- arg$height}else {height_choice <- NULL}

if (!exists("time1_raw") | !exists("time2_raw")) {
  stop("No input successfully collected.")
}

novels <- setdiff(time2_raw$isolate, time1_raw$isolate)
heights <- colnames(time1_raw)[-1]

all_isolates <- c(time1_raw$isolate, time2_raw$isolate) %>% unique() %>% 
  as_tibble() %>% set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

t1_coded <- time1_raw %>% codeIsolates(., "tp1", all_isolates)
t2_coded <- time2_raw %>% codeIsolates(., "tp2", all_isolates)

t2_colnames <- t2_coded$tp2_h %>% unique() %>% sort()

t1_comps <- t1_coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>%
  compsSet(., "TP1", indicate_progress = TRUE)

novels_coded <- setdiff(t2_coded$isolate, t1_coded$isolate)
counting_novels <- t2_coded %>% filter(isolate %in% novels_coded) %>%
  group_by(id) %>% summarise(num_novs = n(), .groups = "drop")

t2_comps <- t2_coded %>% rename(tp_h = tp2_h, tp_cl = tp2_cl) %>%
  compsSet(., "TP2", indicate_progress = TRUE) %>% left_join(., counting_novels, by = "id")
t2_comps$num_novs[is.na(t2_comps$num_novs)] <- 0

# Identifying the last time each cluster was seen
allcomps <- t1_comps$composition %>% unique()
dfx <- lapply(1:length(allcomps), function(i) {
  t1_comps %>% filter(composition == allcomps[i]) %>% slice(1, nrow(.)) %>%
    pull(id) %>% t() %>% data.frame(stringsAsFactors = FALSE)
}) %>% bind_rows() %>% as_tibble() %>% set_colnames(c("first_flag", "last_flag"))


b1 <- time1_raw %>% melt(id = "isolate") %>% as_tibble() %>% rename(tp1_cl = value, tp1_h = variable) %>% 
  mutate(across(tp1_h, as.character)) %>% createID(., "tp1", "tp1_h", "tp1_cl") %>% rename(tp1_id = id)

b2 <- time2_raw %>% melt(id = "isolate") %>% as_tibble() %>% rename(tp2_cl = value, tp2_h = variable) %>% 
  mutate(across(tp2_h, as.character)) %>% createID(., "tp2", "tp2_h", "tp2_cl") %>% rename(tp2_id = id)

if (is.null(data_type)) {
  
  if (is.numeric(height_choice)) {
    c4 <- oneHeight(height_choice, novels)
  }else {
    print("Not numeric")
  }
  
}else if (data_type == "nawc") {
  
  a2 <- read.table("outputs/summary/stats.tsv", sep = "\t", header = TRUE) %>% as_tibble()
  a2$Neighbour.AWC[is.na(a2$Neighbour.AWC)] <- 0
  
  x <- a2$Neighbour.AWC
  # differences x[i+1] - x[i]
  y4 <- which(diff(sign(diff(x))) > 0) + 1
  # find the points that shrank from the left and grew to the right, note the indexing starts from 0
  
  pt <- a2[y4,] %>% filter(Neighbour.AWC < 0.95)
  
  if (length(colnames(a2)) < 4) {
    a2$Selected.for.Analysis <- 0
    a2$Selected.for.Analysis[a2$Threshold %in% pt$Threshold] <- 1
    write.table(a2, file = file.path("outputs/summary/", "stats_with_analysis.tsv"), 
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  p <- ggplot(a2, aes(x = Threshold, y = Neighbour.AWC)) +
    geom_step() +
    geom_point(data = pt, aes(x = Threshold, y = Neighbour.AWC), color = "red") +
    labs(x = "goeBURST Threshold", y = "")
  
  ggsave("outputs/summary/nawc_stability_plot.png",
         device = "png", plot = p, width = 16, height = 9, units = "in", dpi = 320)
  print(pt$Threshold)
  
  Sys.time()
  c5 <- lapply(pt$Threshold, function(hval) {oneHeightV2(hval, novels)}) %>% bind_rows()
  # saveRDS(c5, "tmp1.Rds")
  Sys.time()
  
  # # 7:50-8:24
  # c4 <- lapply(pt$Threshold[10:20], function(hval) {
  #   print(hval)
  #   oneHeight(hval, novels)}) %>% bind_rows()
  # # saveRDS(c4, "tmp2.Rds")
  # Sys.time()
}

# a1 <- c4 %>% arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)
# a1$actual_growth_rate %<>% round(., digits = 3)
# a1$new_growth %<>% round(., digits = 3)
# 
# a2 <- c5 %>% arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)
# 
# a2$actual_growth_rate %<>% round(., digits = 3)
# a2$actual_size_change %<>% as.integer()
# a2$new_growth %<>% round(., digits = 3)
# a2$additional_tp1 %<>% as.integer()
# 
# identical(a1, a2)

clusters_formatted <- c4 %>% set_colnames(c(
    "TP1 ID", "TP1 height", "TP1 cluster", "TP1 cluster size",
    "First time this cluster was seen in TP1", "Last time this cluster was seen in TP1",
    "First time this cluster was seen in TP2", "TP2 height", "TP2 cluster", "TP2 cluster size",
    "Number of additional TP1 strains in the TP2 match", "Number of novels in the TP2 match",
    "Actual cluster size change (TP2 size - TP1 size)",
    "Actual growth rate = (TP2 size - TP1 size) / (TP1 size)",
    "Number of novels / Actual growth rate",
    "Novel growth = (TP2 size) / (TP2 size - number of novels)"))

write.csv(clusters_formatted, paste0("outputs/summary/TP1_cluster_results.csv"), row.names = FALSE)
  

isolates_formatted <- time1_raw %>%
  select(isolate, as.character(pt$Threshold)) %>%
  melt(., id = "isolate") %>% as_tibble() %>%
  rename(tp1_h = variable, tp1_cl = value) %>%
  mutate(across(tp1_h, as.character)) %>%
  mutate(across(tp1_h, as.integer)) %>%
  leadingZeros(., "tp1_h", "h") %>%
  leadingZeros(., "tp1_cl", "c") %>%
  right_join(., c4, by = c("tp1_h", "tp1_cl")) %>%
  arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>%
  select(isolate, colnames(c4)) %>%
  set_colnames(c("Isolates", colnames(clusters_formatted)))

write.csv(isolates_formatted, paste0("outputs/summary/TP1_strain_results.csv"), row.names = FALSE)

