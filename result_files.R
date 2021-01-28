
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

oneHeight <- function(h_i, novels) {
  oneh <- readRDS(paste0("outputs/height_data/h", h_i, ".Rds")) %>%
    left_join(., dfx, by = c("flag" = "first_flag")) %>%
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)

  # reading in the growth data for a particular height, and arranging by first TP1 then by TP2
  # extracting the first TP2_h_c match for each TP1 cluster at this height, and adding a column
  # indicating the actual size change (which may include non-novels)
    c1 <- oneh %>% pull(id) %>% unique() %>% 
      lapply(., function(idx) {
        a1 <- oneh %>% filter(id == idx) %>% slice(1) %>% 
          rename(tp1_id = id) %>% createID(., "tp2", "tp2_h", "tp2_cl") %>% rename(tp2_id = id)
        additionalTP1(b1, b2, a1$tp1_id, a1$tp2_id, novels) %>% 
          left_join(a1, ., by = c("tp1_id", "tp2_id")) %>% return()
      }) %>% bind_rows()
  
  tp1_size <- c1$tp1_cl_size; tp2_size <- c1$tp2_cl_size; novels_hx <- c1$num_novs

  c2 <- c1 %>% add_column(
      actual_size_change = tp2_size - tp1_size, 
      actual_growth_rate = (tp2_size - tp1_size) / tp1_size,
      # actual_growth_rate = tp2_size / tp1_size, 
      new_growth = tp2_size / (tp2_size - novels_hx)
      # b_ov_growth = novels_hx / ((tp2_size - tp1_size) / tp1_size)
      (tp1_size + novels_hx) / tp1_size
      ) %>% 
    
    rename(additional_tp1 = add_tp1) %>% 
    select(tp1_id, tp1_h, tp1_cl, tp1_cl_size, flag, last_flag, 
           tp2_id, tp2_h, tp2_cl, tp2_cl_size, additional_tp1, num_novs, actual_size_change, 
           actual_growth_rate, b_ov_growth, new_growth)
  c2$b_ov_growth[is.na(c2$b_ov_growth)] <- 0
  

  c3 <- c2 %>% arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
    leadingZeros(., "tp1_cl", "c") %>% leadingZeros(., "tp2_cl", "c") %>%
    leadingZeros(., "tp1_h", "h", w = max(nchar(colnames(time1_raw)[-1]))) %>%
    leadingZeros(., "tp2_h", "h", w = max(nchar(colnames(time2_raw)[-1])))

  return(c3)
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

t1_comps <- c(time1_raw$isolate, time2_raw$isolate) %>% unique() %>% as_tibble() %>%
  set_colnames("char_isolate") %>% rowid_to_column("num_isolate") %>%
  codeIsolates(time1_raw, "tp1", .) %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>%
  compsSet(., "TP1", indicate_progress = FALSE) %>% arrange(tp1_h, tp1_cl)

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
  
  c4 <- lapply(pt$Threshold, function(hval) {oneHeight(hval, novels)})
}

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
  
  