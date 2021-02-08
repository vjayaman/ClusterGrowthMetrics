#! /usr/bin/env Rscript

# Note: a line of stars (********) will be placed after verified complicated operations, 
# involved statements that do what they are expected to do, and this has been checked for sure

# SUMMARY: this is the script for generating the cluster metrics, given input files for two time points
# 	- each time point dataset should be a file that can be read with the following statement 
# 	      read.csv(file = tp1_filename, stringsAsFactors = FALSE, numerals = "no.loss")
# 	  and should be composed of a single column with name "isolate" followed by a series of columns with 
# 	  cluster assignments at each of the thresholds

# sending errors and warnings to the log file
msg <- file("outputs/logfile_datacollection.txt", open="wt")
sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("funcs.R")))

option_list <- list(
  make_option(c("-a", "--tp1"), metavar = "file", default = NULL, help = "Time point 1 file name"),
  make_option(c("-b", "--tp2"), metavar = "file", default = NULL, help = "Time point 2 file name"),
  # make_option(c("-t", "--type"), metavar = "character", default = NULL,
  #             help = paste0("Collection type: 'nawc' to use local minima heights of stability ", 
  #                           "plot, otherwise leave empty and provide a particular height in ", 
  #                           "the next argument")), 
  make_option(c("-x", "--heights"), metavar = "character", default = NULL,
              help = paste0("A string of comma-delimited numbers, e.g. '50,75,100' to ", 
                            "use as heights for which to generate cluster and strain tables"))
)

arg <- parse_args(OptionParser(option_list=option_list))

stopwatch <- rep(0,2) %>% set_names(c("start_time", "end_time"))
stopwatch[1] <- Sys.time()

outputDetails(paste0("\n||", paste0(rep("-", 32), collapse = ""), " Cluster metric generation ", 
                     paste0(rep("-", 32), collapse = ""), "||\nStarted process at: ", Sys.time()))

cat(paste0("\nIf at any point the process cuts off abruptly with no success message, ", 
           "please see the log file.\nPreparing data for processing.\n"))

# FOR USER: replace the filename variables with quoted file paths if you don't want to input them each time
# the cluster assignments, in form: || isolates | height 0 | height 1 | ... ||

outputDetails(paste0("\nPART 1 OF 3: Data processing ", paste0(rep(".", 66), collapse = ""), "\n", 
                     "  Preparing data for processing\n  Loading and formatting datafiles..."), newcat = TRUE)

# the coded datasets (where the isolates are now replaced by positive integers)
time1_raw <- readBaseData(arg$tp1, 1)
time2_raw <- readBaseData(arg$tp2, 2)

all_isolates <- c(time1_raw$isolate, time2_raw$isolate) %>% unique() %>% 
  as_tibble() %>% set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

t1_coded <- time1_raw %>% codeIsolates(., "tp1", all_isolates)
t2_coded <- time2_raw %>% codeIsolates(., "tp2", all_isolates)

t2_colnames <- t2_coded$tp2_h %>% unique() %>% sort()
message("  Successfully read in datafiles")

# Note: assumption made that first column is labelled "isolate"
outputDetails("  Processing time point 1 (TP1) clusters for easier data handling", newcat = TRUE)
t1_comps <- t1_coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>% compsSet(., "TP1", indicate_progress = TRUE)

outputDetails("  Collecting and counting novel isolates in each TP2 cluster...", newcat = TRUE)
novels <- setdiff(t2_coded$isolate, t1_coded$isolate)
counting_novels <- t2_coded %>% filter(isolate %in% novels) %>%
  group_by(id) %>% summarise(num_novs = n(), .groups = "drop")

outputDetails("  Processing time point 2 (TP2) clusters for easier data handling...", newcat = TRUE)
t2_comps <- t2_coded %>% rename(tp_h = tp2_h, tp_cl = tp2_cl) %>% 
  compsSet(., "TP2", indicate_progress = TRUE) %>% left_join(., counting_novels, by = "id")
t2_comps$num_novs[is.na(t2_comps$num_novs)] <- 0

### BASE CASE:
outputDetails(paste0("\nPART 2 OF 3: Tracking and flagging clusters for base case ", 
                     paste0(rep(".", 41), collapse = "")), newcat = TRUE)

# this should be '1', the first column is the isolates
outputDetails("  Tracking clusters", newcat = TRUE)
heights <- strsplit(arg$heights, split = ",") %>% unlist()
h_before <- heights[1]     # h_before <- unique(t1_coded$tp1_h)[1]

outputDetails(paste0("  Collecting height data for base case, height ", h_before, "..."), newcat = TRUE)
hdata <- t1_comps %>% filter(tp1_h == h_before) %>% arrange(tp1_h, tp1_cl)

cc <- trackClusters(hdata, t2_comps, t2_colnames, t1_coded, t2_coded, indicate_progress = TRUE) %>%
  createID(., "tp1", "tp1_h", "tp1_cl")

# note: the way these are flagged, we can see how many TP2 clusters were formed from the
# selected TP1 cluster before we get to the "key growth" where a cluster went from its state
# at TP1 to its state at TP2 with a remarkable amount of growth

outputDetails("  Flagging clusters", newcat = TRUE)
hbdata <- cc %>% add_column(flag = cc$id)

formatC(h_before, width = nchar(max(t1_coded$tp1_h)), format = "d", flag = "0") %>% 
  saveData(dtype=3, tmp=hbdata, h=.)

bef_comps <- t1_comps %>% filter(tp1_h == h_before) %>%
  set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))

outputDetails(paste0("\nPART 3 OF 3: Tracking and flagging clusters for the rest of the heights (", 
                     length(heights) - 1, " of them) ........."), newcat = TRUE)

outputDetails(paste0("  This may take some time. \n  Note that for a more detailed ", 
                     "look at progress, you can keep an eye on the outputs \n  ", 
                     "directory, where updated data will be saved after each height.\n"))

# ALL OTHER HEIGHTS
outputDetails("  Collecting data for other heights: ", newcat = TRUE)

fcb <- txtProgressBar(min = 0, max = length(heights[-1]), initial = 0, style = 3)
for (j in 1:length(heights[-1])) {
  setTxtProgressBar(fcb, j)
  
  h_after <- heights[-1][j]
  message(paste0("  Height ", j, " / ", length(heights)))
  
  aft_comps <- t1_comps %>% filter(tp1_h == h_after) %>%
    set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft"))

  # Part 1: identifying clusters that have not changed from the previous height
  ss <- noChange(aft_comps, bef_comps, hbdata)
  
  # Part 2: tracking clusters that changed
  hdata <- aft_comps %>% filter(!(id_aft %in% ss$id)) %>% set_colnames(colnames(t1_comps))
  
  # indicate_progress = FALSE
  cc <- trackClusters(hdata, t2_comps, t2_colnames, t1_coded, t2_coded, indicate_progress = FALSE) %>% 
    createID(., "tp1", "tp1_h", "tp1_cl")
  
  # Part 3: flagging clusters to indicate when they were first seen
  if (nrow(cc) > 0) {
    hnew <- cc %>% add_column(flag = cc$id) %>% bind_rows(., ss) %>% arrange(tp1_h, tp1_cl)
  }else {
    hnew <- ss %>% arrange(tp1_h, tp1_cl)
  }
  formatC(h_after, width = nchar(max(t1_coded$tp1_h)), format = "d", flag = "0") %>% 
    saveData(dtype=3, tmp=hnew, h=.)

  # Part 4 - prepping data for next iteration
  hbdata <- hnew; h_before <- h_after
  bef_comps <- aft_comps %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
}
close(fcb)

stopwatch[2] <- Sys.time(); saveData(dtype = 7, sw = stopwatch)
outputDetails("  Merging height data into a single results file - note, just tracked clusters.", newcat = TRUE)
timeTaken(pt = "data collection", stopwatch) %>% outputDetails(., newcat = TRUE)
outputDetails(paste0("  Successfully collected data for all heights."), newcat = TRUE)

op <- "outputs/summary/"
data_type <- "manual_heights"
if (data_type == "nawc") {
  
  a2 <- read.table(paste0(op, "stats.tsv"), sep = "\t", header = TRUE) %>% as_tibble()
  a2$Neighbour.AWC[is.na(a2$Neighbour.AWC)] <- 0
  
  # differences x[i+1] - x[i]
  y4 <- which(diff(sign(diff(a2$Neighbour.AWC))) > 0) + 1
  # find the points that shrank from the left and grew to the right, note the indexing starts from 0
  
  pt <- a2[y4,] %>% filter(Neighbour.AWC < 0.95)
  
  if (length(colnames(a2)) < 4) {
    a2$Selected.for.Analysis <- 0
    a2$Selected.for.Analysis[a2$Threshold %in% pt$Threshold] <- 1
    write.table(a2, file=paste0(op, "stats_with_analysis.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  }
  p <- ggplot(a2, aes(x = Threshold, y = Neighbour.AWC)) + geom_step() + 
    geom_point(data = pt, aes(x = Threshold, y = Neighbour.AWC), color = "red") +
    labs(x = "goeBURST Threshold", y = "")
  
  ggsave(paste0(op,"nawc_stability_plot.png"), device="png", plot=p, width=16, height=9, units="in", dpi=320)
  print(pt$Threshold)
  
  Sys.time()
  c4 <- lapply(pt$Threshold, function(hval) {
    oneh <- readRDS(paste0("outputs/height_data/h", hval, ".Rds")) %>% 
      left_join(., dfx, by = c("flag" = "first_flag")) %>% 
      arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)
    
    oneHeightV2(hval, novels, t2_comps, t1_comps, oneh) %>% return()
  }) %>% bind_rows()
  
  Sys.time()
}else {
  
  
  ba <- time1_raw %>% melt(id = "isolate") %>% as_tibble() %>% 
    rename("tp1_cl" = value, "tp1_h" = variable) %>% 
    mutate(across("tp1_h", as.character)) %>% 
    createID(., "tp1", "tp1_h", "tp1_cl") %>% rename("tp1_id" = id)

  b2 <- time2_raw %>% melt(id = "isolate") %>% as_tibble() %>% 
    rename(tp2_cl = value, tp2_h = variable) %>% 
    mutate(across(tp2_h, as.character)) %>% 
    createID(., "tp2", "tp2_h", "tp2_cl") %>% rename(tp2_id = id)

  # Identifying the last time each cluster was seen
  allcomps <- t1_comps$composition %>% unique()
  dfx <- lapply(1:length(allcomps), function(i) {
    t1_comps %>% filter(composition == allcomps[i]) %>% slice(1, nrow(.)) %>% 
      pull(id) %>% t() %>% data.frame(stringsAsFactors = FALSE)
    }) %>% bind_rows() %>% as_tibble() %>% set_colnames(c("first_flag", "last_flag"))

  c4 <- lapply(heights, function(hval) {
    oneh <- readRDS(paste0("outputs/height_data/h", hval, ".Rds")) %>% 
      left_join(., dfx, by = c("flag" = "first_flag")) %>% 
      arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)
    
    oneHeightV2(hval, novels, t2_comps, t1_comps, oneh, b1, b2) %>% return()
  }) %>% bind_rows()
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

write.csv(clusters_formatted, paste0(op,"TP1_cluster_results.csv"), row.names = FALSE)

isolates_formatted <- time1_raw %>% select(isolate, as.character(pt$Threshold)) %>%
  melt(., id = "isolate") %>% as_tibble() %>% rename(tp1_h = variable, tp1_cl = value) %>%
  mutate(across(tp1_h, as.character)) %>% mutate(across(tp1_h, as.integer)) %>%
  leadingZeros(., "tp1_h", "h") %>% leadingZeros(., "tp1_cl", "c") %>%
  right_join(., c4, by = c("tp1_h", "tp1_cl")) %>% arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>%
  select(isolate, colnames(c4)) %>% set_colnames(c("Isolates", colnames(clusters_formatted)))

write.csv(isolates_formatted, paste0(op, "TP1_strain_results.csv"), row.names = FALSE)

outputDetails("  Closing all connections...", newcat = TRUE)
closeAllConnections()

