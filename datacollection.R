#! /usr/bin/env Rscript
msg <- file("outputs/logfile_datacollection.txt", open="wt")
sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("functions/tracking_functions.R")))
source("class_definitions.R")

option_list <- list(
  make_option(c("-a", "--tp1"), metavar = "file", default = NULL, help = "Time point 1 file name"),
  make_option(c("-b", "--tp2"), metavar = "file", default = NULL, help = "Time point 2 file name"),
  make_option(c("-d", "--delimiter"), metavar = "character", default = "\t"), 
  make_option(c("-x", "--heights"), metavar = "character", default = NULL,
              help = paste0("A string of comma-delimited numbers, e.g. '50,75,100' to ", 
                            "use as heights for which to generate cluster and strain tables")))

arg <- parse_args(OptionParser(option_list=option_list))

stopwatch <- list(as.character.POSIXt(Sys.time()))

outputDetails(paste0("\n||", paste0(rep("-", 32), collapse = ""), " Cluster metric generation ", 
                     paste0(rep("-", 32), collapse = ""), "||\nStarted process at: ", Sys.time()))
cat(paste0("\nIf at any point the process cuts off abruptly with no success message, ", 
           "please see the log file.\nPreparing data for processing.\n"))
outputDetails(paste0("\nPART 1 OF 3: Data processing ", paste0(rep(".", 66), collapse = ""), "\n", 
                     "  Preparing data for processing\n  Loading and formatting datafiles..."), 
              newcat = TRUE)

# DATA PREPARATION
f1 <- readBaseData(arg$tp1, 1, "\t")#arg$delimiter)
f2 <- readBaseData(arg$tp2, 2, "\t")#arg$delimiter)
heights <- strsplit(arg$heights, split = ",") %>% unlist()

# f1 <- readBaseData("t1_clusters_processed.csv", 1, "\t")
# f2 <- readBaseData("t2_clusters_processed.csv", 2, "\t")
# heights <- "100"

colnames(f1)[1] <- colnames(f2)[1] <- "isolate"

all_isolates <- unique(c(f1$isolate, f2$isolate)) %>% as_tibble() %>% 
  set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

tp1 <- Timedata$new("tp1", raw = f1, all_isolates)
tp2 <- Timedata$new("tp2", raw = f2, all_isolates)

# Note: assumption made that first column is labelled "isolate"
outputDetails("  Processing time point 1 (TP1) clusters for easier data handling", newcat = TRUE)
tp1$coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>% tp1$set_comps(.)

t2_colnames <- tp2$coded$tp2_h %>% unique() %>% sort()

outputDetails("  Collecting and counting novel isolates in each TP2 cluster...", newcat = TRUE)
novels <- setdiff(tp2$coded$isolate, tp1$coded$isolate)
counting_novels <- tp2$coded %>% filter(isolate %in% novels) %>%
  group_by(id) %>% summarise(num_novs = n(), .groups = "drop")

outputDetails("  Processing time point 2 (TP2) clusters for easier data handling...", newcat = TRUE)
tp2$coded %>% rename(tp_h = tp2_h, tp_cl = tp2_cl) %>% tp2$set_comps()
tp2$comps <- tp2$comps %>% left_join(., counting_novels, by = "id")
tp2$comps$num_novs[is.na(tp2$comps$num_novs)] <- 0

### BASE CASE ---------------------------------------------------------------------------------
outputDetails(paste0("\nPART 2 OF 3: Tracking and flagging clusters for base case ", 
                     paste0(rep(".", 37), collapse = "")), newcat = TRUE)
outputDetails(paste0("  Collecting height data for base case, height ", heights[1], "..."), newcat = TRUE)

hx <- Heightdata$new(h_before = heights[1], t1_comps = tp1$comps)

outputDetails("  Tracking and flagging clusters", newcat = TRUE)
hx$clust_tracking(tp2$comps, t2_colnames, tp1$coded, tp2$coded, TRUE)$update_tracking()$prior_data(tp1$comps)

hx$saveTempFile(tp1$coded, hx$h_before, "outputs")

if (length(heights) > 1) {
  outputDetails(paste0("\nPART 3 OF 3: Tracking and flagging clusters for the rest of the heights (", 
                       length(heights) - 1, " of them) ..........."), newcat = TRUE)
  outputDetails(paste0("  This may take some time. \n  Note that for a more detailed ", 
                       "look at progress, you can keep an eye on the outputs \n  ", 
                       "directory, where updated data will be saved after each height.\n"))
  outputDetails("  Collecting data for other heights: ", newcat = TRUE)

  fcb <- txtProgressBar(min = 0, max = length(heights[-1]), initial = 0, style = 3)
  for (j in 1:length(heights[-1])) {
    setTxtProgressBar(fcb, j)
    
    hx$h_after <- heights[-1][j]
    message(paste0("  Height ", j + 1, " / ", length(heights)))
    # unchanged(): identifying clusters that have not changed from the previous height
    hx$post_data(tp1$comps)$unchanged()
    # Part 2: tracking clusters that changed
    hx$comps <- hx$aft %>% filter(!(id_aft %in% hx$same$id)) %>% set_colnames(colnames(tp1$comps))
    hx$clust_tracking(tp2$comps, t2_colnames, tp1$coded, tp2$coded, FALSE)
    
    # Part 3: flagging clusters to indicate when they were first seen
    if (nrow(hx$changed) > 0) {
      # hx$add_flag()
      hx$tracked <- bind_rows(hx$tracked, hx$same) %>% arrange(tp1_h, tp1_cl)
    }else {
      hx$tracked <- hx$same %>% arrange(tp1_h, tp1_cl)
    }
    hx$saveTempFile(tp1$coded, hx$h_after, "outputs")
    hx$update_iteration()
  }
  close(fcb)
}else {
  outputDetails(paste0("\nPART 3 OF 3: Only one threshold provided, so no further tracking necessary"), newcat = TRUE)
}

outputDetails("  Merging height data into a single results file.", newcat = TRUE)  

# Identifying the first and last time each TP1 cluster was seen in TP1
t1_fal <- tp1$comps %>% group_by(composition) %>% slice(1, n()) %>% 
  add_column(type = rep(c("first", "last"), nrow(.)/2)) %>% dplyr::ungroup()
# 
t1_ff <- t1_fal %>% filter(type == "first") %>% select(id, composition) %>% rename("first_tp1_flag" = "id")
t1_lf <- t1_fal %>% filter(type == "last")  %>% select(id, composition) %>% rename("last_tp1_flag" = "id")
t1_ff_lf <- full_join(t1_ff, t1_lf, by = "composition") %>% 
  left_join(tp1$comps, ., by = "composition") %>% 
  select(-composition) # tp1_h, tp1_cl, id, tp1_cl_size, first_flag, last_flag

# Identifying the first and last time each TP2 cluster was seen in TP2
t2_fal <- tp2$comps %>% group_by(composition) %>% slice(1, n()) %>% 
  add_column(type = rep(c("first", "last"), nrow(.)/2)) %>% dplyr::ungroup()
t2_ff <- t2_fal %>% filter(type == "first") %>% select(id, composition) %>% rename("first_tp2_flag" = "id")
t2_lf <- t2_fal %>% filter(type == "last")  %>% select(id, composition) %>% rename("last_tp2_flag" = "id")
t2_ff_lf <- full_join(t2_ff, t2_lf, by = "composition") %>% 
  left_join(tp2$comps, ., by = "composition") %>% 
  rename(tp2_id = id) %>% 
  select(tp2_id, tp2_h, tp2_cl, tp2_cl_size, first_tp2_flag, last_tp2_flag)

hfiles <- lapply(heights, function(h) {
  formatC(as.integer(h), width = nchar(max(tp1$coded$tp1_h)), format = "d", flag = "0") %>%
    paste0("h", ., ".Rds") %>% tibble(h, f = .)
}) %>% bind_rows()

outputDetails("  Saving the data in two separate files, with cluster and strain identifiers.\n", newcat = TRUE)
datafiles <- lapply(1:nrow(hfiles), function(i) {
  tmp <- readRDS(file.path("outputs", hfiles$f[i])) %>%  
    left_join(., t1_ff_lf, by = c("tp1_h", "tp1_cl", "id", "tp1_cl_size")) %>% rename(tp1_id = id) %>% 
    left_join(., t2_ff_lf, by = c("tp2_h", "tp2_cl", "tp2_cl_size")) %>% 
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)
    
  oneHeight(hfiles$h[i], novels, tp1$comps, tp2$comps, tmp, tp1$melted, tp2$melted) %>% return()
}) %>% bind_rows()

resultFiles(datafiles, "outputs", heights, tp1$raw, tp2$raw, tp1$melted)

stopwatch <- append(stopwatch, values = as.character.POSIXt(Sys.time())) %>% 
  set_names(c("start_time", "end_time"))

outputDetails(paste0("  Successfully collected data for all heights."), newcat = TRUE)

timeTaken(pt = "data collection", stopwatch) %>% outputDetails(., newcat = TRUE)
outputDetails(paste0("\n||", paste0(rep("-", 28), collapse = ""), " End of cluster metric generation ", 
                     paste0(rep("-", 29), collapse = ""), "||"))

# outputDetails("  Closing all connections...", newcat = TRUE)
closeAllConnections()








