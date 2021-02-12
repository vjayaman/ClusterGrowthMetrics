#! /usr/bin/env Rscript
msg <- file("outputs/logfile_datacollection.txt", open="wt")
sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("functions/tracking_functions.R")))

setClass("timedata", slots = list(name = "character", raw = "data.frame", coded = "data.frame", 
                                  comps = "data.frame", melted = "data.frame"))

setClass("heightdata", slots = list(h_before = "character", h_after = "character", 
                                    comps = "data.frame", changed = "data.frame", same = "data.frame", 
                                    tracked = "data.frame", bef = "data.frame", aft = "data.frame"))
option_list <- list(
  
  make_option(c("-a", "--tp1"), metavar = "file", default = NULL, help = "Time point 1 file name"),
  
  make_option(c("-b", "--tp2"), metavar = "file", default = NULL, help = "Time point 2 file name"),
  
  make_option(c("-x", "--heights"), metavar = "character", default = NULL,
              help = paste0("A string of comma-delimited numbers, e.g. '50,75,100' to ", 
                            "use as heights for which to generate cluster and strain tables")))

arg <- parse_args(OptionParser(option_list=option_list))

stopwatch <- rep(0,2) %>% set_names(c("start_time", "end_time"))
stopwatch[1] <- Sys.time()

outputDetails(paste0("\n||", paste0(rep("-", 32), collapse = ""), " Cluster metric generation ", 
                     paste0(rep("-", 32), collapse = ""), "||\nStarted process at: ", Sys.time()))
cat(paste0("\nIf at any point the process cuts off abruptly with no success message, ", 
           "please see the log file.\nPreparing data for processing.\n"))
outputDetails(paste0("\nPART 1 OF 3: Data processing ", paste0(rep(".", 66), collapse = ""), "\n", 
                     "  Preparing data for processing\n  Loading and formatting datafiles..."), 
              newcat = TRUE)

# DATA PREPARATION
tpt1 <- new("timedata", name = "tp1", raw = readBaseData(arg$tp1, 1))
tpt2 <- new("timedata", name = "tp2", raw = readBaseData(arg$tp2, 2))

message("  Successfully read in datafiles")
all_isolates <- c(tpt1@raw$isolate, tpt2@raw$isolate) %>% unique() %>% 
  as_tibble() %>% set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

tpt1@coded <- tpt1@raw %>% codeIsolates(., tpt1@name, all_isolates)
tpt1@melted <- meltedIDs(tpt1@raw, "tp1")

tpt2@coded <- tpt2@raw %>% codeIsolates(., tpt2@name, all_isolates)
tpt2@melted <- meltedIDs(tpt2@raw, "tp2")

t2_colnames <- tpt2@coded$tp2_h %>% unique() %>% sort()

# Note: assumption made that first column is labelled "isolate"
outputDetails("  Processing time point 1 (TP1) clusters for easier data handling", newcat = TRUE)

tpt1@comps <- tpt1@coded %>% rename(tp_h = tp1_h, tp_cl = tp1_cl) %>% 
  compsSet(., "TP1", indicate_progress = TRUE)

outputDetails("  Collecting and counting novel isolates in each TP2 cluster...", newcat = TRUE)
novels <- setdiff(tpt2@coded$isolate, tpt1@coded$isolate)
counting_novels <- tpt2@coded %>% 
  filter(isolate %in% novels) %>% group_by(id) %>% 
  summarise(num_novs = n(), .groups = "drop")

outputDetails("  Processing time point 2 (TP2) clusters for easier data handling...", newcat = TRUE)
tpt2@comps <- tpt2@coded %>% rename(tp_h = tp2_h, tp_cl = tp2_cl) %>% 
  compsSet(., "TP2", indicate_progress = TRUE) %>% left_join(., counting_novels, by = "id")
tpt2@comps$num_novs[is.na(tpt2@comps$num_novs)] <- 0

### BASE CASE:
outputDetails(paste0("\nPART 2 OF 3: Tracking and flagging clusters for base case ", 
                     paste0(rep(".", 41), collapse = "")), newcat = TRUE)

# this should be '0', the first column is the isolates
outputDetails("  Tracking clusters", newcat = TRUE)
heights <- strsplit(arg$heights, split = ",") %>% unlist() # "0,5,25"

outputDetails(paste0("  Collecting height data for base case, height ", heights[1], "..."), newcat = TRUE)

hx <- new("heightdata", h_before = heights[1])
hx@comps <- tpt1@comps %>% filter(tp1_h == hx@h_before) %>% arrange(tp1_h, tp1_cl)
hx@changed <- trackClusters(hx@comps, tpt2@comps, t2_colnames, 
                            tpt1@coded, tpt2@coded, indicate_progress = TRUE) %>%
  createID(., "tp1", "tp1_h", "tp1_cl")

outputDetails("  Flagging clusters", newcat = TRUE)
hx@tracked <- hx@changed %>% add_column(flag = hx@changed$id)

formatC(as.integer(hx@h_before), width=nchar(max(tpt1@coded$tp1_h)), 
        format="d", flag="0") %>% saveData(tmp=hx@tracked, h=.)

hx@bef <- tpt1@comps %>% filter(tp1_h == hx@h_before) %>%
  set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))

outputDetails(paste0("\nPART 3 OF 3: Tracking and flagging clusters for the rest of the heights (", 
                     length(heights) - 1, " of them) ........."), newcat = TRUE)
outputDetails(paste0("  This may take some time. \n  Note that for a more detailed ", 
                     "look at progress, you can keep an eye on the outputs \n  ", 
                     "directory, where updated data will be saved after each height.\n"))
outputDetails("  Collecting data for other heights: ", newcat = TRUE)

fcb <- txtProgressBar(min = 0, max = length(heights[-1]), initial = 0, style = 3)
for (j in 1:length(heights[-1])) {
  setTxtProgressBar(fcb, j)
  
  hx@h_after <- heights[-1][j]
  message(paste0("  Height ", j + 1, " / ", length(heights)))
  
  hx@aft <- tpt1@comps %>% filter(tp1_h == hx@h_after) %>%
    set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft"))

  # Part 1: identifying clusters that have not changed from the previous height
  hx@same <- noChange(hx@aft, hx@bef, hx@tracked)
  
  # Part 2: tracking clusters that changed
  hx@comps <- hx@aft %>% filter(!(id_aft %in% hx@same$id)) %>% set_colnames(colnames(tpt1@comps))
  
  # indicate_progress = FALSE
  hx@changed <- trackClusters(hx@comps, tpt2@comps, t2_colnames, 
                              tpt1@coded, tpt2@coded, indicate_progress = FALSE) %>% 
    createID(., "tp1", "tp1_h", "tp1_cl")
  
  # Part 3: flagging clusters to indicate when they were first seen
  if (nrow(hx@changed) > 0) {
    hx@tracked <- hx@changed %>% add_column(flag = hx@changed$id) %>% 
      bind_rows(., hx@same) %>% arrange(tp1_h, tp1_cl)
  }else {
    hx@tracked <- hx@same %>% arrange(tp1_h, tp1_cl)
  }
  formatC(as.integer(hx@h_after), width = nchar(max(tpt1@coded$tp1_h)), format = "d", flag = "0") %>% 
    saveData(tmp = hx@tracked, h = .)

  # Part 4 - prepping data for next iteration
  hx@h_before <- hx@h_after
  hx@bef <- hx@aft %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
}
close(fcb)

outputDetails("  Merging height data into a single results file.", newcat = TRUE)
op <- "outputs/summary/"

# Identifying the last time each cluster was seen
a1 <- tpt1@comps %>% group_by(composition) %>% slice(1, n()) %>% 
  add_column(type = rep(c("first", "last"), nrow(.)/2)) %>% dplyr::ungroup()

a2a <- a1[which(a1$type=="first"), c("id", "composition")] %>% rename("first_flag" = "id")
a2 <- a1[which(a1$type=="last"), ] %>% 
  select(id, composition) %>% rename("last_flag" = "id") %>% 
  full_join(a2a, ., by = "composition") %>% 
  select(composition, first_flag, last_flag) %>% 
  left_join(tpt1@comps, ., by = "composition") %>% 
  select(c("tp1_h", "tp1_cl", "id", "first_flag", "last_flag"))

hfiles <- lapply(heights, function(h) {
  as.integer(h) %>% 
    formatC(., width = nchar(max(tpt1@coded$tp1_h)), format = "d", flag = "0") %>% 
    paste0("h", ., ".Rds") %>% tibble(h, f = .)
}) %>% bind_rows()

outputDetails("  Saving the data in two separate files, with cluster and strain identifiers.", newcat = TRUE)
datafiles <- lapply(1:nrow(hfiles), function(i) {
  readRDS(paste0("outputs/height_data/", hfiles$f[i])) %>% 
    left_join(., a2, by = c("tp1_h", "tp1_cl", "id")) %>% 
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
    oneHeight(hfiles$h[i], novels, tpt1@comps, tpt2@comps, ., tpt1@melted, tpt2@melted) %>% return()
}) %>% bind_rows()

datafiles$actual_growth_rate %<>% format(., digits = 3, nsmall = 3)
datafiles$new_growth %<>% format(., digits = 3, nsmall = 3)

resultFiles(datafiles, op, heights, tpt1@raw, tpt1@melted)

stopwatch[2] <- Sys.time()
timeTaken(pt = "data collection", stopwatch) %>% outputDetails(., newcat = TRUE)
outputDetails(paste0("  Successfully collected data for all heights."), newcat = TRUE)
outputDetails("  Closing all connections...", newcat = TRUE)
closeAllConnections()
