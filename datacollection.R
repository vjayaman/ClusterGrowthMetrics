#! /usr/bin/env Rscript
msg <- file("outputs/logfile_datacollection.txt", open="wt")
sink(msg, type="message")

suppressWarnings(suppressPackageStartupMessages(source("functions/tracking_functions.R")))
source("class_definitions.R")

option_list <- list(
  make_option(c("-a", "--tp1"), metavar = "file", default = NULL, help = "Time point 1 file name (TP1)"),
  make_option(c("-b", "--tp2"), metavar = "file", default = NULL, help = "Time point 2 file name (TP2)"),
  make_option(c("-x", "--heights"), metavar = "character", default = NULL,
              help = paste0("A string of comma-delimited numbers, e.g. '50,75,100' to ", 
                            "use as heights for metric generation (strain table outputs)")))

arg <- parse_args(OptionParser(option_list=option_list))
# arg <- tibble(tp1 = "t1_clusters_processed.csv", tp2 = "t2_clusters_processed.csv", heights = "0,5,25")
outputDetails(paste0("\n||", paste0(rep("-", 32), collapse = ""), " Cluster metric generation ", 
                     paste0(rep("-", 32), collapse = ""), "||\nStarted process at: ", Sys.time()))
cat(paste0("\nIf at any point the process cuts off with no success message, please see the log file.\n"))
outputDetails(paste0("\nPART 1 OF 3: Data processing ", paste0(rep(".", 66), collapse = "")), newcat = TRUE)

stopwatch <- list("start_time" = as.character.POSIXt(Sys.time()), "end_time" = NULL)

# TP DATA PREPARATION ------------------------------------------------------------------------------------------
f1 <- readBaseData(arg$tp1, 1, reader::get.delim(arg$tp1))
f2 <- readBaseData(arg$tp2, 2, reader::get.delim(arg$tp2))
colnames(f1)[1] <- colnames(f2)[1] <- "isolate"

heights <- strsplit(arg$heights, split = ",") %>% unlist()

all_isolates <- unique(c(f1$isolate, f2$isolate)) %>% as_tibble() %>% 
  set_colnames("char_isolate") %>% rowid_to_column("num_isolate")

ph <- max(nchar(colnames(f1)[-1]), nchar(colnames(f2)[-1]))
pc <- f1 %>% select(-isolate) %>% max(., f1 %>% select(-isolate)) %>% nchar()

outputDetails("  Processing timepoint clusters for easier data handling and flagging clusters", newcat = TRUE)
tp1 <- Timedata$new("tp1", raw = f1, all_isolates, pad_height = ph, pad_cluster = pc)$set_comps()$flag_clusters()

tp2 <- Timedata$new("tp2", raw = f2, all_isolates, pad_height = ph, pad_cluster = pc)$set_comps()$set_cnames()

outputDetails("  Collecting and counting novel isolates in TP2 clusters ...", newcat = TRUE)
novels <- setdiff(tp2$coded$isolate, tp1$coded$isolate)

tp2$comps <- tp2$coded %>% filter(isolate %in% novels) %>% 
  group_by(tp2_id) %>% 
  summarise(num_novs = n(), .groups = "drop") %>% 
  left_join(tp2$comps, ., by = "tp2_id") %>% 
  mutate(num_novs = ifelse(is.na(num_novs), 0, num_novs))

tp2$flag_clusters()$coded_status(novels)
tp1$coded_status(novels)

# BASE CASE (FIRST HEIGHT) -------------------------------------------------------------------------------------
outputDetails(paste0("\nPART 2 OF 3: Tracking and flagging clusters for base case ", 
                     paste0(rep(".", 37), collapse = "")), newcat = TRUE)
outputDetails(paste0("  Collecting height data for base case, height ", heights[1], "..."), newcat = TRUE)

hx <- Heightdata$new(starter = heights[1], t1_comps = tp1$comps, hvals = heights)$
  clust_tracking(tp2$comps, tp2$cnames, tp1$coded, tp2$coded, TRUE)$
  update_iteration()

# REST OF THE HEIGHTS ------------------------------------------------------------------------------------------
if (length(heights) > 1) {
  outputDetails(paste0("\nPART 3 OF 3: Tracking and flagging clusters for the rest of the heights (",
                       length(heights) - 1, " of them) ..........."), newcat = TRUE)
  outputDetails(paste0("  This may take some time. \n  For a more detailed look at progress, ", 
                       "see the logfile in the outputs directory.\n"))
  outputDetails("  Collecting data for other heights: ", newcat = TRUE)

  fcb <- txtProgressBar(min = 0, max = length(heights[-1])*2, initial = 0, style = 3)
  for (j in 1:length(heights[-1])) {
    hx$h_after <- heights[-1][j]
    message(paste0("  Height ", j + 1, " / ", length(heights)))

    # Part 1: unchanged(): identifying clusters that have not changed from the previous height
    # Part 2: takes comps for new height, previous height, and the tracked data for the previous height
    #   --> identifies clusters that have not changed since the previous height, and reuses their tracking data
    hx$post_data(tp1$comps)$unchanged()

    # Part 2: tracking clusters that changed, saving to results list, and prepping variable for next height
    hx$comps <- hx$aft %>% filter(!(id_aft %in% hx$same$tp1_id)) %>% set_colnames(colnames(tp1$comps))
    
    setTxtProgressBar(fcb, j*2 - 1)
    
    hx$clust_tracking(tp2$comps, tp2$cnames, tp1$coded, tp2$coded, FALSE)$
      update_iteration()$reset_values()
    
    setTxtProgressBar(fcb, j*2)
  }
  close(fcb)
}else {
  outputDetails(paste0("\nPART 3 OF 3: Only one threshold provided, so no further tracking necessary"), newcat = TRUE)
}

outputDetails("  Saving the data in two separate files, with cluster and strain identifiers.\n", newcat = TRUE)

clusters_just_tp1 <- lapply(heights, function(h) {
  hx$results[[h]] %>% left_join(., tp1$flagged) %>% left_join(., tp2$flagged) %>% 
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
    oneHeight(novels, tp1$status, tp2$status, .) %>% return()
}) %>% bind_rows()

isolates_base <- tp1$melted %>% mutate(across(tp1_h, as.integer)) %>% 
  filter(tp1_h %in% heights) %>% 
  left_join(., clusters_just_tp1, by = c("tp1_id", "tp1_h", "tp1_cl")) %>% 
  arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)

# NOVELS -------------------------------------------------------------------------------------------------------
first_nov_flag <- tp2$status %>% filter(!is.na(status)) %>% group_by(isolate) %>% slice(1) %>% ungroup() %>% pull(tp2_id)

novels_only_tracking <- tp2$flagged %>% filter(tp2_id %in% first_nov_flag) %>% 
  left_join(., tp2$comps[,c("tp2_id", "num_novs")]) %>% arrange(tp2_h, tp2_cl) %>% 
  add_column(tp1_id = NA, tp1_h = NA, tp1_cl = NA, first_tp1_flag = NA, last_tp1_flag = NA)

novel_asmts <- tp2$melted %>% mutate(across(tp2_h, as.integer)) %>% 
  filter(isolate %in% setdiff(tp2$raw$isolate, tp1$raw$isolate))

# Pure novel TP2 cluster
pure_novels <- novels_only_tracking %>% filter(num_novs == tp2_cl_size) %>% 
  mutate(tp1_cl_size = 0, add_TP1 = 0, 
         actual_size_change = tp2_cl_size - tp1_cl_size, 
         actual_growth_rate = ((tp2_cl_size - tp1_cl_size) / tp1_cl_size) %>% round(digits = 3), 
         new_growth = (tp2_cl_size / (tp2_cl_size - num_novs)) %>% round(digits = 3)) %>% 
  right_join(novel_asmts, ., by = c("tp2_h", "tp2_cl", "tp2_id")) %>% select(colnames(isolates_base))

# Mixed novel TP2 clusters
mixed_novels <- novels_only_tracking %>% filter(num_novs != tp2_cl_size) %>% 
  left_join(., novel_asmts, by = c("tp2_id", "tp2_h", "tp2_cl"))

# Mixed novels clusters should inherit data from the tracked TP1 strains
all_mixed <- isolates_base %>% 
  select(-isolate, -tp1_id, -tp1_h, -tp1_cl, -first_tp1_flag, -last_tp1_flag) %>% 
  left_join(mixed_novels, .) %>% unique() %>% select(colnames(isolates_base))

# FULL STRAIN FILE ---------------------------------------------------------------------------------------------
# note: two types of novel clusters, those that are fully novel, and those that are not
isolates_file <- bind_rows(isolates_base, pure_novels) %>% bind_rows(., all_mixed) %>% 
  mutate(novel = ifelse(isolate %in% setdiff(tp2$raw$isolate, tp1$raw$isolate), 1, 0)) %>% 
  rename(strain = isolate)
isolates_file[,c("tp1_h", "tp2_h")] %<>% apply(., 2, padCol, padval = ph, padchr = "h")
isolates_file[,c("tp1_cl", "tp2_cl")] %<>% apply(., 2, padCol, padval = pc, padchr = "c")
  
isolates_file %>% 
  select(strain, novel, first_tp2_flag, tp2_h, tp2_cl, tp2_cl_size, last_tp2_flag, tp1_id, tp1_h, tp1_cl, tp1_cl_size, 
         first_tp1_flag, last_tp1_flag, add_TP1, num_novs, actual_size_change, actual_growth_rate, new_growth) %>% 
  set_colnames(c("Strain", "Novel", "First time this cluster was seen in TP2", "TP2 height", "TP2 cluster", 
                 "TP2 cluster size", "Last time this cluster was seen in TP2", "TP1 ID", "TP1 height", "TP1 cluster",
                 "TP1 cluster size", "First time this cluster was seen in TP1", 
                 "Last time this cluster was seen in TP1", "Number of additional TP1 strains in the TP2 match", 
                 "Number of novels in the TP2 match", "Actual cluster size change (TP2 size - TP1 size)", 
                 "Actual growth rate = (TP2 size - TP1 size) / (TP1 size)", 
                 "Novel growth = (TP2 size) / (TP2 size - number of novels)")) %>% 
  write.table(., file.path("outputs/","CGM_strain_results.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

# WRAPPING THINGS UP -------------------------------------------------------------------------------------------
stopwatch[["end_time"]] <- as.character.POSIXt(Sys.time())

outputDetails(paste0("\nSuccessfully collected data for all heights."), newcat = TRUE)
timeTaken(pt = "data collection", stopwatch) %>% outputDetails(., newcat = TRUE)
outputDetails(paste0("\n||", paste0(rep("-", 28), collapse = ""), " End of cluster metric generation ", 
                     paste0(rep("-", 29), collapse = ""), "||\n"))
closeAllConnections()
