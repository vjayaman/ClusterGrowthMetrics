
# assertthat::see_if(identical(a1, a2), msg = paste0("File ", j, " is not the same in both directories."))
timeTaken <- function(pt, sw) {
  t1 <- (sw[['end_time']] - sw[['start_time']])/60/60
  t2 <- abs(t1 - trunc(t1))*60
  t3 <- abs(t2 - trunc(t2))*60
  
  word1 <- if_else(trunc(t1) == 1, "hour", "hours")
  word2 <- if_else(trunc(t2) == 1, "minute", "minutes")
  word3 <- if_else(trunc(t3) == 1, "second", "seconds")
  
  paste0("  The ", pt, " process took ", trunc(t1), " ", word1, ", ", 
         trunc(t2), " ", word2, " and ", trunc(t3), " ", word3, ".\n") %>% return()
}

outputDetails <- function(msg, newcat = FALSE) {
  cat(msg)
  if (newcat) {cat("\n")}
  message(msg)
}

createID <- function(df, c0, c1, c2) {
  df %>% add_column(id = paste0(toupper(c0), "_h", pull(df, c1), "_c", pull(df, c2))) %>% return()
}

rmIDComp <- function(df) {df %>% select(-id, -composition) %>% return()}

meltData <- function(dataset, id_val) {melt(dataset, id = id_val) %>% as_tibble() %>% return()}

factorToInt <- function(dataset, cname) {
  dataset %>% mutate(across(all_of(cname), as.character)) %>% 
    mutate(across(all_of(cname), as.integer)) %>% return()
}

charToInt <- function(x, v) {
  gsub(v, "", x) %>% as.integer()
}

codeIsolates <- function(df, tpx, all_iso) {
  hx <- paste0(tpx, "_h")
  cx <- paste0(tpx, "_cl")
  
  df %>% right_join(all_iso, ., by = c("char_isolate" = "isolate")) %>% 
    select(-char_isolate) %>% rename(isolate = num_isolate) %>% melt(id = "isolate") %>% 
    as_tibble() %>% set_colnames(c("isolate", hx, cx)) %>% factorToInt(., hx) %>% 
    createID(., tpx, hx, cx) %>% mutate(isolate = paste0("-", isolate, "-")) %>% return()
}

meltedSizing <- function(df, y) {
  tp <- paste0("tp", y, "_")
  # the TPX cluster assignments
  tp_melted <- df %>% melt(id = "isolate") %>% as_tibble() %>% 
    factorToInt(., "variable") %>% 
    createID(., paste0("tp", y), "variable", "value") %>% 
    set_colnames(c("isolate", "tp_h", "tp_cl", "tp_id"))
  
  # adding a column with the TPX cluster sizes
  tp_melted %>% group_by(tp_id) %>% summarise(tp_cl_size = n(), .groups = "drop") %>% 
    left_join(tp_melted, ., by = "tp_id") %>% 
    set_colnames(c("isolate", paste0(tp, "h"), paste0(tp, "cl"), 
                   paste0(tp, "id"), paste0(tp, "cl_size"))) %>% return()
}

leadingZeros <- function(df, cname, lc, w = NULL) {
  if (is.null(w)) {w <- df[,cname] %>% max() %>% nchar()}
  df[,cname] <- pull(df, cname) %>% formatC(., width = w, format = "d", flag = "0") %>% paste0(lc, .)
  return(df)
}

# given the defining filename, read in the data (need the full path from your working directory)
readBaseData <- function(filename, file_number) {
  if (is.na(filename)) {
    stop(paste0("Time point ", file_number, " dataset not found."))
  }else {
    read.csv(file = filename, stringsAsFactors = FALSE, numerals = "no.loss", 
             check.names = FALSE, sep = "\t") %>% as_tibble() %>% return()
  }
}

resultFiles <- function(df, op, heights, time1_raw) {
  clusters_formatted <- df %>% set_colnames(c(
    "TP1 ID", "TP1 height", "TP1 cluster", "TP1 cluster size", "First time this cluster was seen in TP1", 
    "Last time this cluster was seen in TP1", "First time this cluster was seen in TP2", "TP2 height", 
    "TP2 cluster", "TP2 cluster size", "Number of additional TP1 strains in the TP2 match", 
    "Number of novels in the TP2 match", "Actual cluster size change (TP2 size - TP1 size)",
    "Actual growth rate = (TP2 size - TP1 size) / (TP1 size)", "Number of novels / Actual growth rate", 
    "Novel growth = (TP2 size) / (TP2 size - number of novels)"))
  
  write.csv(clusters_formatted, paste0(op,"TP1_cluster_results.csv"), row.names = FALSE)
  
  isolates_formatted <- time1_raw %>% select(isolate, as.character(heights)) %>% 
    melt(., id = "isolate") %>% as_tibble() %>% 
    rename(tp1_h = variable, tp1_cl = value) %>% 
    mutate(across(tp1_h, as.character)) %>% mutate(across(tp1_h, as.integer)) %>% 
    leadingZeros(., "tp1_h", "h") %>% leadingZeros(., "tp1_cl", "c") %>% 
    right_join(., df, by = c("tp1_h", "tp1_cl")) %>% 
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% select(isolate, colnames(df)) %>% 
    set_colnames(c("Isolates", colnames(clusters_formatted)))
  
  write.csv(isolates_formatted, paste0(op, "TP1_strain_results.csv"), row.names = FALSE)
}

saveData <- function(tmp = NULL, h = NULL) {
  if (!dir.exists("outputs")) {dir.create("outputs/height_data", recursive = TRUE)}
  saveRDS(tmp, paste0("outputs/height_data/h", h, ".Rds"))
}

meltedIDs <- function(df, k) {
  cnames <- paste0("tp", k, c("", "_h", "_cl", "_id"))
  df %>% 
    melt(id = "isolate") %>% as_tibble() %>% 
    set_colnames(c("isolate", cnames[2:3])) %>% 
    mutate(across(cnames[2], as.character)) %>% 
    createID(., cnames[1], cnames[2], cnames[3]) %>% 
    set_colnames(c("isolate", cnames[2:4])) %>% return()
}