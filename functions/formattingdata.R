x <- c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress", 
       "stringr", "ggplot2", "plotly", "optparse", "methods", "R6", "rlist")
lapply(x, require, character.only = TRUE)

# Indicates length of a process in hours, minutes, and seconds, when given a name of the process 
# ("pt") and a two-element named vector with Sys.time() values named "start_time" and "end_time"
timeTaken <- function(pt, sw) {
  z <- difftime(sw[['end_time']], sw[['start_time']], units = "secs") %>% as.double()
  m <- 60
  h <- m^2
  
  if (z >= h) {
    hrs <- trunc(z/h)
    mins <- trunc(z/m - hrs*m)
    paste0("The ", pt, " process took ", hrs, " hour(s), ", mins, " minute(s), and ", 
           round(z - hrs*h - mins*m), " second(s).") %>% return()
  }else if (z < h & z >= m) {
    mins <- trunc(z/m)
    paste0("The ", pt, " process took ", mins, " minute(s) and ", round(z - mins*m), " second(s).") %>% return()
  }else {
    paste0("The ", pt, " process took ", round(z), " second(s).") %>% return()
  }
}

# Outputs the same message in two ways, one is directed to stdout and one to a log file
outputDetails <- function(msg, newcat = FALSE) {
  cat(msg)
  if (newcat) {cat("\n")}
  message(msg)
}

# Given a dataframe df, two column names c1 and c2 (height and cluster respectively) and a new
# ID prefix tpx (e.g. "tp1"), creates an ID column and adds to df before returning df
createID <- function(df, tpx, c1, c2) {
  df %>% add_column(id = paste0(toupper(tpx), "_h", pull(df, c1), "_c", pull(df, c2))) %>% return()
}

rmIDComp <- function(df) {
  df %>% select(-id, -composition) %>% return()
}

meltData <- function(dataset, id_val) {
  melt(dataset, id = id_val) %>% as_tibble() %>% return()
}

# Given a dataframe df, converts all the elements of column c1 from factors to integers
factorToInt <- function(df, c1) {
  df %>% mutate(across(all_of(c1), as.character)) %>% 
    mutate(across(all_of(c1), as.integer)) %>% return()
}

charToInt <- function(x, v) {
  gsub(v, "", x) %>% as.integer()
}

# Given a raw time point dataset, the timepoint ID (e.g. "tp1"), and a list of all isolates 
# present at TP1 and TP2, return a dataframe with isolates given numeric codes (e.g. "-1-", "-10-")
codeIsolates <- function(df, tpx, all_iso) {
  hx <- paste0(tpx, "_h")
  cx <- paste0(tpx, "_cl")
  
  df %>% right_join(all_iso, ., by = c("char_isolate" = "isolate")) %>% select(-char_isolate) %>% 
    rename(isolate = num_isolate) %>% 
    melt(id = "isolate") %>% as_tibble() %>% 
    set_colnames(c("isolate", hx, cx)) %>% 
    factorToInt(., hx) %>% 
    createID(., tpx, hx, cx) %>% 
    mutate(isolate = paste0("-", isolate, "-")) %>% return()
}

# Function used in manual_check_results.R, given raw data for a timepoint and a lowercase 
# timepoint ID, return a dataframe with a column containing cluster sizes
meltedSizing <- function(df, y) {
  tp <- paste0(y, "_")
  tp_melted <- df %>% melt(id = "isolate") %>% as_tibble() %>% 
    factorToInt(., "variable") %>% 
    createID(., y, "variable", "value") %>% 
    set_colnames(c("isolate", "tp_h", "tp_cl", "tp_id"))
  
  tp_melted %>% 
    group_by(tp_id) %>% 
    summarise(tp_cl_size = n(), .groups = "drop") %>% 
    left_join(tp_melted, ., by = "tp_id") %>% 
    set_colnames(c("isolate", paste0(tp, "h"), paste0(tp, "cl"), 
                   paste0(tp, "id"), paste0(tp, "cl_size"))) %>% return()
}

# Adds leading zeros to column cname in dataframe df, with a prefix of a leading character (lc) 
# followed by w zeros (if no number is provided, the max number of characters in the column is used)
leadingZeros <- function(df, cname, lc, w = NULL) {
  if (is.null(w)) {
    w <- df[,cname] %>% max() %>% nchar()
  }
  df[,cname] <- pull(df, cname) %>% formatC(., width = w, format = "d", flag = "0") %>% paste0(lc, .)
  return(df)
}

# given the defining filename, read in the data (need the full path from your working directory), 
# indicate to user if file is not found
readBaseData <- function(filename, file_number, separator) {
  if (is.na(filename)) {
    stop(paste0("Time point ", file_number, " dataset not found."))
  }else {
    read.csv(file = filename, stringsAsFactors = FALSE, numerals = "no.loss", 
             check.names = FALSE, sep = separator) %>% as_tibble() %>% return()
  }
}

# Given a dataframe df, an output path op, the heights analyzed, and two TP1 datasets (one with the raw 
# data and one in melted form), save the cluster file and strain file to the output directory
# Note: write permissions required
resultFiles <- function(df, op, heights, time1_raw, t1_melted) {
  clusters_formatted <- df %>% set_colnames(c(
    "TP1 ID", "TP1 height", "TP1 cluster", "First time this cluster was seen in TP1", 
    "Last time this cluster was seen in TP1", "First time this cluster was seen in TP2", 
    "TP2 height", "TP2 cluster", "TP1 cluster size", "TP2 cluster size", 
    "Number of additional TP1 strains in the TP2 match", "Number of novels in the TP2 match", 
    "Actual cluster size change (TP2 size - TP1 size)",
    "Actual growth rate = (TP2 size - TP1 size) / (TP1 size)", 
    "Novel growth = (TP2 size) / (TP2 size - number of novels)"))
  
  write.table(clusters_formatted, file.path(op,"CGM_cluster_results.txt"), 
              row.names = FALSE, quote = FALSE, sep = "\t")
  
  m1 <- t1_melted$tp1_h %>% as.integer() %>% max() %>% nchar()
  m2 <- t1_melted$tp1_cl %>% as.integer() %>% max() %>% nchar()
  
  isolates_formatted <- time1_raw %>% 
    select(isolate, as.character(heights)) %>% 
    melt(., id = "isolate") %>% as_tibble() %>% 
    rename(tp1_h = variable, tp1_cl = value) %>% 
    mutate(across(tp1_h, as.character)) %>%
    mutate(across(tp1_h, as.integer)) %>% 
    leadingZeros(., "tp1_h", "h", m1) %>% 
    leadingZeros(., "tp1_cl", "c", m2) %>% 
    left_join(., df, by = c("tp1_h", "tp1_cl")) %>% 
    arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
    select(isolate, colnames(df)) %>% 
    set_colnames(c("Isolates", colnames(clusters_formatted)))
  
  write.table(isolates_formatted, file.path(op, "CGM_strain_results.txt"), 
              row.names = FALSE, quote = FALSE, sep = "\t")
}

meltedIDs <- function(df, k) {
  cnames <- paste0(k, c("", "_h", "_cl", "_id"))
  df %>% melt(id = "isolate") %>% as_tibble() %>% 
    set_colnames(c("isolate", cnames[2:3])) %>% 
    mutate(across(cnames[2], as.character)) %>% 
    createID(., cnames[1], cnames[2], cnames[3]) %>% 
    set_colnames(c("isolate", cnames[2:4])) %>% return()
}

convertAndSave <- function(ip, op) {
  df <- read.csv(ip, stringsAsFactors = FALSE, sep = ",", numerals = "no.loss") %>% as_tibble()
  m1 <- ncol(df)-2
  df %>% set_colnames(c("isolate", 0:m1)) %>% 
    write.table(., op, row.names = FALSE, quote = FALSE, sep = "\t")
}

# assertthat::see_if(identical(a1, a2), msg = paste0("File ", j, " is not the same in both directories."))
