x <- c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress", "reader", 
       "stringr", "ggplot2", "plotly", "optparse", "methods", "R6", "testit")
lapply(x, require, character.only = TRUE)

flaggingClusters <- function(tp_comps, tpx) {
  # Identifying the first and last time each TP1 cluster was seen in TP1
  t_fal <- tp_comps %>% group_by(composition) %>% slice(1, n()) %>% 
    add_column(type = rep(c("first", "last"), nrow(.)/2)) %>% dplyr::ungroup()

  t_ff <- t_fal %>% filter(type == "first") %>% select(grep("id", colnames(.)), composition)
  colnames(t_ff) <- c(paste0("first_", tpx, "_flag"), "composition")
  
  t_lf <- t_fal %>% filter(type == "last") %>% select(grep("id", colnames(.)), composition)
  colnames(t_lf) <- c(paste0("last_", tpx, "_flag"), "composition")
  
  full_join(t_ff, t_lf, by = "composition") %>% 
    left_join(tp_comps, ., by = "composition") %>% 
    select(-composition) %>% return()
}

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
newID <- function(df, tpx, c1, c2, ph, pc) {
  newh <- df %>% pull(c1) %>% as.character() %>% as.integer() %>% 
    formatC(., width = ph, format = "d", flag = "0") %>% paste0("h", .)
  newc <- df %>% pull(c2) %>% as.character() %>% as.integer() %>% 
    formatC(., width = pc, format = "d", flag = "0") %>% paste0("c", .)
  df %>% add_column(id = paste0(toupper(tpx), "_", newh, "_", newc)) %>% return()
}

meltData <- function(dataset, id_val) {
  melt(dataset, id = id_val) %>% as_tibble() %>% return()
}

# Given a dataframe df, converts all the elements of column c1 from factors to integers
factorToInt <- function(df, c1) {
  df %>% mutate(across(all_of(c1), as.character)) %>% 
    mutate(across(all_of(c1), as.integer)) %>% return()
}

# Given a raw time point dataset, the timepoint ID (e.g. "tp1"), and a list of all isolates 
# present at TP1 and TP2, return a dataframe with isolates given numeric codes (e.g. "-1-", "-10-")
codeIsolates <- function(df, tpx, all_iso, ph, pc) {
  hx <- paste0(tpx, "_h")
  cx <- paste0(tpx, "_cl")
  idx <- paste0(tpx, "_id")
  
  df %>% right_join(all_iso, ., by = c("char_isolate" = "isolate")) %>% 
    select(-char_isolate) %>% 
    rename(isolate = num_isolate) %>% 
    melt(id = "isolate") %>% as_tibble() %>% 
    set_colnames(c("isolate", hx, cx)) %>% 
    factorToInt(., hx) %>% 
    newID(., tpx, hx, cx, ph, pc) %>% 
    mutate(isolate = paste0("-", isolate, "-")) %>% 
    rename_with(., ~gsub("id", idx, .x)) %>% return()
  
}

# Function used in manual_check_results.R, given raw data for a timepoint and a lowercase 
# timepoint ID, return a dataframe with a column containing cluster sizes
meltedSizing <- function(df, y, ph, pc) {
  tp <- paste0(y, "_")
  tp_melted <- df %>% melt(id = "isolate") %>% as_tibble() %>% 
    factorToInt(., "variable") %>% 
    newID(., y, "variable", "value", ph, pc) %>% 
    # createID(., y, "variable", "value") %>% 
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
readBaseData <- function(filename, file_number, delimiter) {
  if (is.na(filename)) {
    stop(paste0("Time point ", file_number, " dataset not found."))
  }else {
    read.table(file = filename, stringsAsFactors = FALSE, check.names = FALSE, 
               header = TRUE, sep = delimiter, allowEscapes = TRUE) %>% as_tibble() %>% return()
  }
}

padCol <- function(colvalues, padval, padchr) {
  ifelse(!is.na(colvalues), 
         formatC(colvalues, width = padval, format = "d", flag = "0") %>% 
           paste0(padchr, .), NA) %>% return()
}

meltedIDs <- function(df, k, ph, pc) {
  cnames <- paste0(k, c("", "_h", "_cl", "_id"))
  df %>% melt(id = "isolate") %>% as_tibble() %>% 
    set_colnames(c("isolate", cnames[2:3])) %>% 
    mutate(across(cnames[2], as.character)) %>% 
    newID(., cnames[1], cnames[2], cnames[3], ph, pc) %>% 
    set_colnames(c("isolate", cnames[2:4])) %>% return()
}

convertAndSave <- function(ip, op) {
  df <- read.csv(ip, stringsAsFactors = FALSE, sep = ",", numerals = "no.loss") %>% as_tibble()
  m1 <- ncol(df)-2
  df %>% set_colnames(c("isolate", 0:m1)) %>% 
    write.table(., op, row.names = FALSE, quote = FALSE, sep = "\t")
}
