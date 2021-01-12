#! /usr/bin/env Rscript
input_args = commandArgs(trailingOnly = TRUE)

source("functions/base_functions.R")
source("functions/processing_functions.R")
library(ggplot2)
# time1_raw <- readBaseData("data/timepoint1_data.csv", 1)
time1_raw <- readBaseData(input_args[1], 1)
heights <- colnames(time1_raw)[-1]
# h_i <- '125'
h_i <- input_args[2]
b1 <- readData(dtype = 4, h = h_i) %>% arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)

c1 <- b1 %>% arrange(tp1_h, tp1_cl) %>% pull(id) %>% unique() %>% 
  lapply(., function(idx) b1 %>% filter(id == idx) %>% slice(1)) %>% bind_rows()
c1$actual_size_change <- c1$tp2_cl_size - c1$tp1_cl_size

c2 <- c1 %>% select(tp1_h, tp1_cl, tp1_cl_size, 
                    tp2_h, tp2_cl, tp2_cl_size, 
                    num_novs, actual_size_change, id, 
                    actual_growth_rate, b_ov_growth)
c3 <- c2 %>% 
  add_column(tp1_id = paste0("TP1_h", h_i, "_c", c2$tp1_cl)) %>% 
  add_column(tp2_id = paste0("TP2_h", c2$tp2_h, "_c", c2$tp2_cl))

c4 <- c3 %>% arrange(-b_ov_growth) %>% 
  select(tp1_id, tp1_h, tp1_cl, tp1_cl_size, 
         tp2_id, tp2_h, tp2_cl, tp2_cl_size, 
         num_novs, actual_size_change, actual_growth_rate, b_ov_growth)

c4 %>% 
  set_colnames(c("TP1 ID", "TP1 height", "TP1 cluster", "TP1 cluster size", 
                 "TP2 ID", "TP2 height", "TP2 cluster", "TP2 cluster size", 
                 "Number of novels in the TP2 match", 
                 "Actual size change (may incl. some TP1 isolates)", 
                 "Growth = (TP2 C.S. - TP1 C.S.) / (TP1 C.S)", 
                 "b / growth = (Number of novels) / (Growth rate)")) %>% 
  write.csv(., paste0("outputs/summary/cluster_results_for_TP1_h", h_i, ".csv"), row.names = FALSE)

c5 <- time1_raw %>% select(isolate, all_of(h_i)) %>% 
  set_colnames(c("isolate", "tp1_cl")) %>% 
  right_join(., c4, by = "tp1_cl") %>% 
  arrange(isolate) %>% 
  select(isolate, tp1_id, tp1_h, tp1_cl, tp1_cl_size, tp2_id, tp2_h, tp2_cl, 
         tp2_cl_size, num_novs, actual_size_change, actual_growth_rate, b_ov_growth)

c5 %>% set_colnames(c("Isolates", "TP1 ID", "TP1 height", "TP1 cluster", "TP1 cluster size",
                      "TP2 ID", "TP2 height", "TP2 cluster", "TP2 cluster size",
                      "Number of novels in the TP2 match",
                      "Actual size change (may incl. some TP1 isolates)",
                      "Growth = (TP2 C.S. - TP1 C.S.) / (TP1 C.S)",
                      "b / growth = (Number of novels) / (Growth rate)")) %>%
  write.csv(., paste0("outputs/summary/strain_results_for_TP1_h", h_i, ".csv"), row.names = FALSE)

text1 <- paste0("TP2 height: ", c2$tp2_h, "\nTP2 cluster: ", 
                c2$tp2_cl, "\nChange in size: ", c2$b_ov_growth)

g <- ggplot(c3, aes(x = tp1_cl, y = b_ov_growth, color = b_ov_growth)) + #, text = text1)) +
  geom_point() + xlab(paste0("TP1 clusters at height ", h_i)) +
  ylab(paste0("Cluster size change"))

ggsave(file.path("outputs/summary/size_plot.png"), device = "png",
       plot = g, width = 16, height = 9, units = "in", dpi = 320)
# ggplotly(g, tooltip = c("text" = "text1"))

# e.g. Rscript tmp.R data/timepoint1_data.csv 125