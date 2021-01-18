#! /usr/bin/env Rscript
input_args = commandArgs(trailingOnly = TRUE)

suppressWarnings(suppressPackageStartupMessages(source("funcs.R")))
library(ggplot2)
# time1_raw <- readBaseData("data/timepoint1_data.csv", 1)
# time2_raw <- readBaseData("data/timepoint2_data.csv", 2)
time1_raw <- readBaseData(input_args[1], 1)
time2_raw <- readBaseData(input_args[2], 1)

heights <- colnames(time1_raw)[-1]
# h_i <- '125'
h_i <- input_args[3]

# reading in the growth data for a particular height, and arranging by first TP1 then by TP2
b1 <- readData(dtype = 4, h = h_i) %>% arrange(tp1_h, tp1_cl, tp2_h, tp2_cl)

# extracting the first TP2_h_c match for each TP1 cluster at this height, and adding a column 
# indicating the actual size change (which may include non-novels)
c1 <- b1 %>% 
  pull(id) %>% unique() %>% 
  lapply(., function(idx) b1 %>% filter(id == idx) %>% slice(1)) %>% bind_rows()

c1 <- c1 %>% add_column(actual_size_change = c1$tp2_cl_size - c1$tp1_cl_size) %>% 
  select(tp1_h, tp1_cl, tp1_cl_size, tp2_h, tp2_cl, tp2_cl_size, 
         num_novs, actual_size_change, id, actual_growth_rate, b_ov_growth) %>% 
  rename(tp1_id = id) %>% 
  createID(., "tp2", "tp2_h", "tp2_cl") %>% 
  rename(tp2_id = id)
                    
c2 <- c1 %>% arrange(tp1_h, tp1_cl, tp2_h, tp2_cl) %>% 
  select(tp1_id, tp1_h, tp1_cl, tp1_cl_size, 
         tp2_id, tp2_h, tp2_cl, tp2_cl_size, 
         num_novs, actual_size_change, actual_growth_rate, b_ov_growth)

# checking a fishy case to see if it matches the generated values
# fishy_height <- '125'
# fishy_cluster <- 9
# cx <- time1_raw %>% 
#   select("isolate", all_of(fishy_height)) %>% 
#   rename(h1 = fishy_height) %>% 
#   filter(h1 == fishy_cluster)
# idx <- cx$isolate
# time2_raw %>% 
#   filter(isolate %in% idx) %>% 
#   melt(id = "isolate") %>% 
#   as_tibble() %>%
#   group_by(variable, value) %>% 
#   summarise(size = n(), .groups = "drop") %>%
#   filter(size >= nrow(cx))

c2 %>% 
  set_colnames(c("TP1 ID", "TP1 height", "TP1 cluster", "TP1 cluster size", 
                 "TP2 ID", "TP2 height", "TP2 cluster", "TP2 cluster size", 
                 "Number of novels in the TP2 match", 
                 "Actual size change (may incl. some TP1 isolates)", 
                 "Growth = (TP2 C.S. - TP1 C.S.) / (TP1 C.S)", 
                 "b / growth = (Number of novels) / (Growth rate)")) %>% 
  write.csv(., paste0("outputs/summary/TP1_h", h_i, "cluster_results.csv"), row.names = FALSE)

c3 <- time1_raw %>% select(isolate, all_of(h_i)) %>% 
  set_colnames(c("isolate", "tp1_cl")) %>% 
  right_join(., c2, by = "tp1_cl") %>% 
  arrange(isolate) %>% 
  select(isolate, tp1_id, tp1_h, tp1_cl, tp1_cl_size, tp2_id, tp2_h, tp2_cl, 
         tp2_cl_size, num_novs, actual_size_change, actual_growth_rate, b_ov_growth)

c3 %>% set_colnames(c("Isolates", "TP1 ID", "TP1 height", "TP1 cluster", "TP1 cluster size",
                      "TP2 ID", "TP2 height", "TP2 cluster", "TP2 cluster size",
                      "Number of novels in the TP2 match",
                      "Actual size change (may incl. some TP1 isolates)",
                      "Growth = (TP2 C.S. - TP1 C.S.) / (TP1 C.S)",
                      "b / growth = (Number of novels) / (Growth rate)")) %>%
  write.csv(., paste0("outputs/summary/TP1_h", h_i, "strain_results.csv"), row.names = FALSE)

text1 <- paste0("TP2 height: ", c3$tp2_h, "\nTP2 cluster: ", 
                c3$tp2_cl, "\nChange in size: ", c3$b_ov_growth)

g <- ggplot(c3, aes(x = tp1_cl, y = b_ov_growth, color = b_ov_growth)) + #, text = text1)) +
  geom_point() + xlab(paste0("TP1 clusters at height ", h_i)) +
  ylab(paste0("Cluster size change"))

ggsave(file.path("outputs/summary/size_plot.png"), device = "png",
       plot = g, width = 16, height = 9, units = "in", dpi = 320)
# ggplotly(g, tooltip = c("text" = "text1"))

# e.g. Rscript tmp.R data/timepoint1_data.csv 125