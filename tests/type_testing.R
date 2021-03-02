suppressWarnings(suppressPackageStartupMessages(source("functions/tracking_functions.R")))
source("class_definitions.R")
cgm_all <- readRDS("outputs/isolates.Rds")

# f1 <- readBaseData(arg$tp1, 1, "\t")#arg$delimiter)
# f2 <- readBaseData(arg$tp2, 2, "\t")#arg$delimiter)
# heights <- strsplit(arg$heights, split = ",") %>% unlist()

f1 <- readBaseData("t1_clusters_processed.csv", 1, "\t")
f2 <- readBaseData("t2_clusters_processed.csv", 2, "\t")

hchars <- cgm_all %>% filter(!is.na(tp1_h)) %>% pull(tp1_h) %>% unique()
heights <- hchars %>% gsub("h", "", .) %>% as.integer() %>% as.character()

colnames(f1)[1] <- colnames(f2)[1] <- "strain"
padding_heights <- max(nchar(colnames(f1)[-1]), nchar(colnames(f2)[-1]))
padding_clusters <- f1 %>% select(-strain) %>% max(., f1 %>% select(-strain)) %>% nchar()

x2 <- f1 %>% melt(id = "strain") %>% as_tibble() %>% 
  newID(., "tp1", "variable", "value", padding_heights, padding_clusters)
x3 <- x2 %>% group_by(id) %>% summarise(size = n()) %>% left_join(x2, ., by = "id") %>% 
  rename(tp2_h = variable, tp2_cl = value) %>% mutate(across(tp2_h, as.character)) %>% 
  mutate(across(tp2_h, as.integer)) %>% arrange(tp2_h, tp2_cl)

# In the following tests, we go through the indicated columns to verify that each contains the expected values
# IMPORTANT NOTE: this test suite only works for a single TP1 threshold at a time (would need some refactoring 
# to work for multiple heights)

hx <- hchars[1]
hy <- hx %>% gsub("h", "", .) %>% as.integer() %>% as.character()
cgm <- cgm_all %>% filter(tp1_h %in% c(NA, hx))

# Checking strains all TP2 strains (for identity with results)
assert("All TP2 strains are found in results and vice versa", identical(sort(cgm$strain), sort(f2$strain)))


# Checking novels are noted properly in results
result_novels <- cgm %>% filter(novel == 1) %>% pull(strain) %>% sort()
actual_novels <- setdiff(f2$strain, f1$strain) %>% sort()
assert("Novels are indicated correctly", identical(result_novels, actual_novels))


# Checking TP1 IDs are appropriately labelled for novels and originals
res_tp1_ids <- cgm$tp1_id %>% str_split_fixed(., "_", 3) %>% data.frame(stringsAsFactors = FALSE) %>% 
  select(2:3) %>% set_colnames(c("tp1_h", "tp1_cl")) %>% as_tibble() %>% 
  add_column(strain = cgm$strain, .before = 1)

res_tp1_ids$tp1_h %<>% gsub("h", "", .) %>% as.integer()
res_tp1_ids$tp1_cl %<>% gsub("c", "", .) %>% as.integer()

actual_tp1_ids <- f1 %>% select(strain, all_of(hy)) %>% 
  melt(id = "strain") %>% as_tibble() %>% rename(tp1_h = variable, tp1_cl = value) %>% 
  mutate(across(tp1_h, as.character)) %>% mutate(across(tp1_h, as.integer)) %>% arrange(strain)

tp1_ids <- res_tp1_ids %>% filter(!is.na(tp1_h)) %>% arrange(strain)


# Checking TP1 IDs, heights, and clusters are correct for non-novels
assert(paste0("For non-novels, TP1 IDs refer to the cluster assignments at height ", hx), 
       identical(actual_tp1_ids, tp1_ids))


# Checking TP1 IDs are correct for novels
res_novel <- cgm %>% filter(strain %in% setdiff(f2$strain, f1$strain))
assert(paste0("For novels, TP1 IDs are NA"), all(is.na(res_novel$tp1_id)))


# Checking TP1 heights are correct for novels
assert(paste0("For novels, TP1 height are NA"), all(is.na(res_novel$tp1_h)))


# Checking TP1 clusters are correct for novels
assert(paste0("For novels, TP1 clusters are NA"), all(is.na(res_novel$tp1_cl)))


# Checking TP1 cluster sizes are correct for novels first found in fully novel TP2 clusters
fully_novel <- res_novel %>% filter(tp2_cl_size == num_novs)
assert(paste0("For novels first found in fully novel TP2 clusters, TP1 size is 0"), 
       all(fully_novel$tp1_cl_size == 0))


# Checking TP1 cluster sizes are correct for novels first found in mixed TP2 clusters
mixed_novel <- res_novel %>% filter(!(strain %in% fully_novel$strain))

cluster_to_int <- lapply(mixed_novel$first_tp2_flag, function(x) {
  strsplit(x, "_") %>% unlist() %>%  extract2(3) %>% gsub("c", "", .) %>% as.integer()
}) %>% unlist()
mixed_tp2 <- tibble(tp2_h = cluster_to_int) %>% 
  add_column(strain = mixed_novel$strain, .before = 1) %>% 
  add_column(res_tp1_cl_size = mixed_novel$tp1_cl_size, .after = 2)

tp1_sizes <- f2 %>% select(strain, all_of(hy)) %>% 
  rename(tp2_h = hy) %>% 
  filter(tp2_h %in% mixed_tp2$tp2_h) %>% 
  filter(!(strain %in% setdiff(f2$strain, f1$strain))) %>% 
  group_by(tp2_h) %>% summarise(tp1_cl_size = n())

size_check1 <- left_join(mixed_tp2, tp1_sizes, by = "tp2_h") %>% 
  mutate(across(c(res_tp1_cl_size, tp1_cl_size), as.integer))
assert(paste0("For novels first found in mixed clusters, they inherit the TP1 cluster's size"), 
       identical(size_check1$res_tp1_cl_size, size_check1$tp1_cl_size))


# Checking that the first and last times a cluster found in TP1 is flagged correctly
# for originals e.g. # idx <- "TP1_h0_c735"
tp1_clusters <- cgm %>% filter(!is.na(tp1_id)) %>% pull(tp1_id) %>% unique()
flag_check_t1 <- lapply(1:length(tp1_clusters), function(i) {
  idx <- tp1_clusters[i]
  c1 <- idx %>% strsplit(., "_") %>% unlist() %>% extract2(3) %>% gsub("c", "", .) %>% as.integer()
  c1strains <- f1 %>% select(strain, all_of(hy)) %>% 
    rename(tp1_h = hy) %>% filter(tp1_h == c1) %>% pull(strain)
  
  t1flags <- x3 %>% filter(strain %in% c1strains & size == length(c1strains)) %>% slice(1,n()) %>% pull(id) %>% 
    t() %>% data.frame(stringsAsFactors = FALSE) %>% set_colnames(c("first_tp1_flag", "last_tp1_flag")) %>% as_tibble()
  
  res_t1_flags <- cgm %>% filter(tp1_id == idx) %>% select(first_tp1_flag, last_tp1_flag) %>% slice(1)
  return(identical(res_t1_flags, t1flags))
}) %>% unlist() %>% set_names(tp1_clusters)
assert(paste0("For original TP1 clusters, the first and last TP1 clusters (i.e. when they formed ", 
              "and just before they changed, the range of stability) match the results"), 
       all(flag_check_t1))


# Checking that no novels are marked as being found in TP1 clusters
nov_first <- cgm %>% filter(novel == 1) %>% pull(first_tp1_flag)
assert(paste0("For novels, the first TP1 cluster each was found in is NA"), all(is.na(nov_first)))
nov_last <- cgm %>% filter(novel == 1) %>% pull(last_tp1_flag)
assert(paste0("For novels, the last TP1 cluster each was found in is NA"), all(is.na(nov_last)))


# Checking that for all novels, the first TP2 cluster they are found in is correctly marked
meltedtp2 <- f2 %>% melt(id = "strain") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% mutate(across(variable, as.integer)) %>% 
  rename(tp2_h = variable, tp2_cl = value) %>% 
  arrange(strain, tp2_h, tp2_cl)

novel_tp2_flags <- cgm %>% filter(novel == 1) %>% arrange(strain) %>% select(strain, first_tp2_flag)
  
novels <- setdiff(f2$strain, f1$strain)

nov_first_flags <- meltedtp2 %>% filter(strain %in% novels) %>% group_by(strain) %>% slice(1) %>% 
  ungroup() %>% newID(., "tp2", "tp2_h", "tp2_cl", padding_heights, padding_clusters) %>% 
  rename(first_tp2_flag = id) %>% select(strain, first_tp2_flag)

assert("For novels, the first TP2 cluster each was found in is correct", 
       identical(nov_first_flags, novel_tp2_flags))


# Checking that for all originals, the first TP2 cluster they are found in is correctly marked
full_tp2 <- f2 %>% melt(id = "strain") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% mutate(across(variable, as.integer)) %>% 
  rename(tp2_h = variable, tp2_cl = value) %>% 
  newID(., "tp2", "tp2_h", "tp2_cl", padding_heights, padding_clusters)
full_tp2 <- full_tp2 %>% group_by(id) %>% summarise(tp2_cl_size = n()) %>% 
  left_join(full_tp2, ., by = "id") %>% arrange(tp2_h, tp2_cl)

original_tp2_flags <- cgm %>% filter(novel == 0) %>% 
  select(strain, tp1_id, tp1_cl_size, first_tp2_flag) %>% arrange(tp1_id)

# SPEED UP THIS PART HERE - use group_by or something
orig_first_flags <- lapply(unique(original_tp2_flags$tp1_id), function(z) {
  strains_z <- original_tp2_flags %>% filter(tp1_id == z) %>% pull(strain)
  full_tp2 %>% filter(strain %in% strains_z) %>% 
    filter(tp2_cl_size >= length(strains_z)) %>% 
    slice(1) %>% select(id) %>% rename(first_tp2_flag = id) %>% 
    add_column(tp1_id = z, .before = 1) %>% return()
}) %>% bind_rows()

res_flags <- original_tp2_flags %>% select(tp1_id, first_tp2_flag) %>% unique()
assert("For all originals, the first TP2 cluster each was found in is correct", 
       identical(res_flags, orig_first_flags))


# Checking that all first_tp2_flags were accounted for
act_tp2_flags <- c(orig_first_flags$first_tp2_flag, novel_tp2_flags$first_tp2_flag) %>% unique() %>% sort()
res_tp2_flags <- cgm$first_tp2_flag %>% unique() %>% sort()
assert("All first TP2 flags are accounted for", identical(act_tp2_flags, res_tp2_flags))


# Checking that the tp2_h, tp2_cl columns are directly from "first_tp2_flag"
res_tp2_data <- cgm %>% pull(first_tp2_flag) %>% unique()
x6 <- paste0("TP2_", cgm$tp2_h, "_", cgm$tp2_cl) %>% unique()
assert("TP2 height and cluster columns directly correspond to the first_tp2_flag column", identical(res_tp2_data, x6))


# Checking that the TP2 cluster size is correct, considering the first TP2 flag as the TP2 ID
dfx <- cgm %>% select(tp2_cl, tp2_cl_size) %>% unique()
dfx$tp2_cl %<>% gsub("c", "", .) %>% as.integer()
dfy <- f2 %>% select(strain, all_of(hy)) %>% rename(tp2_cl = hy) %>% 
  group_by(tp2_cl) %>% summarise(actual_size = n()) %>% left_join(dfx, ., by = "tp2_cl")
assert("TP2 cluster size matches that of the actual TP2 clusters", identical(dfy$tp2_cl_size, dfy$actual_size))


# Checking actual_size_change column (should be TP2 cluster size - TP1 cluster size)
assert("Actual size change == TP2 cluster size - TP1 cluster size", 
       all(equals(cgm$tp2_cl_size - cgm$tp1_cl_size, cgm$actual_size_change)))


# Checking actual growth rate is (TP2 cluster size - TP1 cluster size) / (TP1 cluster size) for originals
originals <- cgm %>% filter(tp1_cl_size != 0)
new_calcs <- ((originals$tp2_cl_size - originals$tp1_cl_size) / originals$tp1_cl_size) %>% round(., digits = 3)
assert("Actual growth rate == (TP2 cluster size - TP1 cluster size) / (TP1 cluster size), for originals", 
       all(equals(originals$actual_growth_rate, new_calcs)))


# Checking actual growth rate is Inf for novels
assert("Growth rate for novels is denoted 'Inf' for infinite, since the denominator TP1 cluster size is 0", 
       cgm %>% filter(tp1_cl_size == 0) %>% pull(actual_growth_rate) %>% unique() %>% is.infinite())

# Checking novel growth was properly calculated for all originals
novel_growth <- (originals$tp2_cl_size / (originals$tp2_cl_size - originals$num_novs)) %>% round(., digits = 3)
assert("Novel growth == (TP2 cluster size / (TP2 cluster size - number of novels))", 
       identical(novel_growth, originals$new_growth))

# Checking novel growth was properly calculated for all novels
pure_novel <- cgm %>% filter(tp1_cl_size == 0)
assert("All novels first found in pure novel clusters at TP2 have 'Inf' for new growth (TP1 cluster size is 0)", 
       all(is.infinite(pure_novel$new_growth)))

# Checking that for all novels, the last TP2 cluster they are found in is correctly marked
# Checking that for all originals, the last TP2 cluster they are found in is correctly marked
# Checking that all last_tp2_flags were accounted for
# this test uses the fact that the first_tp2_flag was correct for all cases
assert("For novels, the first TP2 cluster each was found in is correct", identical(nov_first_flags, novel_tp2_flags))
assert("For all originals, the first TP2 cluster each was found in is correct", identical(res_flags, orig_first_flags))
assert("All first TP2 flags are accounted for", identical(act_tp2_flags, res_tp2_flags))

# given the first TP2 flags, track and find the height just before they change, in TP2
tp2_flags <- cgm %>% select(first_tp2_flag, last_tp2_flag)

tomelt <- f2 %>% melt(id = "strain") %>% as_tibble() %>% 
  mutate(across(variable, as.character)) %>% 
  mutate(across(variable, as.integer)) %>% 
  newID(., "tp2", "variable", "value", padding_heights, padding_clusters)

tp2_sizes <- tomelt %>% group_by(id) %>% 
  summarise(tp2_cl_size = n()) %>% 
  left_join(tomelt, tp2_sizes, by = "id")

# SPEED UP THIS PART
starter_strains <- lapply(1:length(tp2_flags$first_tp2_flag), function(i) {
  tp2_sizes %>% filter(id == tp2_flags$first_tp2_flag[i]) %>% pull(strain) %>% return()
}) %>% set_names(tp2_flags$first_tp2_flag)

# SPEED UP THIS PART
all_cases <- lapply(1:length(tp2_flags$first_tp2_flag), function(i) {
  tp2_sizes %>% filter(strain %in% starter_strains[[tp2_flags$first_tp2_flag[i]]]) %>% 
    filter(tp2_cl_size == length(starter_strains[[tp2_flags$first_tp2_flag[i]]])) %>% 
    slice(n()) %>% select(id) %>% return()
}) %>% bind_rows()

assert("All last_tp2_flags are correct, given that the first_tp2_flag is correct", 
       identical(all_cases$id, tp2_flags$last_tp2_flag))

# Check that the number of novels was correctly identified for every cluster
nov_check <- cgm %>% select(first_tp2_flag, num_novs)
just_novels <- tomelt %>% filter(strain %in% novels)
real_num_novs <- lapply(1:nrow(nov_check), function(i) {
  just_novels %>% filter(id == nov_check$first_tp2_flag[i]) %>% nrow() %>% return()
}) %>% unlist()

nov_check <- nov_check %>% add_column(real_num_novs) %>% 
  mutate(across(c(num_novs, real_num_novs), as.integer))
assert("The number of novels in each TP2 cluster was correctly identified", 
       identical(nov_check$num_novs, real_num_novs))

# Check that the number of sneakers was correctly identified for every cluster
# cgm %>% filter(tp1_cl_size == 0) %>% identical(., pure_novel)
assert(paste0("For all novels first found in TP2 in purely novel clusters, the ", 
              "number of 'sneakers' is 0"), all(pure_novel$add_TP1 == 0))
# for originals, should have real values, or be 0
cgm <- cgm_all %>% filter(tp1_h %in% c(NA, hchars[3]))

atp1 <- cgm %>% filter(novel == 0) %>% 
  select(tp1_id, first_tp2_flag, add_TP1) %>% unique()

part1 <- f1 %>% melt(id = "strain") %>% as_tibble() %>% 
  rename(tp1_h = variable, tp1_cl = value) %>% 
  mutate(across(tp1_h, as.character)) %>% 
  mutate(across(tp1_h, as.integer)) %>% 
  newID(., "tp1", "tp1_h", "tp1_cl", padding_heights, padding_clusters) %>% 
  rename(tp1_id = id) %>% 
  filter(!(strain %in% novels))

part2 <- f2 %>% melt(id = "strain") %>% as_tibble() %>% 
  rename(tp2_h = variable, tp2_cl = value) %>% 
  mutate(across(tp2_h, as.character)) %>% 
  mutate(across(tp2_h, as.integer)) %>% 
  newID(., "tp2", "tp2_h", "tp2_cl", padding_heights, padding_clusters) %>% 
  rename(tp2_id = id)

part3 <- part2 %>% filter(!(strain %in% novels))

part4 <- part1 %>% group_by(tp1_id) %>% summarise(tp1_size = n()) %>% ungroup() %>% 
  left_join(part1, ., by = "tp1_id")

check_sneakers <- lapply(1:nrow(atp1), function(i) {
  print(paste0(i, " / ", nrow(atp1)))
  case1 <- atp1[i,]
  a1 <- part1 %>% filter(tp1_id %in% case1$tp1_id) %>% select(strain)
  a2 <- part3 %>% filter(tp2_id %in% case1$first_tp2_flag) %>% select(strain)
  setdiff(a2, a1) %>% nrow() %>% return()
}) %>% unlist() %>% 
  add_column(atp1, actual = .)
assert("For all originals, the sneakers have been correctly counted", 
       all.equal(check_sneakers$add_TP1, check_sneakers$actual))
# for novels first found in mixed novel TP2 clusters, should inherit from TP1 data
# left_to_check <- cgm %>% filter(!(strain %in% pure_novel$strain)) %>% filter(novel != 0)
# identical(left_to_check, mixed_novel)
mixed_strains <- mixed_novel$strain
first_seen_tp2 <- part2 %>% filter(strain %in% mixed_strains) %>% 
  arrange(strain, tp2_h, tp2_cl) %>% group_by(strain) %>% slice(1) %>% ungroup() %>% 
  left_join(., mixed_novel[,c("strain", "add_TP1")])
# seen before (in TP1 cases)
a3 <- check_sneakers %>% select(-add_TP1)
sb <- first_seen_tp2[which(first_seen_tp2$tp2_id %in% check_sneakers$first_tp2_flag),] %>% 
  left_join(., a3, by = c("tp2_id" = "first_tp2_flag"))
assert(paste0("Novels that are first found in mixed TP2 clusters, whose TP1 clusters have ", 
              "been seen before, have sneakers counted accurately"), 
       all.equal(sb$add_TP1, sb$actual))
# not seen before
ns <- first_seen_tp2[which(!(first_seen_tp2$tp2_id %in% check_sneakers$first_tp2_flag)),]
actual_sizes <- ns$tp2_id %>% unique() %>% 
  lapply(., function(x) {
    stx <- part2 %>% filter(tp2_id %in% x) %>% filter(!(strain %in% novels)) %>% pull(strain)
    part4 %>% filter(strain %in% stx) %>% arrange(tp1_h, tp1_cl) %>% 
      filter(tp1_size == length(stx)) %>% slice(1) %>% pull(tp1_size) %>% return()
  }) %>% unlist() %>% tibble(tp2_id = unique(ns$tp2_id), tracked_tp1_size = .)

leftover_check <- cgm %>% filter(strain %in% ns$strain) %>% 
  select(strain, tp1_cl_size, add_TP1, tp2_cl_size, first_tp2_flag) %>% 
  rename(tp2_id = first_tp2_flag) %>% left_join(., actual_sizes)

assert(paste0("Novels first found in mixed TP2 clusters, whose TP1 clusters have not ", 
              "been seen before, have been tracked back to appropriate TP1 clusters ", 
              "(using current TP2)"), 
       all.equal(leftover_check$tp1_cl_size, leftover_check$tracked_tp1_size))
assert(paste0("Assuming the previous assert passed, we then know there must be ", 
              "no 'sneakers' in these clusters"), all(leftover_check$add_TP1 == 0))

cat("All cases checked.")
