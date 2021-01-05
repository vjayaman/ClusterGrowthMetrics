message("Going through the adaptive threshold development step.")

# ADAPTIVE THRESHOLD DEVELOPMENT

# We can plot what a couple of different threshold models can look like.
# This part will be heavily altered depending on what the actual data looks like.
# Currently being used to model synthetic data, and what we might expect significant
# size change to look like.

message("Running local regression steps for adaptive threshold development.")

x <- c(1, 25, 50, 100, 150)
y <- c(300, 150, 75, 30, 15)

# http://r-statistics.co/Loess-Regression-With-R.html
loessMod1 <- loess(y ~ x, span = 1)
smoothed1 <- predict(loessMod1)

loessMod5 <- loess(y ~ x, span = 5)
smoothed5 <- predict(loessMod5)

lm_df <- tibble(x, y, smoothed1, smoothed5) %>%
  set_colnames(c("x", "Preset values", "Local regression (span 1)", "Local regression (span 5)")) %>%
  meltData(., "x") %>% set_colnames(c("Cluster size", "Function", "Growth"))

message("Comparing models in the adaptive threshold step.")

# The model selection step, and using it to determine what the adaptive threshold
# function values would be for each input of initial cluster size. Anything that
# needs to be extrapolated is set to a pre-determined plateau value of percent increase.

model_used <- loessMod1
# predict.lm(model_used, data.frame(x = 30)) # note there are different types of prediction methods