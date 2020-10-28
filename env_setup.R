

require("tibble")
require("magrittr")
require("dplyr")
require("reshape2")
require("scales")

dir.create("data", showWarnings = FALSE)

# Now move the data you want to use to the newly created "data" directory

checkDirInput <- function() {
  cat(paste0("Please move any data files you\'d like ",
             "to use to the newly created \'data\' directory. ",
             "Hit Enter to continue."))
  msg = readLines(con = "stdin", 1)
}

checkDirInput()
