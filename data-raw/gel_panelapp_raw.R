# Script to keep a raw version of the GEL panelapp NDD genes

library(readr)
library(dplyr)
library(stringr)
library(purrr)

# Function to load and add the panel ID column
load_and_clean_panel <- function(url) {
  data <- read_delim(url, delim = "\t", col_names = TRUE,
                     col_types = cols(.default = "c"))
  # Extract the panel id from the URL using regex
  panel_id <- str_extract(url, "(?<=panels/)[0-9]+")

  # Add the panel id as a new column
  data <- data %>% mutate(panel_id = panel_id)
  data
}

panel_urls <- c(
  "https://panelapp.genomicsengland.co.uk/panels/285/download/01234/",
  "https://panelapp.genomicsengland.co.uk/panels/197/download/01234",
  "https://panelapp.genomicsengland.co.uk/panels/78/download/01234/",
  "https://panelapp.genomicsengland.co.uk/panels/96/download/01234/"
)

# Load all panels and bind them into one data frame
panelapp <- map_dfr(panel_urls, load_and_clean_panel)

write.csv(panelapp, "data-raw/gel_panelapp_raw.csv", row.names = FALSE)



