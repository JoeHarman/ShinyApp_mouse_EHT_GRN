library(shiny)
library(shinymanager)
library(tidyverse)

# Set up user passwords
credentials <- data.frame(
  user = c("jharman"), # mandatory
  password = c("eht_is_cool"), # mandatory
  start = c("2023-04-19"),
  expire = c(NA),
  admin = c(TRUE),
  comment = "Simple and secure authentification mechanism 
      for single ‘Shiny’ applications.",
  stringsAsFactors = FALSE
)

# Def functions
p_formatter <- function(x){
  format(x, digits = 2, scientific = TRUE)
}

# Read in data
ATAC_exprs <- readRDS("./data/ATAC_exprs.rds")
RNA_exprs <- readRDS("./data/RNA_exprs.rds")
ATAC_stats <- readRDS("./data/ATAC_stats.rds") %>%
  mutate_at(c(4, 6, 8, 10, 12), function(x) round(x, 2)) %>%
  mutate_at(c(5, 7, 9, 11, 13), function(x) signif(x, 3))
RNA_stats <- readRDS("./data/RNA_stats.rds") %>%
  mutate_at(c(2, 4, 6, 8, 10), function(x) round(x, 2)) %>%
  mutate_at(c(3, 5, 7, 9, 11), function(x) signif(x, 3))

# Load ui and server scripts
ui <- source("./logic_scripts/ui_logic.R",  local = TRUE)$value
server <- source("./logic_scripts/server_logic.R",  local = TRUE)$value

# Create shiny app
shinyApp(secure_app(ui), server)