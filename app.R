library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinymanager)
library(shinyWidgets)
library(shinythemes)
library(shinycssloaders)
library(shinyBS)
library(htmlwidgets)
library(readr)
library(dplyr)
library(dbplyr)
library(purrr)
library(tibble)
library(ggplot2)
library(DESeq2)
library(visNetwork)
library(igraph)

# Set up user passwords
credentials <- data.frame(
  user = c("jharman", "mdb_group"),
  password = c("eht_is_cool", "shiny_eht"),
  start = c("2023-04-19", "2023-04-19"),
  expire = c(NA, NA),
  admin = c(TRUE, FALSE),
  comment = "Basic password protection",
  stringsAsFactors = FALSE
)

# Source functions
source("./logic_scripts/functions.R")

# Read in data
ATAC_stats <- readRDS("./data/ATAC_stats.rds")
RNA_stats <- readRDS("./data/RNA_stats.rds")
RNA_anno <- readRDS("./data/RNA_anno.rds")
ATAC_anno <- readRDS("./data/ATAC_anno.rds")

pca_tables <- readRDS("./data/PCA_tables.rds")

sql_db <- DBI::dbConnect(RSQLite::SQLite(), dbname = "./data/SQL_DB.sqlite")

# Useful variables
samples <- readRDS("./data/vect_list.rds")$samples %>%
  as.character()
deg <- readRDS("./data/vect_list.rds")$deg %>%
  as.character()
dae <- readRDS("./data/vect_list.rds")$dae %>%
  as.character()

# Colour schemes
col_scheme <- list(
  Population = c(EC = "#4444c0", HE = "#53bf53",
    proHSC = "#daa84c", preI = "#c34747", preII = "#843b3b"),
  Stage = c(E8 = "#4444c0", E9 = "#53bf53", E10 = "#c34747"),
  Group = gg_color_hue(length(unique(pull(tbl(sql_db, "RNA_exprs"), Group))))
)
names(col_scheme$Group) <- unique(pull(tbl(sql_db, "RNA_exprs"), Group))

# Text
help1 <- read_file("./Text/help-1.txt")
help2 <- read_file("./Text/help-2.txt")
help3 <- read_file("./Text/help-3.txt")
help4 <- read_file("./Text/help-4.txt")
tooltips <- read_tsv("Text/tooltips.txt", col_names = FALSE, comment = "#") %>%
  unlist()

# Parameters
font_size <- 18

# Load ui and server scripts
ui <- source("./logic_scripts/ui_logic.R",  local = TRUE)$value
server <- source("./logic_scripts/server_logic.R",  local = TRUE)$value



# Create shiny app
shinyApp(secure_app(ui), server)