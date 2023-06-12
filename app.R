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
library(plotly)

### Retrieve data from figshare
if (!file.exists("./data/import.zip")) {
  download.file(
    url = "https://figshare.com/ndownloader/articles/22905005/versions/1",
    destfile = "./data/import.zip")
  unzip("./data/import.zip", exdir = "data")
}

### Source functions
source("./logic_scripts/functions.R")

### Read in data
ATAC_stats <- readRDS("./data/ATAC_stats.rds")
RNA_stats <- readRDS("./data/RNA_stats.rds")
RNA_anno <- readRDS("./data/RNA_anno.rds")
ATAC_anno <- readRDS("./data/ATAC_anno.rds")
pca_tables <- readRDS("./data/PCA_tables.rds")
coop_stats <- readRDS("./data/Cointeraction_stats.rds")

### Connect to SQLite database
sql_db <- DBI::dbConnect(RSQLite::SQLite(),
  dbname = "./data/RNA_ATAC_expression_and_GRN_tables.sqlite")

### Useful variables
samples <- as.character(unique(RNA_anno$Group))
deg <- as.character(RNA_stats$GeneID)
dae <- as.character(ATAC_stats$peak_coord)
font_size <- 15

### Colour schemes
col_scheme <- list(
  Population = c(EC = "#4444c0", HE = "#53bf53",
    `proHSC/HPC` = "#daa84c", preI = "#c34747", preII = "#843b3b"),
  Stage = c(E8 = "#4444c0", E9 = "#53bf53", E10 = "#c34747"),
  Group = gg_color_hue(length(unique(pull(tbl(sql_db, "RNA_exprs"), Group))))
)
names(col_scheme$Group) <- unique(pull(tbl(sql_db, "RNA_exprs"), Group))

### Text
about_txt <- read_file("./Text/about.txt")
help1 <- read_file("./Text/help-1.txt")
help2 <- read_file("./Text/help-2.txt")
help3 <- read_file("./Text/help-3.txt")
help4 <- read_file("./Text/help-4.txt")
tooltips <- read_tsv("Text/tooltips.txt", col_names = FALSE, comment = "#",
  show_col_types = FALSE) %>%
  unlist()

# Load ui and server scripts
ui <- source("./logic_scripts/ui_logic.R",  local = TRUE)$value
server <- source("./logic_scripts/server_logic.R",  local = TRUE)$value

# Create shiny app
shinyApp(ui, server)