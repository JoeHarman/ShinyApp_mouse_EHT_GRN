library(shiny)
library(shinymanager)
library(shinythemes)
library(shinycssloaders)
library(htmlwidgets)
library(tidyverse)
library(edgeR)
library(DESeq2)

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

# Source functions
source("./logic_scripts/functions.R")

# Read in data
ATAC_exprs <- readRDS("./data/ATAC_exprs.rds")
RNA_exprs <- readRDS("./data/RNA_exprs.rds")
ATAC_stats <- readRDS("./data/ATAC_stats.rds")
RNA_stats <- readRDS("./data/RNA_stats.rds")
RNA_exprs_wide <- readRDS("./data/RNA_exprs_wide.rds")
ATAC_exprs_wide <- readRDS("./data/ATAC_exprs_wide.rds")
RNA_anno <- readRDS("./data/RNA_anno.rds")
ATAC_anno <- readRDS("./data/ATAC_anno.rds")

# Useful variables
samples <- unique(RNA_exprs$Group)
names(samples) <- samples
deg <- RNA_stats$GeneID
dae <- ATAC_stats$peak_coord

# Colour schemes
col_scheme <- list(
  Population = c(EC = "#4444c0", HE = "#53bf53",
    proHSC = "#daa84c", preI = "#c34747", preII = "#843b3b"),
  Stage = c(E8 = "#4444c0", E9 = "#53bf53", E10 = "#c34747"),
  Group = gg_color_hue(length(unique(ATAC_exprs$Group)))
)
names(col_scheme$Group) <- unique(ATAC_exprs$Group)

# Load ui and server scripts
ui <- source("./logic_scripts/ui_logic.R",  local = TRUE)$value
server <- source("./logic_scripts/server_logic.R",  local = TRUE)$value

# Create shiny app
shinyApp(secure_app(ui), server)