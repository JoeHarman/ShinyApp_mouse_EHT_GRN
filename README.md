# ShinyApp_EHT-GRN

This is a Shiny App code to accompany the publication [[publication reference & link]]. This app offers a user interface to browse the related datasets as well as a gene regulatory network (GRN) model of EHT. A running version of this app is available at [[link]], and can be run locally.

Shiny app code author: Joe Harman

## About the data
The underlying data (available at <a href="https://doi.org/10.6084/m9.figshare.22905005.v1" target="_blank">DOI: 10.6084/m9.figshare.22905005.v1</a> and <a href="" target="_blank">GSE......</a>) consist of mini-bulk RNA-seq and ATAC-seq analyses for embryonic mouse EHT populations. These consisted of E8.5 to E10.5 embryos, on endothelial cells (EC), hemogenic endothelium (HE), pro-HSC, and pre-HSC type I and type II. Samples used in the publication, and upon which statistics were calculated, consisted of E8.5 EC, E8.5 HE, E9.5 HE, E9.5 pro-HSC, E10.5 pre-HSC type I, and E10.5 pre-HSC type II.

The GRN model describing EHT was generated through annotating differentially accessible ATAC-seq peaks to differentially expressed genes (RNA-seq), then connecting gene expression to upstream regulators via motif analysis of the cis-regulatory elements. These links were robustly filtered using correlative measures. Please see publication for details.

## Running locally
To run this Shiny app locally, please install R and the packages as listed in app.R, and clone the git repository to your computer. To initiate the app, run the following lines of code in R:

```
library(shiny)
runApp("path/to/shiny_directory/")
```
