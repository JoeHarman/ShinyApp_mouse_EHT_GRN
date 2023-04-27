navbarPage(
  tags$h2("My secure application"),
  verbatimTextOutput("auth_output")
)

ui <- fluidPage(

  titlePanel("EHT GRN app"),

  sidebarLayout(

    ##### SIDE BAR CODE #####
    sidebarPanel(
      checkboxGroupInput("select_groups", "Samples to plot:", samples,
        selected = samples[c(1, 2, 4, 5, 8, 9)]),
      actionButton(inputId = "subsetSamples", label = "Subset samples"),
    width = 2),

    ##### MAIN PANEL CODE #####
    mainPanel(tabsetPanel(

      ### RNA panel code ###
      tabPanel("RNA",

        # Stats row
        fluidRow(
          column(8, h3("Differential expression analysis:")),
          DT::dataTableOutput("RNA_stats_tbl", width = "90%") %>%
            withSpinner(type = 5, color = "#0dc5c1")
        ),

        # Plots row
        fluidRow(
          column(6, h3("PCA plot"),
            plotOutput("RNA_pca") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE)),

          column(6, h3("Expression plot"),
            selectInput("gene", "Choose a gene for expression:",
            choices = unique(deg)),
            plotOutput("RNA_exprs") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE))
        ),

        # Annotations row
        fluidRow(
          sliderInput("RNA_topn", "Top variable genes for PCA calculation:",
            min = 0, max = nrow(RNA_exprs_wide), value = 5000),
          radioButtons(
            "RNA_exprs_anno", "Colour annotation:", inline = TRUE,
            choices = c(
              "Population" = 1,
              "Stage" = 2,
              "Sample" = 3)
          ),
        )
      ),

      ### ATAC panel code ###
      tabPanel("ATAC",
        # Stats row
        fluidRow(
          column(8, h3("Differential accessibility analysis:")),
          DT::dataTableOutput("ATAC_stats_tbl", width = "90%") %>%
            withSpinner(type = 5, color = "#0dc5c1")
        ),

        # Plots row
        fluidRow(
          column(6, h3("PCA plot"),
            plotOutput("ATAC_pca") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE)),
          column(6, h3("Accessibility plot"),
            selectizeInput("enhancer",
              "Choose an enhancer for accessibility plot:",
              choices = NULL),
            plotOutput("ATAC_exprs") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE))
        ),

        # Annotations row
        fluidRow(
          sliderInput("ATAC_topn", "Number of observations:",
            min = 0, max = nrow(ATAC_exprs_wide), value = 1000),
          radioButtons(
            "ATAC_exprs_anno", "Colour annotation:", inline = TRUE,
            choices = c(
              "Population" = 1,
              "Stage" = 2,
              "Sample" = 3)
          )
        )
      )
    ), width = 10)
  )
)
