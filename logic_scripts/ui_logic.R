navbarPage(
  tags$h2("My secure application"),
  verbatimTextOutput("auth_output")
)

ui <- fluidPage(theme = shinytheme("flatly"), shinythemes::themeSelector(),

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

        # Annotations row
        wellPanel(fluidRow(
          column(4,
            selectizeInput("gene", "Choose a gene for expression:",
              choices = NULL)
          ),
          column(4,
            radioButtons(
              "RNA_exprs_anno", "Colour annotation:", inline = TRUE,
              choices = c(
                "Population" = 1,
                "Stage" = 2,
                "Sample" = 3)
            )
          ),
          column(4,
            sliderInput("RNA_topn",
              "Top variable genes for PCA calculation:",
              min = 0,
              max = length(pull(tbl(sql_db, "RNA_exprs_wide"), GeneID)),
              value = length(pull(tbl(sql_db, "RNA_exprs_wide"), GeneID)))
          ),
        )),

        # Plots row
        fluidRow(
          column(6, h3("Expression plot"),
            plotOutput("RNA_exprs") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE)),

          column(6, h3("PCA plot"),
            plotOutput("RNA_pca") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE))
        ),

        # Stats row
        fluidRow(
          column(8, h3("Differential expression analysis:")),
          DT::dataTableOutput("RNA_stats_tbl", width = "90%") %>%
            withSpinner(type = 5, color = "#0dc5c1")
        )
      ),

      ### ATAC panel code ###
      tabPanel("ATAC",

        # Annotations row
        wellPanel(fluidRow(
          column(4,
            selectizeInput("enhancer",
              "Choose an enhancer for accessibility plot:",
              choices = NULL)
          ),
          column(4,
            radioButtons(
              "ATAC_exprs_anno", "Colour annotation:", inline = TRUE,
              choices = c(
                "Population" = 1,
                "Stage" = 2,
                "Sample" = 3)
            )
          ),
          column(4,
            sliderInput("ATAC_topn",
              "Top variable enhancers for PCA calculation:\n
              (Max 10,000 due to memory limits)",
              min = 0,
              max = 10000, #length(pull(tbl(sql_db, "ATAC_exprs_wide"), peak_coord)),
              value = 5000)
          ),
        )),

        # Plots row
        fluidRow(
          column(6, h3("Accessibility plot"),
            plotOutput("ATAC_exprs") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE)),

          column(6, h3("PCA plot"),
            plotOutput("ATAC_pca") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE))
        ),

        # Stats row
        fluidRow(
          column(8, h3("Differential accessibility analysis:")),
          DT::dataTableOutput("ATAC_stats_tbl", width = "90%") %>%
            withSpinner(type = 5, color = "#0dc5c1")
        )
      )
    ), width = 10)
  )
)
