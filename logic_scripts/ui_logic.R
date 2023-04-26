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
    width = 2),

    ##### MAIN PANEL CODE #####
    mainPanel(tabsetPanel(

      ### RNA panel code ###
      tabPanel("RNA",

        # Stats row
        fluidRow(
          column(8, h3("Differential expression analysis:")),
          DT::dataTableOutput("RNA_stats_tbl", width = "90%")
        ),

        # Plots row
        fluidRow(
          column(6, h3("PCA plot"), plotOutput("RNA_pca")),
          column(6, h3("Expression plot"),
            selectInput("gene", "Choose a gene:",
            choices = unique(RNA_exprs$GeneID)),
            plotOutput("RNA_exprs"))
        ),

        # Annotations row
        fluidRow(
         # placeholder,
          radioButtons(
            "RNA_exprs_anno", "Colour annotation:", 
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
          DT::dataTableOutput("ATAC_stats_tbl", width = "90%")
        ),

        # Plots row
        fluidRow(
          column(6, h3("PCA plot"), renderPlot("ATAC_pca")),
          column(6, h3("Accessibility plot"),
            selectizeInput("enhancer", "Choose an enhancer:",
            choices = NULL),
            plotOutput("ATAC_exprs"))
        ),

        # Annotations row
        fluidRow(
        #  placeholder,
          radioButtons(
            "ATAC_exprs_anno", "Colour annotation:", 
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
