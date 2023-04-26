navbarPage(
  tags$h2("My secure application"),
  verbatimTextOutput("auth_output")
)

ui <- fluidPage(

  titlePanel("EHT GRN app"),

  sidebarLayout(

    ##### SIDE BAR CODE #####
    sidebarPanel(
        sliderInput("bins", label = "Number of bins:", min = 1, value = 30, max = 50),
    width=2),

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
        #  column(6,h3("Dimensionality Reduction:"),plotlyOutput("plot_graph", height=500)),
        #  column(6,h3("Bin. Expression BoxPlots:"),plotlyOutput("plot_graph2", height=500))
        ),

        # Annotations row
        fluidRow(
        #  placeholder,
        #  placeholder,
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
        #  placeholder,
        #  placeholder,
        ),

        # Annotations row
        fluidRow(
        #  column(6,h3("Dimensionality Reduction:"),plotlyOutput("plot_graph", height=500)),
        #  column(6,h3("Bin. Expression BoxPlots:"),plotlyOutput("plot_graph2", height=500))
        )
      )
    ), width = 10)
  )
)
