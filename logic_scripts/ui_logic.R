navbarPage(
  tags$h2("My secure application"),
  verbatimTextOutput("auth_output")
)

ui <- fluidPage(theme = shinytheme("flatly"), # shinythemes::themeSelector(),

  titlePanel("Gene regulatory network analysis of embryonic mouse EHT"),

  sidebarLayout(

    ##### SIDE BAR CODE #####
    sidebarPanel(
      checkboxGroupInput("select_groups", "Samples to plot:", samples,
        selected = samples[c(1, 2, 4, 5, 8, 9)]),
      br(),
      "Click for instructions:", br(),
      actionBttn("rna_help", "RNA tab", style = "stretch",
        color = "warning", size = "sm"), br(),
      actionBttn("atac_help", "ATAC tab", style = "stretch",
        color = "warning", size = "sm"), br(),
      actionBttn("net_help", "Network tab", style = "stretch",
        color = "warning", size = "sm"), br(),
      actionBttn("coop_help", "Cooperation tab", style = "stretch",
        color = "warning", size = "sm"),
      bsModal("rna_help_b", "RNA-seq analysis", "rna_help",
        size = "large", HTML(help1)),
      bsModal("atac_help_b", "ATAC-seq analysis", "atac_help",
        size = "large", HTML(help2)),
      bsModal("net_help_b", "Network analysis", "net_help",
        size = "large", HTML(help3)),
      bsModal("coop_help_b", "Cooperation analysis", "coop_help",
        size = "large", HTML(help4)),
    width = 2),

    ##### MAIN PANEL CODE #####
    mainPanel(tabsetPanel(
      ### RNA panel code ###
      tabPanel("About",
        br(),
        HTML(about_txt)
      ),
      ### RNA panel code ###
      tabPanel("RNA",

        # Annotations row
        wellPanel(fluidRow(
          column(4,
            selectizeInput("gene", "Choose a gene for expression:",
              choices = NULL)
          ),
          column(2,
            radioButtons(
              "RNA_exprs_anno", "Colour annotation:", inline = TRUE,
              choices = c(
                "Cell population" = 1,
                "Embryo stage" = 2,
                "Sample name" = 3)
            )
          ),
          column(4,
            sliderInput("RNA_topn",
              "Top variable genes for PCA calculation:",
              min = 0,
              max = length(pull(tbl(sql_db, "RNA_exprs_wide"), GeneID)),
              value = length(pull(tbl(sql_db, "RNA_exprs_wide"), GeneID))),
            bsTooltip("RNA_topn", tooltips[19])
          ),
          column(2,
            actionButton(inputId = "runPCA_RNA", label = "Run PCA (slow!)")
          )
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
          column(2,
            radioButtons(
              "ATAC_exprs_anno", "Colour annotation:", inline = TRUE,
              choices = c(
                "Cell population" = 1,
                "Embryo stage" = 2,
                "Sample name" = 3)
            )
          ),
          column(4,
            sliderInput("ATAC_topn",
              "Top variable enhancers for PCA calculation:\n
              (Max 10,000 due to memory limits)",
              min = 0,
              max = 10000,
              #length(pull(tbl(sql_db, "ATAC_exprs_wide"), peak_coord)),
              value = 5000),
            bsTooltip("ATAC_topn", tooltips[19])
          ),
          column(2,
            actionButton(inputId = "runPCA_ATAC", label = "Run PCA (slow!)")
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
      ),

      ### ATAC panel code ###
      tabPanel("Network",

        # Annotations row
        wellPanel(fluidRow(

          column(3,
            actionButton(inputId = "makeGRN", label = "Make GRN"),
            downloadButton("downGRN", "Export"),
            br(),
            radioGroupButtons(inputId = "grn_mode", label = "Plot mode:",
              choices = c("Central TFs", "Upstream", "Downstream")),
            radioTooltip("grn_mode", "Central TFs", tooltips[15]),
            radioTooltip("grn_mode", "Upstream", tooltips[16]),
            radioTooltip("grn_mode", "Downstream", tooltips[17]),
            bsTooltip("makeGRN", tooltips[18])
          ),
          column(3,
            radioButtons(
              "grn_subset",
              "GRN subset:",
              c("Full", "RNA1", "RNA2", "RNA3", "RNA4", "RNA5"),
              selected = "Full",
              inline = TRUE
            ),
            shinyWidgets::prettySwitch("grn_tf", label = "Show TFs only",
              value = TRUE, fill = TRUE),
            radioTooltip("grn_subset", "Full", tooltips[8]),
            radioTooltip("grn_subset", "RNA1", tooltips[9]),
            radioTooltip("grn_subset", "RNA2", tooltips[10]),
            radioTooltip("grn_subset", "RNA3", tooltips[11]),
            radioTooltip("grn_subset", "RNA4", tooltips[12]),
            radioTooltip("grn_subset", "RNA5", tooltips[13]),
            bsTooltip("grn_tf", tooltips[14])
          ),
          column(3,
            sliderInput("top_n_centrality",
              "# top-most central:",
              min = 1, max = 316, value = 25),
            radioButtons(
              "grn_centrality",
              "Centrality measure:",
              c(Degree = "degree", Betweenness = "betweenness",
                Closeness = "closeness", Eigenvector = "eigenvector",
                "Hub score" = "hub_score",
                "Authority score" = "authority_score"),
              selected = "degree",
              inline = TRUE
            ),
            # Tooltips
            radioTooltip("grn_centrality", "degree", tooltips[1]),
            radioTooltip("grn_centrality", "betweenness", tooltips[2]),
            radioTooltip("grn_centrality", "closeness", tooltips[3]),
            radioTooltip("grn_centrality", "eigenvector", tooltips[4]),
            radioTooltip("grn_centrality", "hub_score", tooltips[5]),
            radioTooltip("grn_centrality", "authority_score", tooltips[6]),
            bsTooltip("top_n_centrality", tooltips[7])

          ),
          column(3,
            selectizeInput("core_node", "Reference node:",
              choices = NULL),
            radioGroupButtons(inputId = "corr_filt",
              label = "Side plot options:",
              choices = c("Top-20 correlated", "Top-20 anti-correlated"))
          )
        )),

        # Plots row
        fluidRow(
          column(9, h3("Network"),
            visNetworkOutput("mynetworkid", height = "70vh") %>%
              withSpinner(type = 5, color = "#0dc5c1")
          ),

          column(3, h3("Top nodes"),
            plotOutput("networkSidePlot", height = "70vh") %>%
              withSpinner(type = 5, color = "#0dc5c1")
          )
        )
      ),

      ### Cooperation panel code ###
      tabPanel("Cooperation",

        # Annotations row
        wellPanel(fluidRow(
          column(3,
            selectizeInput("TF_A", "TF A:",
              choices = NULL), br(),
            selectizeInput("TF_B", "TF B:",
              choices = NULL),
            bsTooltip("TF_A", tooltips[20]),
            bsTooltip("TF_B", tooltips[21])
          ),
          column(3,
          ),
          column(3,
          ),
          column(3,
          )
        )),

        # Plots row
        fluidRow(
          column(9, h3("Network"),
            plotOutput("coop_exprs") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE)
          ),

          column(3, h3("Top nodes"),
            #plotOutput("networkSidePlot", height = "70vh") %>%
            #  withSpinner(type = 5, color = "#0dc5c1")
          )
        )
      )
    ), width = 10)
  )
)
