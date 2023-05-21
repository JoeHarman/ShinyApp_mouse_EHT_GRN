ui <- fluidPage(theme = shinytheme("flatly"),

  titlePanel("Gene regulatory network analysis of embryonic mouse EHT"),

  sidebarLayout(

    ##### SIDE BAR CODE #####
    sidebarPanel(
      ### Selecting samples for RNA/ATAC exprs plots
      checkboxGroupInput("select_groups", "Samples to plot:", samples,
        selected = samples[c(1, 2, 4, 5, 8, 9)]),
      br(),

      ### Instructions buttons
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

        wellPanel(fluidRow(
          ### Gene for expression plot
          column(4,
            selectizeInput("gene", "Choose a gene for expression:",
              choices = NULL)
          ),
          ### Options for colour annotation in exprs plot
          column(2,
            radioButtons(
              "RNA_exprs_anno", "Colour annotation:", inline = TRUE,
              choices = c(
                "Cell population" = 1,
                "Embryo stage" = 2,
                "Sample name" = 3)
            )
          ),
          ### Top_n items for PCA calculation
          column(4,
            sliderInput("RNA_topn",
              "Top variable genes for PCA calculation:",
              min = 0,
              max = length(pull(tbl(sql_db, "RNA_exprs_wide"), GeneID)),
              value = length(pull(tbl(sql_db, "RNA_exprs_wide"), GeneID))),
            bsTooltip("RNA_topn", tooltips[19])
          ),
          ### PCA execution button
          column(2,
            actionButton(inputId = "runPCA_RNA", label = "Run PCA")
          )
        )),

        # Plots row
        fluidRow(
          ### ATAC accessibility plot
          column(6, h3("Expression plot"),
            plotOutput("RNA_exprs") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE)),

          ### PCA plot
          column(6, h3("PCA plot"),
            plotOutput("RNA_pca") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE))
        ),

        ### Stats table
        fluidRow(
          column(8, h3("Differential expression analysis:")),
          DT::dataTableOutput("RNA_stats_tbl", width = "90%") %>%
            withSpinner(type = 5, color = "#0dc5c1")
        )
      ),

      ### ATAC panel code ###
      tabPanel("ATAC",

        wellPanel(fluidRow(
          ### Enhancer for accessibility plot
          column(4,
            selectizeInput("enhancer",
              "Choose an enhancer for accessibility plot:",
              choices = NULL)
          ),
          ### Options for colour annotation in access. plot
          column(2,
            radioButtons(
              "ATAC_exprs_anno", "Colour annotation:", inline = TRUE,
              choices = c(
                "Cell population" = 1,
                "Embryo stage" = 2,
                "Sample name" = 3)
            )
          ),
          ### Top_n items for PCA calculation
          column(4,
            sliderInput("ATAC_topn",
              "Top variable enhancers for PCA calculation:",
              min = 0,
              max = 10000,
              value = 5000),
            bsTooltip("ATAC_topn", tooltips[19])
          ),
          ### PCA execution button
          column(2,
            actionButton(inputId = "runPCA_ATAC", label = "Run PCA")
          ),
        )),

        fluidRow(
          ### ATAC accessibility plot
          column(6, h3("Accessibility plot"),
            plotOutput("ATAC_exprs") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE)),

          ### PCA plot
          column(6, h3("PCA plot"),
            plotOutput("ATAC_pca") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE))
        ),

        ### Stats table
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

          ### GRN mode/execute buttons
          column(3,
            # Execute GRN code
            actionButton(inputId = "makeGRN", label = "Make GRN"),
            # Download GRN
            downloadButton("downGRN", "Export"),
            br(),
            # Select plot mode
            radioGroupButtons(inputId = "grn_mode", label = "Plot mode:",
              choices = c("Central TFs", "Upstream", "Downstream")),
            # Tooltips
            radioTooltip("grn_mode", "Central TFs", tooltips[15]),
            radioTooltip("grn_mode", "Upstream", tooltips[16]),
            radioTooltip("grn_mode", "Downstream", tooltips[17]),
            bsTooltip("makeGRN", tooltips[18])
          ),

          ### GRN subset options
          column(3,
            # Module subset
            radioButtons(
              "grn_subset",
              "GRN subset:",
              c("Full", "RNA1", "RNA2", "RNA3", "RNA4", "RNA5"),
              selected = "Full",
              inline = TRUE
            ),

            # Transcription factor toggle
            shinyWidgets::prettySwitch("grn_tf", label = "Show TFs only",
              value = TRUE, fill = TRUE),

            # Tooltips
            radioTooltip("grn_subset", "Full", tooltips[8]),
            radioTooltip("grn_subset", "RNA1", tooltips[9]),
            radioTooltip("grn_subset", "RNA2", tooltips[10]),
            radioTooltip("grn_subset", "RNA3", tooltips[11]),
            radioTooltip("grn_subset", "RNA4", tooltips[12]),
            radioTooltip("grn_subset", "RNA5", tooltips[13]),
            bsTooltip("grn_tf", tooltips[14])
          ),

          ### Centrality options
          column(3,
            # Top-n central nodes to plot
            sliderInput("top_n_centrality",
              "# top-most central:",
              min = 1, max = 316, value = 25),

            # Method of centrality calculation
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

          ### Up/downstream network options
          column(3,
            # Reference node for up/downstream links
            selectizeInput("core_node", "Reference node:",
              choices = NULL),

            # Side-plot options - top correlated or anti-correlated
            radioGroupButtons(inputId = "corr_filt",
              label = "Side plot options:",
              choices = c("Top-20 correlated", "Top-20 anti-correlated"))
          )
        )),

        fluidRow(
          ### Plot network with visNetwork
          column(9, h3("Network"),
            visNetworkOutput("mynetworkid", height = "70vh") %>%
              withSpinner(type = 5, color = "#0dc5c1")
          ),

          ### Side plot (differs by GRN mode)
          column(3, h3("Top nodes"),
            plotOutput("networkSidePlot", height = "70vh") %>%
              withSpinner(type = 5, color = "#0dc5c1")
          )
        )
      ),

      ### Cooperation panel code ###
      tabPanel("Cooperation",

        wellPanel(fluidRow(
          ### TF A for cooperative analysis
          ### Acts as reference TF for cointeraction plot
          column(4,
            selectizeInput("TF_A", "TF A:",
              choices = NULL),
            bsTooltip("TF_A", tooltips[20])
          ),

          ### TF B for cooperative analysis
          column(4,
            selectizeInput("TF_B", "TF B:",
              choices = NULL),
            bsTooltip("TF_B", tooltips[21])
          ),

          ### Control RNA modules for cointeraction plot
          column(4,
            radioButtons(
              "coop_subset",
              "RNA module:",
              c("All", "RNA1", "RNA2", "RNA3", "RNA4", "RNA5"),
              selected = "All",
              inline = TRUE
            ),
            # Tooltips
            radioTooltip("coop_subset", "RNA1", tooltips[9]),
            radioTooltip("coop_subset", "RNA2", tooltips[10]),
            radioTooltip("coop_subset", "RNA3", tooltips[11]),
            radioTooltip("coop_subset", "RNA4", tooltips[12]),
            radioTooltip("coop_subset", "RNA5", tooltips[13])
          )
        )),

        fluidRow(
          ### Co-interaction statistics plot
          column(6, h3("Top TFs co-interacting with TF A"),
            plotlyOutput("coop_top", height = "40vh", width = "90%") %>%
              withSpinner(type = 5, color = "#0dc5c1")
          ),

          ### Gene expression plot comparing TF A and B
          column(6, h3("TF A + TF B comparison"),
            plotOutput("coop_exprs", height = "40vh", width = "90%") %>%
              withSpinner(type = 5, color = "#0dc5c1", hide.ui = FALSE)
          )
        ),
        fluidRow(
          ### Table of full co-interaction stats
          column(6, h3("Co-interaction stats (full table)"),
            DT::dataTableOutput("coop_stats_tbl", width = "90%") %>%
              withSpinner(type = 5, color = "#0dc5c1")
          ),

          ### Table of TF A & B targets
          column(6, h3("TF A & B target genes"),
            DT::dataTableOutput("coop_target_tbl", width = "90%") %>%
              withSpinner(type = 5, color = "#0dc5c1")
          )
        )
      )
    ), width = 10)
  )
)
