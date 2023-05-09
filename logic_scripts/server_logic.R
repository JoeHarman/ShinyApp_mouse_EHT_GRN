function(input, output, session) {

  # Check_credentials returns a function to authenticate users
  # Note - this should be removed in final release
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )

  updateSelectizeInput(session, "gene",
    choices = unique(deg), server = TRUE)

  updateSelectizeInput(session, "enhancer",
    choices = unique(ATAC_stats$peak_coord), server = TRUE)

  ### Action button event handling
  ### On button click, subset samples and calculate PCA
  data <- eventReactive(input$subsetSamples, {

    # Set selected groups to variable. Required step, as
    # SQL database queries don't recognise input variables.
    selected_groups <- unlist(input$select_groups)

    ### Init list
    processedData <- list()

    ### Process RNA expression data
    processedData$RNA_exprs <- tbl(sql_db, "RNA_exprs") %>%
      dplyr::filter(Group %in% selected_groups)

    ### Process ATAC expression data
    processedData$ATAC_exprs <- tbl(sql_db, "ATAC_exprs") %>%
      dplyr::filter(Group %in% selected_groups)

    ### Calculate RNA PCA
    select_ind <- c(1, grep(
      paste(selected_groups, collapse = "|"),
      colnames(tbl(sql_db, "RNA_exprs_wide"))))

    vsd <- tbl(sql_db, "RNA_exprs_wide") %>%
      dplyr::select(c(1, select_ind)) %>%
      collect() %>%
      column_to_rownames("GeneID") %>%
      as.matrix() %>%
      round() %>%
      varianceStabilizingTransformation()

    rv <- rowVars(vsd)
    select <- order(rv, decreasing = TRUE)[
      seq_len(min(input$RNA_topn, length(rv)))]
    vsd <- vsd[select, ]

    RNA_pca <- prcomp(t(vsd))
    rm(vsd)
    RNA_pca$x <- data.frame(RNA_pca$x[, 1:2]) %>%
      rownames_to_column("Sample") %>%
      left_join(RNA_anno)
    RNA_pca$rotation <- data.frame(RNA_pca$rotation[, 1:2]) %>%
      rownames_to_column("GeneID")

    processedData$RNA_pca <- RNA_pca

    ### Calculate ATAC PCA
    select_ind <- c(1, grep(
      paste(selected_groups, collapse = "|"),
      colnames(tbl(sql_db, "ATAC_exprs_wide"))))

    vsd <- tbl(sql_db, "ATAC_exprs_wide") %>%
      dplyr::select(c(1, select_ind)) %>%
      collect() %>%
      column_to_rownames("peak_coord") %>%
      na.omit() %>%
      as.matrix() %>%
      round() %>%
      varianceStabilizingTransformation()

    rv <- rowVars(vsd)
    select <- order(rv, decreasing = TRUE)[
      seq_len(min(input$ATAC_topn, length(rv)))]
    vsd <- vsd[select, ]

    ATAC_pca <- prcomp(t(vsd))
    rm(vsd)
    ATAC_pca$x <- data.frame(ATAC_pca$x[, 1:2]) %>%
      rownames_to_column("Sample") %>%
      left_join(ATAC_anno)
    ATAC_pca$rotation <- data.frame(ATAC_pca$rotation[, 1:2]) %>%
      rownames_to_column("peak_coord")

    processedData$ATAC_pca <- ATAC_pca

    ### Return data
    return(processedData)

  }, ignoreNULL = FALSE) # Allows processing on app startup

  # RNA data table
  output$RNA_stats_tbl <- DT::renderDataTable(
    RNA_stats,
    options = list(
      paging = TRUE,    ## paginate the output
      pageLength = 10,  ## number of rows to output for each page
      lengthMenu = c(10, 50, 100),
      scrollX = TRUE,   ## enable scrolling on X axis
      scrollY = TRUE,   ## enable scrolling on Y axis
      autoWidth = TRUE, ## use smart column width handling
      server = TRUE,   ## use client-side processing
      dom = "lfrtiBp", # Alternative: dom = "Bfrtlip"
      buttons = c("csv", "excel"),
      formatter = list(
        `E8-EC_E8-preHE_FDR` = p_formatter,
        `E8-preHE_E9-HE_FDR` = p_formatter,
        `E9-HE_E9-proHSC_FDR` = p_formatter,
        `E9-proHSC_E10-preI_FDR` = p_formatter,
        `E10-preI_E10-preII_FDR` = p_formatter
      )
    ),
    extensions = "Buttons",
    selection = "single",
    filter = "top",
    rownames = FALSE
  )

  # ATAC data table
  output$ATAC_stats_tbl <- DT::renderDataTable(
    ATAC_stats,
    server = FALSE,
    options = list(
      paging = TRUE,    ## paginate the output
      pageLength = 10,  ## number of rows to output for each page
      lengthMenu = c(10, 50, 100),
      scrollX = TRUE,   ## enable scrolling on X axis
      scrollY = TRUE,   ## enable scrolling on Y axis
      autoWidth = TRUE, ## use smart column width handling
      server = TRUE,   ## use client-side processing
      dom = "lfrtiBp", # Alternative: dom = "Bfrtlip"
      buttons = c("csv", "excel"),
      formatter = list(
        `E8-EC_E8-preHE_FDR` = p_formatter,
        `E8-preHE_E9-HE_FDR` = p_formatter,
        `E9-HE_E9-proHSC_FDR` = p_formatter,
        `E9-proHSC_E10-preI_FDR` = p_formatter,
        `E10-preI_E10-preII_FDR` = p_formatter
      )
    ),
    extensions = "Buttons",
    selection = "single",
    filter = "top",
    rownames = FALSE
  )

  # RNA expression plot
  output$RNA_exprs <- renderPlot({
    gene_select <- input$gene

    data()$RNA_exprs %>%
      dplyr::filter(GeneID == gene_select) %>%
      collect() %>%
      mutate(Group = factor(Group, levels = unique(Group))) %>%
      ggplot(aes(x = Group, y = CPM)) +
        {if (input$RNA_exprs_anno == "1") {
          stat_summary(
            aes(fill = Population), fun.y = mean, geom = "bar", col = "black")
        }else if (input$RNA_exprs_anno == "2") {
          stat_summary(
            aes(fill = Stage), fun.y = mean, geom = "bar", col = "black")
        }else if (input$RNA_exprs_anno == "3") {
          stat_summary(
            aes(fill = Group), fun.y = mean, geom = "bar", col = "black")
        }} +
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
        geom_point(position = position_jitter(width = 0.2, seed = 12345)) +
        ggtitle(input$gene) +
        scale_fill_manual(
          values = col_scheme[[as.numeric(input$RNA_exprs_anno)]]) +
        theme_classic()
  })

  # ATAC expression plot
  output$ATAC_exprs <- renderPlot({
    enhancer_select <- input$enhancer

    data()$ATAC_exprs %>%
      dplyr::filter(peak_coord == enhancer_select) %>%
      collect() %>%
      mutate(Group = factor(Group, levels = unique(Group))) %>%
      ggplot(aes(x = Group, y = CPM)) +
        {if (input$ATAC_exprs_anno == "1") {
          stat_summary(
            aes(fill = Population), fun.y = mean, geom = "bar", col = "black")
        }else if (input$ATAC_exprs_anno == "2") {
          stat_summary(
            aes(fill = Stage), fun.y = mean, geom = "bar", col = "black")
        }else if (input$ATAC_exprs_anno == "3") {
          stat_summary(
            aes(fill = Group), fun.y = mean, geom = "bar", col = "black")
        }} +
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
        geom_point(position = position_jitter(width = 0.2, seed = 12345)) +
        ggtitle(input$enhancer) +
        scale_fill_manual(
          values = col_scheme[[as.numeric(input$ATAC_exprs_anno)]]) +
        theme_classic()
  })

  # RNA PCA plot
  output$RNA_pca <- renderPlot({
    pca_res <- data()$RNA_pca
    pca_res$rotation <- dplyr::filter(pca_res$rotation, GeneID == input$gene)

    ggplot(pca_res$x, aes(x = PC1, y = PC2)) +
      {if (input$RNA_exprs_anno == "1") {
        geom_point(aes(col = Population), size = 3)
      }else if (input$RNA_exprs_anno == "2") {
        geom_point(aes(col = Stage), size = 3)
      }else if (input$RNA_exprs_anno == "3") {
        geom_point(aes(col = Group), size = 3)
      }} +
      geom_point(
        data = pca_res$rotation,
        mapping = aes(x = PC1 * 1000, y = PC2 * 1000),
          col = "red", size = 3) +
      geom_text(data = pca_res$rotation,
        mapping = aes(x = PC1 * 1000, y = PC2 * 1000, label = GeneID),
        col = "red") +
      geom_segment(data = pca_res$rotation,
        mapping = aes(x = 0, y = 0, xend = PC1 * 1000 * 0.85,
        yend = PC2 * 1000 * 0.85), arrow = arrow(length = unit(0.5, "cm"))) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_colour_manual(
        values = col_scheme[[as.numeric(input$RNA_exprs_anno)]]) +
      theme_classic()
  })

  # ATAC PCA plot
  output$ATAC_pca <- renderPlot({
    pca_res <- data()$ATAC_pca
    pca_res$rotation <- dplyr::filter(
      pca_res$rotation, peak_coord == input$enhancer)

    ggplot(pca_res$x, aes(x = PC1, y = PC2)) +
      {if (input$ATAC_exprs_anno == "1") {
        geom_point(aes(col = Population), size = 3)
      }else if (input$ATAC_exprs_anno == "2") {
        geom_point(aes(col = Stage), size = 3)
      }else if (input$ATAC_exprs_anno == "3") {
        geom_point(aes(col = Group), size = 3)
      }} +
      geom_point(
        data = pca_res$rotation,
        mapping = aes(x = PC1 * 2000, y = PC2 * 2000),
          col = "red", size = 3) +
      geom_text(data = pca_res$rotation,
        mapping = aes(x = PC1 * 2000, y = PC2 * 2000, label = peak_coord),
        col = "red") +
      geom_segment(data = pca_res$rotation,
        mapping = aes(x = 0, y = 0, xend = PC1 * 2000 * 0.85,
        yend = PC2 * 2000 * 0.85), arrow = arrow(length = unit(0.5, "cm"))) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      scale_colour_manual(
        values = col_scheme[[as.numeric(input$ATAC_exprs_anno)]]) +
      theme_classic()
  })

  # NETWORK PLOT
  # Need a data grabber function and plot function
  # Grabber handles SQL queries (maybe loads full nodes, subsets edges)...
  # ...then possibly initialises visNetwork
  # Plot function handles annotation adjustments, i.e, colour options.
  # Maybe best for visNetwork to initialise here then... See how it works.

  # Process network button
  grn_list <- eventReactive(input$makeGRN, {

    # Set selected groups to variable. Required step, as
    # SQL database queries don't recognise input variables.
    grn_subset <- unlist(input$grn_subset)
    grn_centrality <- unlist(input$grn_centrality)

    ### Init list
    grn_list <- list()
    node_filt <- "TF"

    # Init SQL tables
    edges <- tbl(sql_db, "GRN_edges")
    nodes <- tbl(sql_db, "GRN_nodes")
    centrality <- tbl(sql_db, "GRN_centrality")

    # Filter for TFs
    if (node_filt == "TF") {
      nodes <- filter(nodes, is_TF == 1)
      edges <- filter(edges,
        to %in% !!pull(nodes, id) & from %in% !!pull(nodes, id))
      centrality <- filter(centrality, is_TF == 1)
    }

    # Attach centrality stats and filter
    # Note: Need to collect early due to dbplyr bug
    # (see https://github.com/tidyverse/dbplyr/issues/1206)
    col_select <- paste(grn_subset, grn_centrality, sep = "_")

    centrality <- centrality %>%
      filter(ifelse(grn_subset == "Full", TRUE, RNA_module == grn_subset)) %>%
      select_if(grepl(paste0(grn_subset, "|name|RNA_module"), names(.)))

    top_n_genes <- centrality %>%
      arrange(desc(.data[[col_select]])) %>%
      head(input$top_n_centrality) %>%
      collect() %>%
      pull(name)

    nodes <- filter(nodes, id %in% top_n_genes) %>%
      left_join(centrality, by = c("id" = "name"))
    edges <- filter(edges, to %in% top_n_genes & from %in% top_n_genes)

    # Ensure there are no isolated nodes, and remove duplicate/self edges
    net <- igraph::graph_from_data_frame(
      collect(edges), directed = TRUE, vertices = collect(nodes)) %>%
      remove_isolated() %>%
      igraph::simplify()

    # Prepare visNetwork tables and export
    n <- igraph::as_data_frame(net, what = "vertices") %>%
      mutate(label = name, value = .data[[col_select]])
    colnames(n)[1] <- "id"
    e <- igraph::as_data_frame(net, what = "edges")
    colnames(e)[1:2] <- c("from", "to")

    grn_list$edges <- e
    grn_list$nodes <- n

    ### Return data
    return(grn_list)

  }, ignoreNULL = FALSE) # Allows processing on app startup

  output$mynetworkid <- renderVisNetwork({

    # Render visNetwork

    visNetwork(grn_list()$nodes, grn_list()$edges) %>%
      # Node specification
      visNodes(
        font = list(size = 30),
        borderWidth = 2,
        borderWidthSelected = 2,
        color = list(
          background = "#b8b3e1",
          border = "#6e64bd",
          highlight = list(background = "#d65a5a",
            border = "black")
        )) %>%

      # Edge specification
      visEdges(physics = FALSE,
        arrows = list(to = list(enabled = TRUE, scaleFactor = 1)),
        color = list(highlight = "#b72929")) %>%

      # Highlight options
      visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1),
        nodesIdSelection = TRUE) %>%

      # Layout options
      visIgraphLayout(smooth = FALSE, randomSeed = 12345) %>%
      visInteraction(navigationButtons = TRUE) %>%

      # Save image option
      visExport("png", float = "left")
      # Ideally would include legend

  })

}