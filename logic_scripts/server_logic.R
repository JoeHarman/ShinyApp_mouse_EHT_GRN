function(input, output, session) {

  updateSelectizeInput(session, "gene",
    choices = unique(deg), server = TRUE)

  updateSelectizeInput(session, "enhancer",
    choices = unique(ATAC_stats$peak_coord), server = TRUE)

  updateSelectizeInput(session, "core_node",
    choices = tbl(sql_db, "GRN_nodes") %>% pull(id),
    selected = "Runx1", server = TRUE)

  updateSelectizeInput(session, "TF_A",
    choices = tbl(sql_db, "GRN_nodes") %>% filter(is_TF == 1) %>% pull(id),
    selected = "Runx1", server = TRUE)

  updateSelectizeInput(session, "TF_B",
    choices = tbl(sql_db, "GRN_nodes") %>% filter(is_TF == 1) %>% pull(id),
    selected = "Ikzf1", server = TRUE)


  ### Action button event handling
  ### On button click, subset samples and calculate PCA
  expression_data <- reactive({

    # Set selected groups to variable. Required step, as
    # SQL database queries don't recognise input variables.
    selected_groups <- unlist(input$select_groups)

    ### Init list
    exprs_list <- list()

    ### Process RNA expression data
    exprs_list$RNA_exprs <- tbl(sql_db, "RNA_exprs") %>%
      dplyr::filter(Group %in% selected_groups)

    ### Process ATAC expression data
    exprs_list$ATAC_exprs <- tbl(sql_db, "ATAC_exprs") %>%
      dplyr::filter(Group %in% selected_groups)

    ### Return data
    return(exprs_list)

  })

  run_pca_rna <- eventReactive(input$runPCA_RNA, {

    # Set selected groups to variable. Required step, as
    # SQL database queries don't recognise input variables.
    selected_groups <- unlist(input$select_groups)

    if (identical(selected_groups, samples[c(1, 2, 4, 5, 8, 9)])) {
      return(pca_tables$RNA)
    }

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

    rna_pca <- prcomp(t(vsd))
    rm(vsd)
    rna_pca$x <- data.frame(rna_pca$x[, 1:2]) %>%
      rownames_to_column("Sample") %>%
      left_join(RNA_anno)
    rna_pca$rotation <- data.frame(rna_pca$rotation[, 1:2]) %>%
      rownames_to_column("GeneID")

    return(rna_pca)

  }, ignoreNULL = FALSE) # Allows processing on app startup

  run_pca_atac <- eventReactive(input$runPCA_ATAC, {

    # Set selected groups to variable. Required step, as
    # SQL database queries don't recognise input variables.
    selected_groups <- unlist(input$select_groups)

    if (identical(selected_groups, samples[c(1, 2, 4, 5, 8, 9)])) {
      return(pca_tables$ATAC)
    }

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

    atac_pca <- prcomp(t(vsd))
    rm(vsd)
    atac_pca$x <- data.frame(atac_pca$x[, 1:2]) %>%
      rownames_to_column("Sample") %>%
      left_join(ATAC_anno)
    atac_pca$rotation <- data.frame(atac_pca$rotation[, 1:2]) %>%
      rownames_to_column("peak_coord")

    return(atac_pca)

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

    expression_data()$RNA_exprs %>%
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
        xlab("") +
        scale_fill_manual(
          values = col_scheme[[as.numeric(input$RNA_exprs_anno)]]) +
        theme_classic() +
        theme(text = element_text(size = font_size),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  })

  # ATAC expression plot
  output$ATAC_exprs <- renderPlot({
    enhancer_select <- input$enhancer

    expression_data()$ATAC_exprs %>%
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
        xlab("") +
        scale_fill_manual(
          values = col_scheme[[as.numeric(input$ATAC_exprs_anno)]]) +
        theme_classic() +
        theme(text = element_text(size = font_size),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  })

  # RNA PCA plot
  output$RNA_pca <- renderPlot({
    pca_res <- run_pca_rna()
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
      theme_classic() +
      theme(text = element_text(size = font_size))
  })

  # ATAC PCA plot
  output$ATAC_pca <- renderPlot({
    pca_res <- run_pca_atac()
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
      theme_classic() +
      theme(text = element_text(size = font_size))
  })

  # NETWORK PLOT

  # Process network button
  grn_list <- eventReactive(input$makeGRN, {

    # Set selected groups to variable. Required step, as
    # SQL database queries don't recognise input variables.
    grn_subset <- unlist(input$grn_subset)
    grn_centrality <- unlist(input$grn_centrality)
    core_node <- unlist(input$core_node)

    ### Init list
    grn_list <- list()
    node_filt <- "TF"

    # Init SQL tables
    edges <- tbl(sql_db, "GRN_edges")
    nodes <- tbl(sql_db, "GRN_nodes")
    centrality <- tbl(sql_db, "GRN_centrality")

    # Filter for TFs
    if (input$grn_tf) {
      nodes <- filter(nodes, is_TF == 1)
      edges <- filter(edges,
        to %in% !!pull(nodes, id) & from %in% !!pull(nodes, id))
      centrality <- filter(centrality, is_TF == 1)
    }

    # Filter centrality stats
    # Note: Need to collect early due to dbplyr bug
    # (see https://github.com/tidyverse/dbplyr/issues/1206)
    col_select <- paste(grn_subset, grn_centrality, sep = "_")

    centrality <- centrality %>%
      filter(ifelse(grn_subset == "Full", TRUE, RNA_module == grn_subset)) %>%
      select_if(grepl(paste0(grn_subset, "|name|RNA_module"), names(.)))

    # Process network modes
    if (input$grn_mode == "Central TFs") {

      # Filter for subset (RNA1-5 or full)
      nodes <- nodes %>%
        filter(ifelse(grn_subset == "Full", TRUE, RNA_module == grn_subset))
      edges <- filter(edges,
        to %in% !!pull(nodes, id) & from %in% !!pull(nodes, id))

      top_n_genes <- centrality %>%
        arrange(desc(.data[[col_select]])) %>%
        head(input$top_n_centrality) %>%
        collect() %>%
        pull(name)

      nodes <- filter(nodes, id %in% top_n_genes) %>%
        left_join(centrality, by = c("id" = "name", "RNA_module"))
      edges <- filter(edges, to %in% top_n_genes & from %in% top_n_genes)

    } else if (input$grn_mode == "Upstream") {

      # Filter for subset (RNA1-5 or full)
      nodes <- nodes %>%
        filter(ifelse(grn_subset == "Full",
          TRUE,
          RNA_module == grn_subset | id == core_node))
      edges <- filter(edges,
        to %in% !!pull(nodes, id) & from %in% !!pull(nodes, id))

      # Process upstream network
      edges <- filter(edges, to == core_node)
      nodes <- filter(nodes,
        id %in% !!pull(edges, from) | id %in% !!pull(edges, to)) %>%
        left_join(centrality, by = c("id" = "name", "RNA_module"))

    } else if (input$grn_mode == "Downstream") {

      # Filter for subset (RNA1-5 or full)
      nodes <- nodes %>%
        filter(ifelse(grn_subset == "Full",
          TRUE,
          RNA_module == grn_subset | id == core_node))
      edges <- filter(edges,
        to %in% !!pull(nodes, id) & from %in% !!pull(nodes, id))

      # Process downstream network
      edges <- filter(edges, from == core_node)
      nodes <- filter(nodes,
        id %in% !!pull(edges, from) | id %in% !!pull(edges, to)) %>%
        left_join(centrality, by = c("id" = "name", "RNA_module"))
    }

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

    grn_list$edges <- e %>%
      left_join(collect(edges))
    grn_list$nodes <- n
    grn_list$col_select <- col_select
    grn_list$plot_mode <- input$grn_mode


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

  output$downGRN <- downloadHandler(
    filename = function() {
      "grn_interactions.csv"
    },
    content = function(fname) {
      write.csv(grn_list()$edges, fname)
    }
  )

  output$networkSidePlot <- renderPlot({

    n <- grn_list()$nodes
    e <- grn_list()$edges
    col_select <- grn_list()$col_select
    plot_mode <- grn_list()$plot_mode

    if (plot_mode == "Central TFs") {
      net_topn <- n %>%
        select(id:RNA_module, Centrality = col_select) %>%
        top_n(20, Centrality) %>%
        arrange(Centrality) %>%
        mutate(id = factor(id, levels = id))

      net_topn_plot <- ggplot(net_topn,
        aes(x = Centrality, y = id, fill = RNA_module)) +
        xlab(gsub("_", " ", col_select))

    } else if (plot_mode == "Upstream") {
      net_topn <- e %>%
        select(id = from, RNA_correlation) %>%
        left_join(n) %>%
        top_n(if (input$corr_filt == "Top-20 correlated") {
            20
          } else {
            -20
          }, RNA_correlation) %>%
        arrange(RNA_correlation) %>%
        mutate(id = factor(id, levels = (id)))

        net_topn_plot <- ggplot(net_topn,
          aes(x = RNA_correlation, y = id, fill = RNA_module)) +
          xlab("RNA correlation")

    } else {
      net_topn <- e %>%
        select(id = to, RNA_correlation) %>%
        left_join(n) %>%
        top_n(if (input$corr_filt == "Top-20 correlated") {
            20
          } else {
            -20
          }, RNA_correlation) %>%
        arrange(RNA_correlation) %>%
        mutate(id = factor(id, levels = (id)))

        net_topn_plot <- ggplot(net_topn,
          aes(x = RNA_correlation, y = id, fill = RNA_module)) +
          xlab("RNA correlation")

    }

    net_topn_plot <- net_topn_plot +
      geom_point(pch = 21, size = 3) +
      scale_fill_manual(values = c(RNA1 = "yellow", RNA2 = "blue",
        RNA3 = "green", RNA4 = "orange", RNA5 = "red")) +
      ylab("") +
      theme_bw() +
      theme(text = element_text(size = font_size), legend.position = "bottom",
        legend.direction = "vertical",
        legend.text = element_text(size = font_size - 2)) +
      guides(fill = guide_legend(ncol = 2))

    return(net_topn_plot)

  })

  # Cooperation expression plot
  output$coop_exprs <- renderPlot({
    tf_a <- input$TF_A
    tf_b <- input$TF_B
    exprs <- expression_data()$RNA_exprs %>%
      dplyr::filter(GeneID %in% c(tf_a, tf_b)) %>%
      collect() %>%
      mutate(Group = factor(Group, levels = unique(Group)))

      ggplot(exprs, aes(x = Group, y = CPM, fill = GeneID)) +
        stat_summary(fun.y = mean, geom = "bar", col = "black",
          position = position_dodge()) +
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
          position = position_dodge(width = 0.9)) +
        geom_point(position = position_jitterdodge(
          jitter.width = 0.2, dodge.width = 0.9, seed = 12345)) +
        ggtitle(paste0(tf_a, " and ", tf_b)) +
        xlab("") +
        theme_classic() +
        theme(text = element_text(size = font_size),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  })

  # Cooperation expression plot
  output$coop_top <- renderPlotly({
    tf_a <- input$TF_A
    tf_b <- input$TF_B
    module <- input$coop_subset

    if (module != "All") {
      gene_filt <- tbl(sql_db, "GRN_nodes") %>%
        filter(RNA_module == module) %>%
        pull(id)
    } else {
      gene_filt <- tbl(sql_db, "GRN_nodes") %>%
        pull(id)
    }

    # TF A selection must always be in filter
    gene_filt <- unique(c(tf_a, gene_filt))

    df <- coop_stats %>%
      filter(TF_A %in% gene_filt & TF_B %in% gene_filt) %>%
      filter(TF_A == tf_a | TF_B == tf_a) %>%
      # Need to ensure label is the partner TF, not selected TF
      mutate(GeneID = paste0(TF_A, TF_B)) %>%
      mutate(GeneID = gsub(tf_a, "", GeneID))

    p <- ggplot(df, aes(y = -log10(Cointeraction_Pvalue), x = RNA_correlation,
      label = GeneID)) +
      geom_point(pch = 21, fill = "#996cd4") +
      # Filter for partner TF_B and highlight in red
      geom_point(data = filter(df, TF_A == tf_b | TF_B == tf_b),
        pch = 21, fill = "red", size = 3) +
      xlab(paste0("RNA correlation with ", tf_a)) +
      ylab(paste("-log10 P-value for\nco-interaction with", tf_a)) +
      theme_bw() +
      theme(text = element_text(size = font_size),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

    ggplotly(p) %>%
      layout(
        xaxis = list(
          title = list(font = list(size = font_size)),
          tickfont = list(size = font_size - 2)),
        yaxis = list(
          title = list(font = list(size = font_size)),
          tickfont = list(size = font_size - 2)))

  })

  # Cointeraction stats table
  output$coop_stats_tbl <- DT::renderDataTable(
    coop_stats,
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

  # Cointeraction stats table

  coop_target_tbl <- reactive({

    tf_a <- input$TF_A
    tf_b <- input$TF_B
    nodes <- tbl(sql_db, "GRN_nodes") %>%
      select(GeneID = id, RNA_module) %>%
      collect()

    tf_a_targets <- tbl(sql_db, "GRN_edges") %>%
      filter(from == tf_a) %>%
      pull(to)
    tf_b_targets <- tbl(sql_db, "GRN_edges") %>%
      filter(from == tf_b) %>%
      pull(to)

    joint_targets <- intersect(tf_a_targets, tf_b_targets)

    tf_a_target_tbl <- tbl(sql_db, "GRN_edges") %>%
      filter(to %in% joint_targets & from == tf_a) %>%
      select(to, RNA_correlation) %>%
      collect() %>%
      purrr::set_names(c("GeneID", paste0("Correlation with ", tf_a)))

    tf_b_target_tbl <- tbl(sql_db, "GRN_edges") %>%
      filter(to %in% joint_targets & from == tf_b) %>%
      select(to, RNA_correlation) %>%
      collect() %>%
      purrr::set_names(c("GeneID", paste0("Correlation with ", tf_b)))

    inner_join(tf_a_target_tbl, tf_b_target_tbl) %>%
      arrange(GeneID) %>%
      left_join(nodes)

  })

  output$coop_target_tbl <- DT::renderDataTable(
    coop_target_tbl(),
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



}
