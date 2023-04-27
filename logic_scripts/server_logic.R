function(input, output, session) {

  # Check_credentials returns a function to authenticate users
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )

  updateSelectizeInput(session, "enhancer",
    choices = unique(dae), server = TRUE)

  data <- eventReactive(input$subsetSamples, {

    ### Init list
    processedData <- list()

    ### Process RNA expression data
    processedData$RNA_exprs <- filter(
      RNA_exprs, Group %in% unlist(input$select_groups))

    ### Process ATAC expression data
    processedData$ATAC_exprs <- filter(
      ATAC_exprs, Group %in% unlist(input$select_groups))

    ### Calculate RNA PCA

    RNA_exprs_wide_filt <- RNA_exprs_wide[,
      grepl(paste(unlist(input$select_groups), collapse = "|"),
      colnames(RNA_exprs_wide))]

    vsd <- varianceStabilizingTransformation(round(
      as.matrix(RNA_exprs_wide_filt)))
    rv <- rowVars(vsd)
    select <- order(rv, decreasing = TRUE)[
      seq_len(min(input$RNA_topn, length(rv)))]

    RNA_pca <- prcomp(t(vsd[select, ]))
    RNA_pca$x <- data.frame(RNA_pca$x[, 1:2]) %>%
      rownames_to_column("Sample") %>%
      left_join(RNA_anno)
    RNA_pca$rotation <- data.frame(RNA_pca$rotation[, 1:2]) %>%
      rownames_to_column("GeneID")

    processedData$RNA_pca <- RNA_pca

    ### Calculate ATAC PCA

    ATAC_exprs_wide_filt <- ATAC_exprs_wide[,
      grepl(paste(unlist(input$select_groups), collapse = "|"),
      colnames(ATAC_exprs_wide))]

    vsd <- varianceStabilizingTransformation(round(
      as.matrix(na.omit(ATAC_exprs_wide_filt))))
    rv <- rowVars(vsd)
    select <- order(rv, decreasing = TRUE)[
      seq_len(min(input$ATAC_topn, length(rv)))]

    ATAC_pca <- prcomp(t(vsd[select, ]))
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
    server = FALSE,
    options = list(
      paging = TRUE,    ## paginate the output
      pageLength = 10,  ## number of rows to output for each page
      scrollX = TRUE,   ## enable scrolling on X axis
      scrollY = TRUE,   ## enable scrolling on Y axis
      autoWidth = TRUE, ## use smart column width handling
      server = FALSE,   ## use client-side processing
      dom = "Bfrtip",
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
   # , '#c1d9da', '#fff'
  )

  # ATAC data table
  output$ATAC_stats_tbl <- DT::renderDataTable(
    ATAC_stats,
    server = FALSE,
    options = list(
      paging = TRUE,    ## paginate the output
      pageLength = 10,  ## number of rows to output for each page
      scrollX = TRUE,   ## enable scrolling on X axis
      scrollY = TRUE,   ## enable scrolling on Y axis
      autoWidth = TRUE, ## use smart column width handling
      server = FALSE,   ## use client-side processing
      dom = "Bfrtip",
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
    #filter(RNA_exprs, GeneID == input$gene) %>%
    filter(data()$RNA_exprs, GeneID == input$gene) %>%
    #filter(Group %in% unlist(input$select_groups)) %>%
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
    #filter(ATAC_exprs, peak_coord == input$enhancer) %>%
    filter(data()$ATAC_exprs, peak_coord == input$enhancer) %>%
    #filter(Group %in% unlist(input$select_groups)) %>%
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
    pca_res$rotation <- filter(pca_res$rotation, GeneID == input$gene)

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
    pca_res$rotation <- filter(pca_res$rotation, peak_coord == input$enhancer)

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

}