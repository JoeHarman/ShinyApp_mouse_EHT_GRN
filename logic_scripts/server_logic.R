function(input, output, session) {

  # Check_credentials returns a function to authenticate users
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )

    updateSelectizeInput(session, "enhancer",
      choices = unique(ATAC_exprs$peak_coord), server = TRUE)

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
    filter(RNA_exprs, GeneID == input$gene) %>%
    filter(Group %in% unlist(input$select_groups)) %>%
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
    filter(ATAC_exprs, peak_coord == input$enhancer) %>%
    filter(Group %in% unlist(input$select_groups)) %>%
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

}
