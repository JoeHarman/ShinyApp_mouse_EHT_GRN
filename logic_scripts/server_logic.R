function(input, output, session){

  # Check_credentials returns a function to authenticate users
  res_auth <- secure_server(
    check_credentials = check_credentials(credentials)
  )

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
      dom = 'Bfrtip',
      buttons = c('csv', 'excel'),
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
      dom = 'Bfrtip',
      buttons = c('csv', 'excel'),
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
