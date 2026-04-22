# ─── Tab Visibility Management ─────────────────────────────────────────────────

# Hide Correlation Analysis tab on initial load
observe({
    hideTab(inputId = "NAVTABS", target = "Correlation Analysis")
})

# Show after custom data submit
observeEvent(input$submit, {
    mat <- MatrixReactive()
    if (!is.null(mat) && nrow(mat) > 0) {
        showTab(inputId = "NAVTABS", target = "Correlation Analysis")
        updateTabsetPanel(session, inputId = "NAVTABS", selected = "Correlation Analysis")
        shinyjs::delay(300, shinyjs::runjs("$(window).trigger('resize');"))
    }
})

# Show after demo data submit
observeEvent(input$demo_submit, {
    showTab(inputId = "NAVTABS", target = "Correlation Analysis")
    updateTabsetPanel(session, inputId = "NAVTABS", selected = "Correlation Analysis")
    shinyjs::delay(300, shinyjs::runjs("$(window).trigger('resize');"))
})

# ─── Pair Viewer Sub-tab Visibility (Mode A only) ────────────────────────────

observe({
    mode <- AnalysisMode()
    if (mode == "single") {
        showTab(inputId = "CorrelationTabs", target = "Pair Viewer")
    } else {
        hideTab(inputId = "CorrelationTabs", target = "Pair Viewer")
    }
})

# ─── Populate Target Gene Selector ───────────────────────────────────────────

observe({
    mat <- ProcessedMatrix()
    req(mat)
    genes <- rownames(mat)
    updateSelectizeInput(session, "target_gene", choices = genes, server = TRUE,
        options = list(placeholder = 'Type gene name...'))
})
