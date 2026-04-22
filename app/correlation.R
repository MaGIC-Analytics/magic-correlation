# ═══════════════════════════════════════════════════════════════════════════════
# correlation.R — Mode A (Single Gene vs All) + Mode B (Gene List Matrix)
# ═══════════════════════════════════════════════════════════════════════════════

# ─── MODE A: Single Gene vs All ──────────────────────────────────────────────

SingleCorResult <- eventReactive(input$run_single, {
    mat <- ProcessedMatrix()
    req(mat, input$target_gene)

    shiny::validate(need(input$target_gene %in% rownames(mat),
        "Selected gene not found in the expression matrix."))

    withProgress(message = "Computing correlations...", value = 0.3, {

        target_vec <- as.numeric(mat[input$target_gene, ])
        n <- length(target_vec)
        other_genes <- setdiff(rownames(mat), input$target_gene)

        # Apply min expression filter
        if (input$min_expr_filter > 0) {
            gene_means <- rowMeans(mat[other_genes, , drop=FALSE], na.rm=TRUE)
            other_genes <- other_genes[gene_means >= input$min_expr_filter]
        }

        shiny::validate(need(length(other_genes) >= 1, "No genes remain after filtering."))

        # Vectorized correlation
        mat_other <- mat[other_genes, , drop=FALSE]
        r_vals <- as.numeric(cor(t(mat_other), target_vec, method = input$cor_method, use = "pairwise.complete.obs"))

        setProgress(0.6)

        # Analytical p-values
        t_stat <- r_vals * sqrt((n - 2) / (1 - r_vals^2))
        t_stat[is.nan(t_stat)] <- 0
        p_vals <- 2 * pt(-abs(t_stat), df = n - 2)
        padj   <- p.adjust(p_vals, method = "BH")

        setProgress(0.8)

        result <- data.frame(
            Gene        = other_genes,
            Correlation = round(r_vals, 4),
            P_value     = signif(p_vals, 4),
            Adj_P_value = signif(padj, 4),
            Abs_Cor     = round(abs(r_vals), 4),
            stringsAsFactors = FALSE
        )

        # Sort by absolute correlation descending
        result <- result[order(-result$Abs_Cor), ]
        rownames(result) <- NULL

        setProgress(1.0)
    })

    result
})

# Full table (before filtering by thresholds) for download
SingleCorFull <- reactive({
    SingleCorResult()
})

# Filtered table for display
SingleCorFiltered <- reactive({
    df <- SingleCorResult()
    req(df)
    df <- df[df$Abs_Cor >= input$min_abs_cor, ]
    df <- df[df$Adj_P_value <= input$pval_cutoff, ]
    df
})

# ─── Mode A: Scatter Grid Plotter ────────────────────────────────────────────

scatter_grid_plotter <- reactive({
    df <- SingleCorFiltered()
    req(df)
    mat <- ProcessedMatrix()
    req(mat, input$target_gene)

    top_genes <- head(df$Gene, input$top_n)
    shiny::validate(need(length(top_genes) >= 1, "No genes pass the current filters."))

    target_vec <- as.numeric(mat[input$target_gene, ])

    # Build long-format data for faceting
    plot_data <- do.call(rbind, lapply(top_genes, function(g) {
        data.frame(
            Target  = target_vec,
            Partner = as.numeric(mat[g, ]),
            Gene    = g,
            stringsAsFactors = FALSE
        )
    }))

    # Add correlation label per gene
    cor_labels <- df[df$Gene %in% top_genes, c("Gene", "Correlation", "Adj_P_value")]
    plot_data$Gene <- factor(plot_data$Gene, levels = top_genes)

    pt_size  <- input$point_size  %||% 2
    pt_alpha <- input$point_alpha %||% 0.6
    ax_size  <- input$axis_size   %||% 8
    ttl_size <- input$title_size  %||% 14
    tk_size  <- input$tick_size   %||% 7

    p <- ggplot(plot_data, aes(x = Target, y = Partner)) +
        geom_point(size = pt_size, alpha = pt_alpha, color = "#377EB8") +
        geom_smooth(method = "lm", se = TRUE, color = "#d73027", linewidth = 0.8) +
        stat_cor(method = input$cor_method, size = tk_size * 0.35, color = "#333333") +
        facet_wrap(~ Gene, scales = "free_y", ncol = min(4, ceiling(sqrt(length(top_genes))))) +
        labs(
            x = paste(input$target_gene, "expression"),
            y = "Partner gene expression",
            title = paste0("Top ", length(top_genes), " genes correlated with ", input$target_gene)
        ) +
        theme_bw() +
        theme(
            plot.title   = element_text(size = ttl_size, face = "bold", hjust = 0.5),
            axis.title   = element_text(size = ax_size),
            axis.text    = element_text(size = tk_size),
            strip.text   = element_text(size = ax_size, face = "bold")
        )

    if (input$legend_pos %||% "right" == "none") {
        p <- p + theme(legend.position = "none")
    }

    p
})

# ─── Mode A: Pair Viewer ─────────────────────────────────────────────────────

output$pair_viewer_ui <- renderUI({
    df <- SingleCorResult()
    if (is.null(df) || nrow(df) == 0) {
        return(p(em("Run a single gene correlation first to use the pair viewer."), style="color:#888;"))
    }

    meta <- ProcessedMeta()
    meta_choices <- if (!is.null(meta)) c("(none)" = "none", colnames(meta)[-1]) else c("(none)" = "none")

    tagList(
        fluidRow(
            column(4,
                selectizeInput("pair_gene", "Select Gene to Plot:",
                    choices = df$Gene, selected = df$Gene[1],
                    options = list(placeholder = 'Select a gene...'))
            ),
            column(4,
                selectInput("pair_color_by", "Color by Metadata:",
                    choices = meta_choices, selected = "none")
            ),
            column(4, style="padding-top:25px;",
                actionButton("plot_pair", "Plot Pair", class="btn btn-info", icon=icon("chart-line"))
            )
        ),
        hr()
    )
})

PairPlotData <- eventReactive(input$plot_pair, {
    mat  <- ProcessedMatrix()
    meta <- ProcessedMeta()
    req(mat, input$target_gene, input$pair_gene)

    target_vec  <- as.numeric(mat[input$target_gene, ])
    partner_vec <- as.numeric(mat[input$pair_gene, ])

    df <- data.frame(
        Target  = target_vec,
        Partner = partner_vec,
        stringsAsFactors = FALSE
    )

    # Add metadata color variable if selected
    if (!is.null(input$pair_color_by) && input$pair_color_by != "none" && !is.null(meta)) {
        df$ColorVar <- as.character(meta[[input$pair_color_by]])
    }

    df
})

pair_scatter_plotter <- reactive({
    df <- PairPlotData()
    req(df)

    pt_size  <- input$point_size  %||% 2
    pt_alpha <- input$point_alpha %||% 0.6
    ax_size  <- input$axis_size   %||% 8
    ttl_size <- input$title_size  %||% 14
    tk_size  <- input$tick_size   %||% 7

    has_color <- "ColorVar" %in% colnames(df)

    if (has_color) {
        p <- ggplot(df, aes(x = Target, y = Partner, color = ColorVar)) +
            geom_point(size = pt_size * 1.5, alpha = pt_alpha) +
            labs(color = input$pair_color_by)
    } else {
        p <- ggplot(df, aes(x = Target, y = Partner)) +
            geom_point(size = pt_size * 1.5, alpha = pt_alpha, color = "#377EB8")
    }

    p <- p +
        geom_smooth(method = "lm", se = TRUE, color = "#d73027", linewidth = 1) +
        stat_cor(method = input$cor_method, size = tk_size * 0.4, color = "#333333") +
        labs(
            x = paste(input$target_gene, "expression"),
            y = paste(input$pair_gene, "expression"),
            title = paste0(input$target_gene, " vs ", input$pair_gene)
        ) +
        theme_bw() +
        theme(
            plot.title   = element_text(size = ttl_size, face = "bold", hjust = 0.5),
            axis.title   = element_text(size = ax_size),
            axis.text    = element_text(size = tk_size),
            legend.text  = element_text(size = tk_size),
            legend.title = element_text(size = ax_size)
        )

    lp <- input$legend_pos %||% "right"
    if (lp == "none") {
        p <- p + theme(legend.position = "none")
    } else {
        p <- p + theme(legend.position = lp)
    }

    if (!isTRUE(input$legend_title_toggle)) {
        p <- p + theme(legend.title = element_blank())
    }

    p
})

output$pair_scatter_out <- renderPlot({
    pair_scatter_plotter()
}, height = function() (input$plot_height %||% 700) * 0.6,
   width  = function() (input$plot_width  %||% 900) * 0.7)


# ─── MODE B: Gene List Correlation Matrix ────────────────────────────────────

GeneListReactive <- reactive({
    mat <- ProcessedMatrix()
    req(mat)
    all_genes <- rownames(mat)
    src <- input$genelist_source

    genes <- character(0)

    if (src == "custom") {
        raw <- input$custom_genes
        shiny::validate(need(nchar(trimws(raw %||% "")) > 0, "Please enter at least one gene name."))
        genes <- trimws(unlist(strsplit(raw, "[,\n\r]+")))
        genes <- genes[genes != ""]
    } else if (src == "file") {
        req(input$genelist_file)
        lines <- readLines(input$genelist_file$datapath)
        genes <- trimws(unlist(strsplit(lines, "[,\n\r\t]+")))
        genes <- genes[genes != ""]
    } else if (src %in% c("hallmark", "go_bp")) {
        req(input$msigdb_geneset)
        tryCatch({
            species_sel <- input$msigdb_species %||% "Homo sapiens"
            if (src == "hallmark") {
                gs_df <- msigdbr(species = species_sel, category = "H")
            } else {
                gs_df <- msigdbr(species = species_sel, category = "C5", subcategory = "GO:BP")
            }
            genes <- gs_df$gene_symbol[gs_df$gs_name == input$msigdb_geneset]
        }, error = function(e) {
            showNotification(paste("Database error:", e$message), type='error', duration=NULL)
        })
    }

    # Intersect with available genes
    found   <- genes[genes %in% all_genes]
    missing <- genes[!genes %in% all_genes]

    if (length(missing) > 0 && length(missing) <= 20) {
        showNotification(
            paste("Genes not found:", paste(missing, collapse=", ")),
            type='warning', duration=8
        )
    } else if (length(missing) > 20) {
        showNotification(
            paste0(length(missing), " genes not found in the expression matrix."),
            type='warning', duration=8
        )
    }

    shiny::validate(need(length(found) >= 3, "Need at least 3 genes found in the matrix."))
    shiny::validate(need(length(found) <= 500, "Gene list too large (max 500). Please reduce your list."))

    found
})

MatrixCorResult <- eventReactive(input$run_matrix, {
    mat   <- ProcessedMatrix()
    genes <- GeneListReactive()
    req(mat, genes)

    # Apply min expression filter
    if (input$min_expr_filter > 0) {
        gene_means <- rowMeans(mat[genes, , drop=FALSE], na.rm=TRUE)
        genes <- genes[gene_means >= input$min_expr_filter]
        shiny::validate(need(length(genes) >= 3, "Fewer than 3 genes remain after expression filter."))
    }

    withProgress(message = "Computing correlation matrix...", value = 0.3, {
        mat_sub <- mat[genes, , drop=FALSE]
        result <- Hmisc::rcorr(t(mat_sub), type = input$cor_method)
        setProgress(1.0)
    })

    showNotification(
        paste0("Correlation matrix: ", length(genes), " genes x ", length(genes), " genes."),
        type='message', duration=4
    )

    list(r = result$r, P = result$P, genes = genes)
})

# ─── Mode B: Heatmap Plotter ─────────────────────────────────────────────────

heatmap_matrix_plotter <- reactive({
    res <- MatrixCorResult()
    req(res)

    r_mat <- res$r
    p_mat <- res$P

    # Color function
    cmin <- input$color_min %||% -1
    cmax <- input$color_max %||% 1

    col_scale <- input$color_scale %||% "RdBu"
    if (col_scale == "custom") {
        col_fun <- colorRamp2(
            c(cmin, 0, cmax),
            c(input$color_neg %||% "#4575b4", "white", input$color_pos %||% "#d73027")
        )
    } else if (col_scale == "viridis") {
        viridis_cols <- c("#440154","#482878","#3E4989","#31688E","#26828E",
                          "#1F9E89","#35B779","#6DCD59","#B4DE2C","#FDE725")
        col_fun <- colorRamp2(seq(cmin, cmax, length.out = length(viridis_cols)), viridis_cols)
    } else {
        pal_colors <- tryCatch(
            rev(brewer.pal(11, col_scale)),
            error = function(e) colorRampPalette(c("#4575b4","white","#d73027"))(11)
        )
        col_fun <- colorRamp2(seq(cmin, cmax, length.out = length(pal_colors)), pal_colors)
    }

    # Cell function for values and significance markers
    show_vals <- isTRUE(input$show_values)
    show_sig  <- isTRUE(input$show_significance)
    val_size  <- input$tick_size %||% 7

    cell_fun <- function(j, i, x, y, w, h, fill) {
        if (show_vals && i != j) {
            grid::grid.text(
                sprintf("%.2f", r_mat[i, j]),
                x, y,
                gp = grid::gpar(fontsize = val_size * 0.8, col = "black")
            )
        }
        if (show_sig && i != j && !is.na(p_mat[i, j])) {
            sig_label <- ""
            if (p_mat[i, j] < 0.001) {
                sig_label <- "***"
            } else if (p_mat[i, j] < 0.01) {
                sig_label <- "**"
            } else if (p_mat[i, j] < 0.05) {
                sig_label <- "*"
            }
            if (nchar(sig_label) > 0) {
                y_offset <- if (show_vals) unit(as.numeric(y) + 0.015, "npc") else y
                grid::grid.text(
                    sig_label, x, y_offset,
                    gp = grid::gpar(fontsize = val_size, col = "#d73027", fontface = "bold")
                )
            }
        }
    }

    ax_size  <- input$axis_size  %||% 8
    ttl_size <- input$title_size %||% 14
    bw       <- input$cell_border_width %||% 0.5

    cluster_on <- isTRUE(input$cluster_toggle)

    ht <- Heatmap(
        r_mat,
        col              = col_fun,
        name             = "Correlation",
        cell_fun         = cell_fun,
        cluster_rows     = cluster_on,
        cluster_columns  = cluster_on,
        show_row_names   = TRUE,
        show_column_names = TRUE,
        row_names_gp     = grid::gpar(fontsize = ax_size),
        column_names_gp  = grid::gpar(fontsize = ax_size),
        column_names_rot = 45,
        border           = TRUE,
        rect_gp          = grid::gpar(col = "white", lwd = bw),
        column_title     = "Gene List Correlation Matrix",
        column_title_gp  = grid::gpar(fontsize = ttl_size, fontface = "bold"),
        heatmap_legend_param = list(
            title     = if (isTRUE(input$legend_title_toggle)) "r" else NULL,
            title_gp  = grid::gpar(fontsize = 10, fontface = "bold"),
            labels_gp = grid::gpar(fontsize = 9)
        )
    )

    ht
})

# ─── Mode B: Pairwise Table ──────────────────────────────────────────────────

pairwise_table_reactive <- reactive({
    res <- MatrixCorResult()
    req(res)

    genes <- res$genes
    r_mat <- res$r
    p_mat <- res$P
    n <- length(genes)

    # Extract upper triangle
    pairs <- which(upper.tri(r_mat), arr.ind = TRUE)
    df <- data.frame(
        Gene1       = genes[pairs[, 1]],
        Gene2       = genes[pairs[, 2]],
        Correlation = round(r_mat[pairs], 4),
        P_value     = signif(p_mat[pairs], 4),
        stringsAsFactors = FALSE
    )
    df$Adj_P_value <- signif(p.adjust(df$P_value, method = "BH"), 4)
    df$Abs_Cor     <- round(abs(df$Correlation), 4)
    df <- df[order(-df$Abs_Cor), ]
    rownames(df) <- NULL

    df
})

# ─── Shared Rendering Dispatch ───────────────────────────────────────────────

output$cor_plot_out <- renderPlot({
    mode <- AnalysisMode()
    if (mode == "single") {
        # req() inside scatter_grid_plotter keeps output in spinner
        # state until SingleCorResult fires from the button click
        scatter_grid_plotter()
    } else {
        # req() inside heatmap_matrix_plotter gates on MatrixCorResult
        ht <- heatmap_matrix_plotter()
        lp <- input$legend_pos %||% "right"
        if (lp == "none") {
            draw(ht, show_heatmap_legend = FALSE, merge_legend = TRUE)
        } else {
            draw(ht, heatmap_legend_side = lp, merge_legend = TRUE)
        }
    }
}, height = function() input$plot_height %||% 700,
   width  = function() input$plot_width  %||% 900)

# ─── Correlation Table ───────────────────────────────────────────────────────

output$cor_table <- DT::renderDataTable({
    mode <- AnalysisMode()
    if (mode == "single") {
        df <- SingleCorFull()
        req(df)
        DT::datatable(df, style='bootstrap',
            options = list(pageLength = 20, scrollX = TRUE),
            filter = 'top',
            selection = 'single'
        )
    } else {
        df <- pairwise_table_reactive()
        req(df)
        DT::datatable(df, style='bootstrap',
            options = list(pageLength = 20, scrollX = TRUE),
            filter = 'top'
        )
    }
})

# ─── Download Handlers ───────────────────────────────────────────────────────

output$download_cor_plot <- downloadHandler(
    filename = function() {
        paste0("correlation_plot.", input$download_format)
    },
    content = function(file) {
        mode   <- AnalysisMode()
        h_px   <- input$plot_height %||% 700
        w_px   <- input$plot_width  %||% 900
        h_in   <- h_px / 96
        w_in   <- w_px / 96
        fmt    <- input$download_format

        if (fmt == "png") {
            png(file, height = h_px, width = w_px, res = 96)
        } else if (fmt == "jpeg") {
            jpeg(file, height = h_px, width = w_px, res = 96)
        } else if (fmt == "tiff") {
            tiff(file, height = h_px, width = w_px, res = 96)
        } else if (fmt == "pdf") {
            pdf(file, height = h_in, width = w_in)
        } else if (fmt == "svg") {
            svg(file, height = h_in, width = w_in)
        } else if (fmt == "eps") {
            setEPS()
            postscript(file, height = h_in, width = w_in)
        }

        if (mode == "single") {
            print(scatter_grid_plotter())
        } else {
            ht <- heatmap_matrix_plotter()
            draw(ht, merge_legend = TRUE)
        }

        dev.off()
    }
)

output$download_cor_csv <- downloadHandler(
    filename = function() {
        mode <- AnalysisMode()
        if (mode == "single") "single_gene_correlations.csv" else "pairwise_correlations.csv"
    },
    content = function(file) {
        mode <- AnalysisMode()
        if (mode == "single") {
            write.csv(SingleCorFull(), file, row.names = FALSE)
        } else {
            write.csv(pairwise_table_reactive(), file, row.names = FALSE)
        }
    }
)

output$download_pair_plot <- downloadHandler(
    filename = function() {
        paste0("pair_", input$target_gene, "_vs_", input$pair_gene, ".png")
    },
    content = function(file) {
        png(file, height = 500, width = 600, res = 96)
        print(pair_scatter_plotter())
        dev.off()
    }
)

# ─── Reproducible Code Modal ─────────────────────────────────────────────────

build_code_string <- function(mode, inp) {
    if (mode == "single") {
        target <- inp$target_gene %||% "GENE"
        paste0(
            "library(ggplot2)\nlibrary(ggpubr)\n\n",
            "# ── Load your expression matrix ─────────────────────────────\n",
            "# mat <- as.matrix(read.csv('expression.csv', row.names=1))\n\n",
            "# ── Single gene vs all correlation ─────────────────────────\n",
            sprintf('target_gene <- "%s"\n', target),
            "target_vec <- as.numeric(mat[target_gene, ])\n",
            "n <- ncol(mat)\n\n",
            sprintf('r_vals <- cor(t(mat), target_vec, method = "%s", use = "pairwise.complete.obs")\n', inp$cor_method %||% "pearson"),
            "t_stat <- r_vals * sqrt((n - 2) / (1 - r_vals^2))\n",
            "p_vals <- 2 * pt(-abs(t_stat), df = n - 2)\n",
            "padj   <- p.adjust(p_vals, method = \"BH\")\n\n",
            "results <- data.frame(\n",
            "    Gene = rownames(mat),\n",
            "    Correlation = round(r_vals, 4),\n",
            "    P_value = signif(p_vals, 4),\n",
            "    Adj_P_value = signif(padj, 4)\n",
            ")\n",
            "results <- results[order(-abs(results$Correlation)), ]\n\n",
            sprintf("# ── Scatter plot (top %d genes) ─────────────────────────────\n", inp$top_n %||% 20),
            sprintf("top_genes <- head(results$Gene, %d)\n\n", inp$top_n %||% 20),
            "plot_data <- do.call(rbind, lapply(top_genes, function(g) {\n",
            "    data.frame(Target = target_vec, Partner = as.numeric(mat[g, ]), Gene = g)\n",
            "}))\n",
            "plot_data$Gene <- factor(plot_data$Gene, levels = top_genes)\n\n",
            "ggplot(plot_data, aes(x = Target, y = Partner)) +\n",
            sprintf("    geom_point(size = %s, alpha = %s, color = \"#377EB8\") +\n", inp$point_size %||% 2, inp$point_alpha %||% 0.6),
            "    geom_smooth(method = \"lm\", se = TRUE, color = \"#d73027\") +\n",
            sprintf("    stat_cor(method = \"%s\") +\n", inp$cor_method %||% "pearson"),
            "    facet_wrap(~ Gene, scales = \"free_y\") +\n",
            sprintf("    labs(x = \"%s expression\", y = \"Partner gene expression\") +\n", target),
            "    theme_bw()\n"
        )
    } else {
        paste0(
            "library(Hmisc)\nlibrary(ComplexHeatmap)\nlibrary(circlize)\n\n",
            "# ── Load your expression matrix ─────────────────────────────\n",
            "# mat <- as.matrix(read.csv('expression.csv', row.names=1))\n\n",
            "# ── Define gene list ────────────────────────────────────────\n",
            "# genes <- c(\"TP53\", \"BRCA1\", \"MYC\", ...)\n",
            "# mat_sub <- mat[genes, , drop=FALSE]\n\n",
            "# ── Compute correlation matrix ──────────────────────────────\n",
            sprintf("result <- Hmisc::rcorr(t(mat_sub), type = \"%s\")\n", inp$cor_method %||% "pearson"),
            "r_mat <- result$r\n",
            "p_mat <- result$P\n\n",
            "# ── Color scale ────────────────────────────────────────────\n",
            sprintf("col_fun <- colorRamp2(c(%s, 0, %s), c(\"%s\", \"white\", \"%s\"))\n\n",
                inp$color_min %||% -1, inp$color_max %||% 1,
                inp$color_neg %||% "#4575b4", inp$color_pos %||% "#d73027"),
            "# ── Heatmap ────────────────────────────────────────────────\n",
            "ht <- Heatmap(\n",
            "    r_mat,\n",
            "    col = col_fun,\n",
            "    name = \"Correlation\",\n",
            sprintf("    cluster_rows = %s,\n", tolower(as.character(isTRUE(inp$cluster_toggle)))),
            sprintf("    cluster_columns = %s,\n", tolower(as.character(isTRUE(inp$cluster_toggle)))),
            "    column_names_rot = 45,\n",
            "    border = TRUE,\n",
            "    rect_gp = gpar(col = \"white\", lwd = 0.5)\n",
            ")\n\n",
            "draw(ht, merge_legend = TRUE)\n"
        )
    }
}

observeEvent(input$show_code_modal, {
    mode <- isolate(AnalysisMode())

    code <- tryCatch(
        build_code_string(mode, input),
        error = function(e) paste0("# Code generation error: ", conditionMessage(e))
    )

    showModal(modalDialog(
        title     = tagList(icon("file-code"), " Reproducible R Code"),
        size      = "l",
        easyClose = TRUE,
        footer    = modalButton("Close"),
        p("Copy this code to reproduce your current analysis in an offline R session.",
          style = "color:#555; margin-bottom:12px;"),
        tags$pre(
            style = paste(
                "background:#1e1e1e; color:#d4d4d4; border-radius:6px;",
                "padding:16px; font-size:12px; max-height:520px; overflow-y:auto;",
                "white-space:pre; font-family:'Courier New', monospace;"
            ),
            code
        )
    ))
})

# ─── Help Modal ───────────────────────────────────────────────────────────────

show_correlation_help_ui <- function() {
    showModal(modalDialog(
        title     = tagList(icon("circle-question"), " Correlation Tool Help"),
        size      = "l",
        easyClose = TRUE,
        footer    = modalButton("Close"),
        tabsetPanel(
            tabPanel("Overview",
                br(),
                h4("What is Gene Correlation Analysis?"),
                p("Gene correlation analysis measures the statistical relationship between expression patterns of genes across samples. Highly correlated genes tend to be co-expressed, suggesting shared regulation, functional relationships, or pathway membership."),
                h4("Correlation Methods"),
                tags$ul(
                    tags$li(strong("Pearson:"), " Measures linear relationships. Best for normally distributed data. Most commonly used."),
                    tags$li(strong("Spearman:"), " Rank-based, measures monotonic relationships. Robust to outliers and non-normal distributions."),
                    tags$li(strong("Kendall:"), " Rank-based, measures concordance. Most robust but slowest. Good for small sample sizes.")
                ),
                h4("Significance"),
                p("P-values are computed for each correlation coefficient. Adjusted p-values (Benjamini-Hochberg) correct for multiple testing. Common thresholds: * p<0.05, ** p<0.01, *** p<0.001.")
            ),
            tabPanel("Mode A: Single Gene",
                br(),
                h4("Single Gene vs All Genes"),
                p("Select a target gene and compute its correlation with every other gene in the matrix."),
                h4("Workflow"),
                tags$ol(
                    tags$li("Select a target gene from the searchable dropdown."),
                    tags$li("Set the number of top genes, minimum correlation, and p-value cutoffs."),
                    tags$li("Click 'Run Correlation' to compute."),
                    tags$li("The scatter grid shows the target gene vs each top correlated gene with regression lines."),
                    tags$li("The Correlation Table shows all genes ranked by correlation."),
                    tags$li("Use the Pair Viewer to inspect any specific gene pair in detail, optionally colored by metadata.")
                )
            ),
            tabPanel("Mode B: Gene List",
                br(),
                h4("Gene List Correlation Matrix"),
                p("Provide a gene list and compute the full pairwise correlation matrix. View as a heatmap."),
                h4("Gene List Sources"),
                tags$ul(
                    tags$li(strong("Custom:"), " Paste gene symbols separated by commas or newlines."),
                    tags$li(strong("File:"), " Upload a text file with one gene per line."),
                    tags$li(strong("MSigDB Hallmark:"), " 50 well-defined biological states and processes."),
                    tags$li(strong("MSigDB GO BP:"), " Gene Ontology Biological Process terms.")
                ),
                h4("Heatmap Options"),
                p("Toggle clustering, significance markers (*//**/***), and correlation values in cells. Choose from multiple color scales or set custom diverging colors.")
            ),
            tabPanel("Controls",
                br(),
                h4("Data Transformation"),
                tags$ul(
                    tags$li(strong("None:"), " Use raw expression values."),
                    tags$li(strong("log2(x+1):"), " Log2 transformation with pseudocount. Standard for RNA-seq counts."),
                    tags$li(strong("log10(x+1):"), " Log10 transformation with pseudocount."),
                    tags$li(strong("Z-score per gene:"), " Center and scale each gene independently.")
                ),
                h4("Min Expression Filter"),
                p("Exclude genes whose mean expression across samples falls below this threshold. Helps remove noise from lowly-expressed genes."),
                h4("Resize & Download"),
                p("Use the Resize panel to adjust plot dimensions. Download the plot in PNG, PDF, SVG, TIFF, JPEG, or EPS format. Download the correlation table as CSV.")
            )
        )
    ))
}

observeEvent(input$show_help_float, { show_correlation_help_ui() })
