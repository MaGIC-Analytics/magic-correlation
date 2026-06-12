# ─── Utilities ─────────────────────────────────────────────────────────────────

`%||%` <- function(a, b) {
    if (is.null(a)) return(b)
    if (length(a) == 0) return(b)
    if (length(a) == 1 && is.character(a) && !nzchar(a)) return(b)
    a
}

read_delim_auto <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("tsv", "txt")) {
        fread(path, sep="\t")
    } else {
        fread(path, sep=",")
    }
}

# ── msigdbr version-compatibility shim ─────────────────────────────────────────
# msigdbr >= 10 renamed category/subcategory -> collection/subcollection and split
# KEGG into CP:KEGG_LEGACY/CP:KEGG_MEDICUS. Take the legacy argument names, try the
# new API first, fall back to the old one. gs_name / gene_symbol are stable.
msigdbr_compat <- function(species, category, subcategory = NULL) {
    sub_new <- if (identical(subcategory, "CP:KEGG")) "CP:KEGG_LEGACY" else subcategory
    tryCatch(
        do.call(msigdbr, c(list(species = species, collection = category),
                           if (!is.null(sub_new)) list(subcollection = sub_new))),
        error = function(e)
            do.call(msigdbr, c(list(species = species, category = category),
                               if (!is.null(subcategory)) list(subcategory = subcategory)))
    )
}

# ─── Demo Data Generator ───���─────────────────────────���───────────────────────

DemoDataCache <- reactiveVal(NULL)

generate_demo_data <- function() {
    set.seed(42)
    n_samples <- 50
    n_modules <- 5
    genes_per_module <- 20
    n_background <- 2900

    # Sample grouping
    groups <- rep(c("Control", "TreatA", "TreatB"), times = c(17, 17, 16))
    batches <- rep(c("Batch1", "Batch2"), length.out = n_samples)
    sex <- rep(c("Male", "Female"), length.out = n_samples)
    sample_names <- paste0("Sample_", sprintf("%02d", 1:n_samples))

    # Module latent factors (n_modules x n_samples)
    z <- matrix(rnorm(n_modules * n_samples), nrow = n_modules)
    # Module 1 (Cell Cycle) elevated in TreatA
    z[1, groups == "TreatA"] <- z[1, groups == "TreatA"] + 2.5
    # Module 2 (Immune) elevated in TreatB
    z[2, groups == "TreatB"] <- z[2, groups == "TreatB"] + 2.5
    # Module 3 (Metabolism) slightly elevated in Control
    z[3, groups == "Control"] <- z[3, groups == "Control"] + 1.5
    # Module 4 (Anti-cell-cycle) anti-correlated with Module 1
    z[4, ] <- -z[1, ] + rnorm(n_samples, sd = 0.3)
    # Module 5 (ECM) independent, slightly elevated in TreatB
    z[5, groups == "TreatB"] <- z[5, groups == "TreatB"] + 1.0

    # Real gene symbols for each module
    module_genes <- list(
        c("CDK1","CCNB1","AURKA","PLK1","TOP2A","MKI67","BUB1","CCNA2","CDC20","BIRC5",
          "KIF11","CENPA","NDC80","TPX2","FOXM1","MYBL2","MCM2","MCM5","PCNA","RFC4"),
        c("TNF","IL6","CCL2","CXCL8","TLR4","NFKB1","IL1B","ICAM1","VCAM1","CCL5",
          "CXCL10","IRF1","STAT1","MYD88","CD14","LY96","IRAK1","RIPK2","CASP1","IL18"),
        c("GAPDH","PKM","LDHA","HK2","ENO1","PGK1","PFKL","GPI","ALDOA","TPI1",
          "PGAM1","MDH2","IDH1","CS","ACO2","SDHA","SDHB","FH","OGDH","DLST"),
        c("RB1","CDKN1A","CDKN2A","CDKN1B","CDKN2B","GADD45A","GADD45B","TP53","BTG2","PTEN",
          "ATM","ATR","BRCA1","BRCA2","RAD51","CHEK1","CHEK2","MDM2","BAX","BBC3"),
        c("COL1A1","COL1A2","COL3A1","COL5A1","FN1","VIM","ACTA2","TGFB1","TGFB2","CTGF",
          "MMP2","MMP9","TIMP1","TIMP2","LOX","LOXL2","SPARC","POSTN","THBS1","BGN")
    )

    all_module_genes <- unlist(module_genes)
    n_total <- length(all_module_genes) + n_background

    # Build expression matrix
    expr_mat <- matrix(0, nrow = n_total, ncol = n_samples)
    gene_names <- character(n_total)

    # Module genes: base + loading * z_module + noise
    row_idx <- 1
    for (m in 1:n_modules) {
        for (g in 1:genes_per_module) {
            loading <- runif(1, 0.6, 0.95)
            base <- abs(rnorm(1, mean = 8, sd = 2))
            expr_mat[row_idx, ] <- base + loading * z[m, ] + rnorm(n_samples, sd = 0.3)
            gene_names[row_idx] <- module_genes[[m]][g]
            row_idx <- row_idx + 1
        }
    }

    # Background genes: no module structure
    bg_prefixes <- paste0("GENE", seq_len(n_background))
    for (g in seq_len(n_background)) {
        base <- abs(rnorm(1, mean = 7, sd = 2.5))
        expr_mat[row_idx, ] <- base + rnorm(n_samples, sd = 0.5)
        gene_names[row_idx] <- bg_prefixes[g]
        row_idx <- row_idx + 1
    }

    # Ensure non-negative values
    expr_mat[expr_mat < 0] <- 0.01

    # Assemble data.tables
    mat_dt <- data.table(Gene = gene_names, as.data.table(expr_mat))
    setnames(mat_dt, c("Gene", sample_names))

    meta_dt <- data.table(
        Sample = sample_names,
        Group  = groups,
        Batch  = batches,
        Sex    = sex
    )

    list(matrix = mat_dt, metadata = meta_dt)
}

# ─── Data Loading Reactives ────────���──────────────────────────────────────────

MatrixReactive <- reactive({
    if (input$DemoData == FALSE) {
        # Demo mode
        cached <- DemoDataCache()
        if (is.null(cached)) {
            cached <- generate_demo_data()
            DemoDataCache(cached)
        }
        cached$matrix
    } else {
        shiny::validate(need(!is.null(input$matrix_file), "Please upload an expression matrix file."))
        tryCatch(
            read_delim_auto(input$matrix_file$datapath),
            error = function(e) { showNotification(paste("Matrix parse error:", e), type='error', duration=NULL); NULL }
        )
    }
})

MetadataReactive <- reactive({
    if (input$DemoData == FALSE) {
        # Demo mode
        cached <- DemoDataCache()
        if (is.null(cached)) {
            cached <- generate_demo_data()
            DemoDataCache(cached)
        }
        cached$metadata
    } else {
        # Custom upload - metadata is optional
        if (is.null(input$metadata_file)) return(NULL)
        tryCatch(
            read_delim_auto(input$metadata_file$datapath),
            error = function(e) { showNotification(paste("Metadata parse error:", e), type='error', duration=NULL); NULL }
        )
    }
})

# ─── Gene Column Selector ────────��───────────────────────────────────────────

output$gene_col_selector <- renderUI({
    req(input$matrix_file)
    mat <- MatrixReactive()
    req(mat)
    tagList(
        selectInput("gene_col", "Which column contains gene names?",
            choices = colnames(mat), selected = colnames(mat)[1]),
        hr()
    )
})

# ─── Sample Subsetting UI ─────────────────────────────────────���──────────────

output$subset_ui <- renderUI({
    meta <- MetadataReactive()
    if (is.null(meta)) return(NULL)
    meta_cols <- colnames(meta)[-1]
    if (length(meta_cols) == 0) return(NULL)
    tagList(
        hr(),
        h5(strong("Sample Subsetting (optional)"), style="color:#0F344C;"),
        selectInput("subset_col", "Subset by metadata column:",
            choices = c("(none)" = "none", meta_cols),
            selected = "none"
        ),
        uiOutput("subset_levels_ui")
    )
})

output$subset_levels_ui <- renderUI({
    req(input$subset_col)
    if (input$subset_col == "none") return(NULL)
    meta <- MetadataReactive()
    req(meta)
    vals <- unique(as.character(meta[[input$subset_col]]))
    vals <- vals[!is.na(vals)]
    checkboxGroupInput("subset_levels", "Include levels:",
        choices = vals, selected = vals)
})

# ─── Preview Tables ──��─────────────────────────────────────���─────────────────

output$matrix_table <- DT::renderDataTable({
    mat <- MatrixReactive()
    req(mat)
    DT::datatable(mat, style='bootstrap', options=list(pageLength=15, scrollX=TRUE))
})

output$metadata_table <- DT::renderDataTable({
    meta <- MetadataReactive()
    req(meta)
    DT::datatable(meta, style='bootstrap', options=list(pageLength=15, scrollX=TRUE))
})

# ���── Processed Matrix (numeric, transformed, filtered) ───────────────────────

ProcessedMatrix <- reactive({
    mat <- MatrixReactive()
    req(mat)

    # Determine gene column
    gene_col <- if (!is.null(input$gene_col) && input$DemoData) input$gene_col else colnames(mat)[1]

    gene_names <- as.character(mat[[gene_col]])
    num_cols   <- setdiff(colnames(mat), gene_col)

    m <- as.matrix(mat[, ..num_cols])
    rownames(m) <- make.unique(gene_names)
    storage.mode(m) <- "numeric"

    # Apply sample subsetting
    if (!is.null(input$subset_col) && input$subset_col != "none" && !is.null(input$subset_levels)) {
        meta <- MetadataReactive()
        if (!is.null(meta)) {
            sample_col   <- colnames(meta)[1]
            keep_samples <- as.character(meta[[sample_col]][ meta[[input$subset_col]] %in% input$subset_levels ])
            keep_samples <- intersect(keep_samples, colnames(m))
            if (length(keep_samples) >= 3) {
                m <- m[, keep_samples, drop=FALSE]
            } else {
                showNotification("Subsetting would leave fewer than 3 samples. Ignoring subset.", type='warning', duration=5)
            }
        }
    }

    # Apply transformation
    trans <- if (input$DemoData) input$transform else (input$transform_demo %||% "none")
    if (trans == "log2") {
        m <- log2(m + 1)
    } else if (trans == "log10") {
        m <- log10(m + 1)
    } else if (trans == "zscore") {
        m <- t(scale(t(m)))
        m[is.nan(m)] <- 0
    }

    m
})

# ─── Processed Metadata (aligned to matrix columns) ──────────────────────────

ProcessedMeta <- reactive({
    meta  <- MetadataReactive()
    mat_m <- ProcessedMatrix()
    if (is.null(meta)) return(NULL)
    req(mat_m)

    sample_col   <- colnames(meta)[1]
    sample_names <- as.character(meta[[sample_col]])

    col_order <- colnames(mat_m)
    idx       <- match(col_order, sample_names)

    if (any(is.na(idx))) {
        showNotification(
            paste("Warning: Some samples not found in metadata:",
                  paste(col_order[is.na(idx)], collapse=", ")),
            type='warning', duration=8
        )
    }

    df <- as.data.frame(meta)[idx, , drop=FALSE]
    rownames(df) <- col_order
    df
})

# ─── Analysis Mode Reactive ───────────────────────────────────────────────���──

AnalysisMode <- reactive({
    # Single mode control, now on the Correlation Analysis tab (not split demo/upload).
    input$analysis_mode %||% "single"
})

output$current_mode <- reactive({ AnalysisMode() })
outputOptions(output, "current_mode", suspendWhenHidden = FALSE)

# ─── MSigDB Gene Set Selector ────────────────────────────────────────────────

observe({
    req(input$genelist_source %in% c("hallmark", "go_bp"))
    tryCatch({
        species_sel <- input$msigdb_species %||% "Homo sapiens"
        src <- input$genelist_source

        if (src == "hallmark") {
            gs_df <- msigdbr_compat(species_sel, "H")
        } else if (src == "go_bp") {
            gs_df <- msigdbr_compat(species_sel, "C5", "GO:BP")
        }

        gene_sets <- sort(unique(gs_df$gs_name))
        updateSelectizeInput(session, "msigdb_geneset", choices = gene_sets, server = TRUE,
            options = list(placeholder = 'Type to search gene sets...'))
    }, error = function(e) {
        showNotification(paste("Database query error:", e$message), type='error', duration=8)
    })
})
