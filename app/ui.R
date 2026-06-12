library(shiny)
require(shinyjs)
library(shinythemes)
require(shinycssloaders)
library(shinyWidgets)

library(DT)
library(tidyverse)
library(data.table)
library(colourpicker)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(ggpubr)
library(msigdbr)


tagList(
    tags$head(
        includeHTML(("www/GA.html")),
        tags$style(type = 'text/css','.navbar-brand{display:none;}'),
        tags$style(HTML("
            .control-group-panel {
                border: 1px solid #ddd;
                border-radius: 6px;
                padding: 10px 12px;
                margin-bottom: 10px;
                background-color: #f9f9f9;
            }
            .control-group-title {
                font-weight: bold;
                font-size: 14px;
                color: #0F344C;
                margin-bottom: 8px;
            }
            #show_help_float {
                position: fixed;
                bottom: 28px;
                right: 28px;
                z-index: 9999;
                border-radius: 50%;
                width: 46px;
                height: 46px;
                font-size: 20px;
                padding: 0;
                box-shadow: 0 3px 8px rgba(0,0,0,0.25);
            }
        "))
    ),
    ## Global always-visible help button (fixed bottom-right)
    actionButton("show_help_float", label=NULL,
        icon=icon("circle-question"),
        title="Help & documentation",
        class="btn btn-info"
    ),
    fluidPage(theme = shinytheme('yeti'),
            windowTitle = "MaGIC Gene Correlation Tool",
            useShinyjs(),
            titlePanel(
                fluidRow(
                column(2, tags$a(href='http://www.bioinformagic.io/', tags$img(height=75, src="MaGIC_Icon_0f344c.svg")), align='center'),
                column(10, fluidRow(
                    column(10, h1(strong('MaGIC Gene Correlation Tool'), align='center', style="color:#0F344C;"))
                ))
                ),
                windowTitle = "MaGIC Gene Correlation Tool"),
                tags$style(type='text/css', '.navbar{font-size:20px;}'),
                tags$style(type='text/css', '.nav-tabs{padding-bottom:20px;}'),
                tags$style(type='text/css', '.navbar-default{background-color:#0F344C;}'),
                tags$style(type='text/css', HTML('.navbar { background-color: #0F344C;}
                          .tab-panel{ background-color: #0F344C;}
                          .navbar-default .navbar-nav > .active > a,
                           .navbar-default .navbar-nav > .active > a:focus,
                           .navbar-default .navbar-nav > .active > a:hover {
                                color: white;
                                background-color: #008cba;
                            }')
                          ),
                tags$head(tags$style(".modal-dialog{ width:1300px}")),

        navbarPage(title="", id='NAVTABS',

        ## Intro Page
##########################################################################################################################################################
            tabPanel('Introduction',
                fluidRow(
                    column(2),
                    column(8,
                        column(12, align="center", style="margin-bottom:25px;",
                            h3(markdown("Welcome to the Gene Correlation & Co-expression Tool by the
                            [Molecular and Genomics Informatics Core (MaGIC)](http://www.bioinformagic.io)."))),
                        hr(),
                        h4("How to Use This Tool", style="color:#0F344C;"),
                        tags$ol(
                            tags$li(strong("Navigate to the Data Input tab."),
                                " Upload your expression matrix and optional sample metadata, or click 'Load Demo Data' to explore with a built-in example."),
                            tags$li(strong("Configure your analysis."),
                                " Choose a data transformation, optionally subset samples by metadata, and select an analysis mode."),
                            tags$li(strong("Submit your data."),
                                " Click the Submit button. The Correlation Analysis tab will become visible once data is loaded."),
                            tags$li(strong("Run correlation analysis."),
                                " Use the controls in the sidebar to configure parameters and click 'Run Correlation' or 'Run Correlation Matrix'."),
                            tags$li(strong("Explore and download."),
                                " Customize the visualizations, explore the results table, and download plots or data.")
                        ),
                        hr(),
                        h4("Two Analysis Modes", style="color:#0F344C;"),
                        fluidRow(
                            column(6,
                                div(class="control-group-panel",
                                    h5(strong("Mode A: Single Gene vs All Genes"), style="color:#0F344C;"),
                                    p("Pick one target gene and find its top correlated partners across all samples."),
                                    tags$ul(
                                        tags$li("Faceted scatter plot grid showing the target gene vs each top correlated gene"),
                                        tags$li("Regression lines and correlation coefficients per facet"),
                                        tags$li("Full ranked table of all genes by correlation"),
                                        tags$li("Single pair viewer for detailed inspection, optionally colored by metadata")
                                    )
                                )
                            ),
                            column(6,
                                div(class="control-group-panel",
                                    h5(strong("Mode B: Gene List Correlation Matrix"), style="color:#0F344C;"),
                                    p("Upload or paste a gene set to view the full pairwise correlation matrix as a heatmap."),
                                    tags$ul(
                                        tags$li("Correlation heatmap with hierarchical clustering"),
                                        tags$li("Significance markers (*, **, ***) and optional correlation values"),
                                        tags$li("Support for MSigDB Hallmark and GO gene sets"),
                                        tags$li("Pairwise correlation table with filtering")
                                    )
                                )
                            )
                        ),
                        hr(),
                        h4("Required Input Data", style="color:#0F344C;"),
                        fluidRow(
                            column(6,
                                div(class="control-group-panel",
                                    h5(strong("Expression Matrix (required)"), style="color:#0F344C;"),
                                    tags$ul(
                                        tags$li("File format: CSV or TSV"),
                                        tags$li("Rows: Genes (one gene per row)"),
                                        tags$li("Columns: Samples (one sample per column)"),
                                        tags$li("First column: Gene identifiers (symbols or IDs)"),
                                        tags$li("All remaining columns: Numeric expression values")
                                    ),
                                    tags$pre("Gene,   Sample1, Sample2, Sample3\nBRCA1,  6.5,     7.1,     5.9\nTP53,   8.2,     7.9,     8.5")
                                )
                            ),
                            column(6,
                                div(class="control-group-panel",
                                    h5(strong("Sample Metadata (optional)"), style="color:#0F344C;"),
                                    tags$ul(
                                        tags$li("File format: CSV or TSV"),
                                        tags$li("Rows: Samples (one sample per row)"),
                                        tags$li("First column: Sample names (must match matrix columns)"),
                                        tags$li("Additional columns: Metadata variables (used for sample subsetting and scatter plot coloring)")
                                    ),
                                    tags$pre("Sample,  Group,    Batch\nSample1, Control,  Batch1\nSample2, Treated,  Batch1\nSample3, Treated,  Batch2")
                                )
                            )
                        ),
                        hr()
                    ),
                    column(2)
                )
            ),


        ## Data Input Page
##########################################################################################################################################################
            tabPanel('Data Input',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('Input Data', align='center'),
                            hr(),
                            materialSwitch("DemoData", label="Upload custom data", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.DemoData",
                                h4("Expression Matrix", style="color:#0F344C;"),
                                fileInput('matrix_file', 'Upload Matrix File (CSV/TSV)',
                                    accept=c('text/csv', 'text/comma-separated-values, text/plain', '.csv',
                                             'text/tsv', 'text/tab-separated-values, text/plain', '.tsv'),
                                    multiple=FALSE
                                ),
                                h4("Sample Metadata (optional)", style="color:#0F344C;"),
                                fileInput('metadata_file', 'Upload Metadata File (CSV/TSV)',
                                    accept=c('text/csv', 'text/comma-separated-values, text/plain', '.csv',
                                             'text/tsv', 'text/tab-separated-values, text/plain', '.tsv'),
                                    multiple=FALSE
                                ),
                                hr(),
                                uiOutput('gene_col_selector'),
                                radioButtons("transform", "Data Transformation:",
                                    choices=c("None"="none", "log2(x+1)"="log2", "log10(x+1)"="log10", "Z-score per gene"="zscore"),
                                    selected="none"
                                ),
                                uiOutput('subset_ui'),
                                hr(),
                                actionButton('submit', "Submit Data", class='btn btn-info btn-block')
                            ),
                            conditionalPanel("input.DemoData==false",
                                p("Use the pre-loaded synthetic expression dataset to explore the tool's features."),
                                p(em("Demo dataset: ~3000 genes with co-expression modules across 50 samples in 3 groups.")),
                                hr(),
                                radioButtons("transform_demo", "Data Transformation:",
                                    choices=c("None"="none", "log2(x+1)"="log2", "log10(x+1)"="log10", "Z-score per gene"="zscore"),
                                    selected="none"
                                ),
                                actionButton('demo_submit', "Load Demo Data", class='btn btn-success btn-block')
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='InputTables',
                            tabPanel(title='Expression Matrix', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('matrix_table')
                                )
                            ),
                            tabPanel(title='Sample Metadata', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('metadata_table')
                                )
                            )
                        )
                    )
                )
            ),


        ## Correlation Analysis Page (hidden until data submitted)
##########################################################################################################################################################
            tabPanel('Correlation Analysis',
                fluidRow(
                    column(3,
                        wellPanel(

                            radioButtons("analysis_mode", "Analysis Mode:",
                                choices=c("Single Gene vs All"="single", "Gene List Correlation Matrix"="matrix"),
                                selected="single"
                            ),
                            hr(),

                            ## Shared Controls
                            h5(strong("Correlation Settings"), style="color:#0F344C; margin-top:4px;"),
                            hr(),
                            radioButtons("cor_method", "Correlation Method:", inline=TRUE,
                                choices=c("Pearson"="pearson", "Spearman"="spearman", "Kendall"="kendall"),
                                selected="pearson"
                            ),
                            sliderInput("min_expr_filter", "Min mean expression filter:",
                                min=0, max=10, step=0.5, value=0),
                            hr(),

                            ## Mode A Controls
                            conditionalPanel("output.current_mode == 'single'",
                                h5(strong("Single Gene vs All"), style="color:#0F344C;"),
                                selectizeInput("target_gene", "Select Target Gene:", choices=NULL,
                                    options=list(placeholder='Type gene name...')),
                                sliderInput("top_n", "Top N correlated genes:", min=5, max=50, step=5, value=20),
                                sliderInput("min_abs_cor", "Min |correlation|:", min=0, max=0.95, step=0.05, value=0.3),
                                sliderInput("pval_cutoff", "Adj. p-value cutoff:", min=0.001, max=0.1, step=0.001, value=0.05),
                                actionButton("run_single", "Run Correlation",
                                    class="btn btn-info btn-block", icon=icon("play")),
                                hr()
                            ),

                            ## Mode B Controls
                            conditionalPanel("output.current_mode == 'matrix'",
                                h5(strong("Gene List Correlation Matrix"), style="color:#0F344C;"),
                                radioButtons("genelist_source", label="Gene List Source:",
                                    choices=c(
                                        "Custom Gene List"="custom",
                                        "Upload Gene List File"="file",
                                        "MSigDB Hallmark"="hallmark",
                                        "MSigDB GO BP"="go_bp"
                                    ),
                                    selected="custom"
                                ),
                                conditionalPanel("input.genelist_source == 'custom'",
                                    textAreaInput("custom_genes", "Enter gene names (one per line or comma-separated):",
                                        rows=5, placeholder="TP53\nBRCA1\nMYC\nEGFR\n...")
                                ),
                                conditionalPanel("input.genelist_source == 'file'",
                                    fileInput('genelist_file', 'Upload gene list (one per line)',
                                        accept=c('text/plain', '.txt', '.csv'),
                                        multiple=FALSE
                                    )
                                ),
                                conditionalPanel("input.genelist_source == 'hallmark' || input.genelist_source == 'go_bp'",
                                    selectInput("msigdb_species", "Species:",
                                        choices=c("Homo sapiens"="Homo sapiens", "Mus musculus"="Mus musculus"),
                                        selected="Homo sapiens"
                                    ),
                                    selectizeInput("msigdb_geneset", "Gene Set:", choices=NULL,
                                        options=list(placeholder='Type to search gene sets...', maxOptions=5000))
                                ),
                                actionButton("run_matrix", "Run Correlation Matrix",
                                    class="btn btn-info btn-block", icon=icon("play")),
                                hr()
                            ),

                            ## Heatmap Options (Mode B)
                            conditionalPanel("output.current_mode == 'matrix'",
                                materialSwitch("show_heatmap_opts", label="Heatmap Options", value=FALSE, right=TRUE, status='info'),
                                conditionalPanel("input.show_heatmap_opts",
                                    hr(),
                                    materialSwitch("cluster_toggle", label="Enable clustering",
                                        value=TRUE, right=TRUE, status='info'),
                                    materialSwitch("show_values", label="Show correlation values",
                                        value=FALSE, right=TRUE, status='info'),
                                    materialSwitch("show_significance", label="Show significance markers",
                                        value=TRUE, right=TRUE, status='info'),
                                    selectInput("color_scale", "Color Scale:",
                                        choices=c("RdBu (diverging)"="RdBu", "RdYlBu"="RdYlBu",
                                                  "PuOr"="PuOr", "PRGn"="PRGn", "PiYG"="PiYG",
                                                  "viridis"="viridis", "Custom"="custom"),
                                        selected="RdBu"
                                    ),
                                    conditionalPanel("input.color_scale == 'custom'",
                                        column(6, colourInput("color_pos", "Positive", "#d73027")),
                                        column(6, colourInput("color_neg", "Negative", "#4575b4")),
                                        div(style="clear:both;")
                                    ),
                                    numericInput("color_min", "Color scale min:", value=-1, min=-1, max=0, step=0.1),
                                    numericInput("color_max", "Color scale max:", value=1, min=0, max=1, step=0.1),
                                    sliderInput("cell_border_width", "Cell border width:", min=0, max=3, step=0.5, value=0.5)
                                )
                            ),

                            ## Fonts
                            materialSwitch("show_fonts", label="Fonts", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_fonts",
                                hr(),
                                sliderInput("title_size", "Title size:", min=8, max=28, step=1, value=14),
                                sliderInput("axis_size", "Axis label size:", min=4, max=20, step=1, value=8),
                                sliderInput("tick_size", "Value text size:", min=4, max=16, step=1, value=7)
                            ),

                            ## Points (Mode A)
                            conditionalPanel("output.current_mode == 'single'",
                                materialSwitch("show_points", label="Points", value=FALSE, right=TRUE, status='info'),
                                conditionalPanel("input.show_points",
                                    hr(),
                                    sliderInput("point_size", "Point size:", min=0.5, max=8, step=0.5, value=2),
                                    sliderInput("point_alpha", "Point opacity:", min=0.1, max=1, step=0.1, value=0.6)
                                )
                            ),

                            ## Legend
                            materialSwitch("show_legend_opts", label="Legend", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_legend_opts",
                                hr(),
                                radioButtons("legend_pos", "Position:", inline=TRUE,
                                    choices=c("Right"="right", "Bottom"="bottom", "None"="none"),
                                    selected="right"
                                ),
                                materialSwitch("legend_title_toggle", label="Show legend title",
                                    value=TRUE, right=TRUE, status='info')
                            ),

                            ## Resize
                            materialSwitch("show_resize", label="Resize Plot", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_resize",
                                hr(),
                                sliderInput("plot_height", "Plot height (px):", min=200, max=2000, step=50, value=700),
                                sliderInput("plot_width",  "Plot width (px):",  min=200, max=2000, step=50, value=900)
                            )

                        )# end wellPanel sidebar
                    ),
                    column(9,
                        tabsetPanel(id='CorrelationTabs',
                            tabPanel(title='Correlation Plot', hr(),
                                fluidRow(style="margin: 0 8px 4px 0;",
                                    column(12, align="right",
                                        actionButton("show_code_modal", label=NULL,
                                            icon=icon("file-code"),
                                            title="View R code to reproduce this plot",
                                            class="btn btn-default btn-sm",
                                            style="border-radius:6px; font-size:16px; padding:4px 8px;"
                                        )
                                    )
                                ),
                                hr(),
                                div(style="overflow-x:auto; width:100%;",
                                    withSpinner(type=6, color='#5bc0de',
                                        plotOutput("cor_plot_out", height='100%')
                                    )
                                ),
                                div(style="margin-top:30px; text-align:center; padding-bottom:50px;",
                                    div(style="display:inline-block; width:250px; margin-bottom:10px;",
                                        selectInput("download_format", "Download format:",
                                            choices=c('png','pdf','svg','tiff','jpeg','eps'))
                                    ),
                                    br(),
                                    downloadButton('download_cor_plot', 'Download Plot')
                                )
                            ),
                            tabPanel(title='Correlation Table', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('cor_table')
                                ),
                                div(style="margin-top:20px; text-align:center; padding-bottom:30px;",
                                    downloadButton('download_cor_csv', 'Download Table (CSV)')
                                )
                            ),
                            tabPanel(title='Pair Viewer', hr(),
                                uiOutput("pair_viewer_ui"),
                                div(style="overflow-x:auto; width:100%;",
                                    withSpinner(type=6, color='#5bc0de',
                                        plotOutput("pair_scatter_out", height='100%')
                                    )
                                ),
                                div(style="margin-top:20px; text-align:center; padding-bottom:30px;",
                                    downloadButton('download_pair_plot', 'Download Pair Plot (PNG)')
                                )
                            )
                        )
                    )
                )
            ),


        ## Footer
##########################################################################################################################################################
            tags$footer(
                wellPanel(
                    fluidRow(
                        column(4, align='center',
                        tags$a(href="https://github.com/MaGIC-Analytics/magic-correlation", icon("github", "fa-3x")),
                        tags$h4('GitHub to submit issues/requests')
                        ),
                        column(4, align='center',
                        tags$a(href="http://www.bioinformagic.io/", icon("magic", "fa-3x")),
                        tags$h4('MaGIC Home Page')
                        ),
                        column(4, align='center',
                        tags$a(href="https://github.com/MaGIC-Analytics", icon("address-card", "fa-3x")),
                        tags$h4("Developer's Page")
                        )
                    ),
                    fluidRow(
                        column(12, align='center',
                            HTML('<a href="https://www.youtube.com/watch?v=dQw4w9WgXcQ">
                            <p>&copy;
                                <script language="javascript" type="text/javascript">
                                var today = new Date()
                                var year = today.getFullYear()
                                document.write(year)
                                </script>
                            </p>
                            </a>
                            ')
                        )
                    )
                )
            )
        )# Ends navbarPage
    )# Ends fluidPage
)# Ends tagList
