#' Launch EnrichGT Web Interface
#'
#' Launch an interactive Shiny web application that provides a user-friendly interface
#' for performing enrichment analysis using EnrichGT package functions.
#'
#' @param port Port number for the application. If NULL (default), a random available port will be used.
#' @param host IP address that the application should listen on. Default is "127.0.0.1" (localhost).
#' @param launch.browser Logical. If TRUE (default), the web browser will be launched automatically.
#' @param LLM An optional LLM chat object created by the \code{ellmer} package.
#'   If provided, enables LLM-powered summarization of recluster analysis results.
#'   If NULL (default), LLM features are disabled.
#'
#' @details
#' The web interface includes three main analysis modules:
#'
#' **Enrichment Analysis (ORA):**
#' - Input: Gene list (one gene per line or comma-separated)
#' - Analysis: Over-representation analysis using \code{egt_enrichment_analysis}
#' - Output: Interactive plot and downloadable table
#'
#' **GSEA Analysis:**
#' - Input: Ranked gene list with weights (tab-separated: gene\ttab\tweight)
#' - Analysis: Gene Set Enrichment Analysis using \code{egt_gsea_analysis}
#' - Output: Interactive plot and downloadable table
#'
#' **Recluster Analysis:**
#' - Input: Results from ORA or GSEA analysis
#' - Analysis: Hierarchical clustering of enrichment results using \code{egt_recluster_analysis}
#' - Output: Interactive gt table with clustered results and optional LLM summary
#'
#' **Available Databases:**
#' - GO Biological Process (GO_BP)
#' - GO Cellular Component (GO_CC)
#' - GO Molecular Function (GO_MF)
#' - KEGG Pathways
#' - Reactome Pathways
#'
#' Each analysis module allows users to adjust parameters such as p-value cutoffs,
#' geneset size limits, and visualization options.
#'
#' @return Starts a Shiny application
#'
#' @examples
#' \dontrun{
#' # Launch the web interface with default settings
#' egt_web_interface()
#'
#' # Launch with LLM support
#' library(ellmer)
#' chat <- chat_deepseek(api_key = "your_api_key", model = "deepseek-chat")
#' egt_web_interface(LLM = chat)
#'
#' # Launch on a specific port
#' egt_web_interface(port = 3838)
#'
#' # Launch without opening browser automatically
#' egt_web_interface(launch.browser = FALSE)
#' }
#'
#' @seealso
#' \code{\link{egt_enrichment_analysis}}, \code{\link{egt_gsea_analysis}},
#' \code{\link{egt_recluster_analysis}}, \code{\link{egt_plot_results}}, \code{\link{egt_llm_summary}}
#'
#' @author Zhiming Ye
#' @export
egt_web_interface <- function(port = NULL, host = "127.0.0.1", launch.browser = TRUE, LLM = NULL) {

  # Check required packages
  required_packages <- c("shiny", "shinydashboard", "DT", "plotly", "gt")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
         "\nPlease install them using: install.packages(c('shiny', 'shinydashboard', 'DT', 'plotly', 'gt'))")
  }

  suppressPackageStartupMessages({
    library(shiny)
    library(shinydashboard)
    library(DT)
    library(plotly)
    library(gt)
    library(EnrichGT)
    library(dplyr)
  })

  # UI Definition
  ui <- dashboardPage(
    dashboardHeader(title = "EnrichGT Web Interface"),

    dashboardSidebar(
      sidebarMenu(
        menuItem("Enrichment Analysis (ORA)", tabName = "ora", icon = icon("search")),
        menuItem("GSEA Analysis", tabName = "gsea", icon = icon("chart-line")),
        menuItem("Recluster Analysis", tabName = "recluster", icon = icon("project-diagram")),
        menuItem("Help", tabName = "help", icon = icon("question-circle"))
      )
    ),

    dashboardBody(
      tags$head(
        tags$style(HTML("
          .content-wrapper, .right-side {
            background-color: #f4f4f4;
          }
          .box {
            border-radius: 5px;
          }
          .nav-tabs-custom > .nav-tabs > li.active {
            border-top-color: #3c8dbc;
          }
        "))
      ),

      tabItems(
        # ORA Analysis Tab
        tabItem(tabName = "ora",
          fluidRow(
            box(
              title = "Over-Representation Analysis (ORA)", status = "primary", solidHeader = TRUE,
              width = 4, height = "800px",

              h4("Input Parameters"),

              textAreaInput("ora_genes",
                           label = "Gene List (one per line or comma-separated):",
                           placeholder = "TP53\nBRCA1\nEGFR\nCDK2",
                           rows = 8),

              selectInput("ora_species", "Species:",
                         choices = list(
                           "Human" = "org.Hs.eg.db",
                           "Mouse" = "org.Mm.eg.db"
                         ),
                         selected = "org.Hs.eg.db"),

              selectInput("ora_database", "Database:",
                         choices = list(
                           "GO Biological Process" = "GO_BP",
                           "GO Molecular Function" = "GO_MF",
                           "GO Cellular Component" = "GO_CC",
                           "KEGG Pathways" = "KEGG",
                           "Reactome" = "Reactome"
                         ),
                         selected = "GO_BP"),

              numericInput("ora_pval", "P-value Cutoff:", value = 0.05, min = 0.001, max = 1, step = 0.01),
              numericInput("ora_min_size", "Min Gene Set Size:", value = 10, min = 1, max = 100),
              numericInput("ora_max_size", "Max Gene Set Size:", value = 500, min = 10, max = 2000),

              br(),
              actionButton("run_ora", "Run ORA Analysis", class = "btn-primary", width = "100%")
            ),

            box(
              title = "ORA Results", status = "success", solidHeader = TRUE,
              width = 8, height = "600px",

              tabsetPanel(
                tabPanel("Plot",
                         div(style = "height: 500px; overflow: auto;",
                             plotOutput("ora_plot", height = "480px")
                         )
                ),
                tabPanel("Data Table",
                         div(style = "height: 500px; overflow: auto;",
                             DT::dataTableOutput("ora_table")
                         )
                )
              )
            )
          )
        ),

        # GSEA Analysis Tab
        tabItem(tabName = "gsea",
          fluidRow(
            box(
              title = "Gene Set Enrichment Analysis (GSEA)", status = "primary", solidHeader = TRUE,
              width = 4, height = "850px",

              h4("Input Parameters"),

              textAreaInput("gsea_genes",
                           label = "Ranked Gene List (Gene\tScore format):",
                           placeholder = "TP53\t2.5\nBRCA1\t-1.8\nEGFR\t1.2\nCDK2\t-0.9",
                           rows = 8),

              selectInput("gsea_species", "Species:",
                         choices = list(
                           "Human" = "org.Hs.eg.db",
                           "Mouse" = "org.Mm.eg.db"
                         ),
                         selected = "org.Hs.eg.db"),

              selectInput("gsea_database", "Database:",
                         choices = list(
                           "GO Biological Process" = "GO_BP",
                           "GO Molecular Function" = "GO_MF",
                           "GO Cellular Component" = "GO_CC",
                           "KEGG Pathways" = "KEGG",
                           "Reactome" = "Reactome"
                         ),
                         selected = "GO_BP"),

              numericInput("gsea_pval", "P-value Cutoff:", value = 0.05, min = 0.001, max = 1, step = 0.01),
              numericInput("gsea_min_size", "Min Gene Set Size:", value = 10, min = 1, max = 100),
              numericInput("gsea_max_size", "Max Gene Set Size:", value = 500, min = 10, max = 2000),
              numericInput("gsea_param", "GSEA Parameter:", value = 1, min = 0, max = 2, step = 0.1),

              br(),
              actionButton("run_gsea", "Run GSEA Analysis", class = "btn-primary", width = "100%")
            ),

            box(
              title = "GSEA Results", status = "success", solidHeader = TRUE,
              width = 8, height = "600px",

              tabsetPanel(
                tabPanel("Plot",
                         div(style = "height: 500px; overflow: auto;",
                             plotOutput("gsea_plot", height = "480px")
                         )
                ),
                tabPanel("Data Table",
                         div(style = "height: 500px; overflow: auto;",
                             DT::dataTableOutput("gsea_table")
                         )
                )
              )
            )
          )
        ),

        # Recluster Analysis Tab
        tabItem(tabName = "recluster",
          fluidRow(
            box(
              title = "Recluster Analysis", status = "primary", solidHeader = TRUE,
              width = 4, height = "700px",

              h4("Input Parameters"),

              radioButtons("recluster_input_type", "Input Type:",
                          choices = list(
                            "Use ORA results from above" = "ora",
                            "Use GSEA results from above" = "gsea"
                          ),
                          selected = "ora"),

              numericInput("recluster_num", "Number of Clusters:", value = 10, min = 2, max = 50),
              numericInput("recluster_pval", "P-value Cutoff:", value = 0.05, min = 0.001, max = 1, step = 0.01),
              numericInput("recluster_ntop", "Top Terms per Cluster:", value = 10, min = 1, max = 50),

              selectInput("recluster_method", "Clustering Method:",
                         choices = list(
                           "Ward D2" = "ward.D2",
                           "Ward D" = "ward.D",
                           "Complete" = "complete",
                           "Average" = "average",
                           "Single" = "single"
                         ),
                         selected = "ward.D2"),

              # Show LLM options only if LLM object is provided
               if (!is.null(LLM)) {
                 list(
                   checkboxInput("recluster_llm_summary", "Generate LLM Summary",
                                value = FALSE),
                   helpText("Note: LLM summarization may take additional time to complete.")
                 )
               },

              br(),
              actionButton("run_recluster", "Run Recluster Analysis", class = "btn-primary", width = "100%")
            ),

            box(
              title = "Recluster Results", status = "success", solidHeader = TRUE,
              width = 8, height = "600px",

              tabsetPanel(
                tabPanel("Clustered Table",
                         div(style = "height: 500px; overflow: auto;",
                             gt::gt_output("recluster_table")
                         )
                ),
                tabPanel("LLM Summary",
                         div(style = "height: 500px; overflow: auto; padding: 15px; background-color: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px;",
                             htmlOutput("recluster_llm_summary")
                         )
                )
              )
            )
          )
        ),

        # Help Tab
        tabItem(tabName = "help",
          fluidRow(
            box(
              title = "How to Use EnrichGT Web Interface", status = "info", solidHeader = TRUE,
              width = 12,

              h3("Over-Representation Analysis (ORA)"),
              p("Input a list of genes (one per line or comma-separated) to identify enriched pathways."),
              tags$ul(
                tags$li("Gene format: Gene symbols (e.g., TP53, BRCA1)"),
                tags$li("Select appropriate database for your analysis"),
                tags$li("Adjust p-value cutoff and gene set size filters as needed")
              ),

              h3("Gene Set Enrichment Analysis (GSEA)"),
              p("Input ranked genes with scores to identify enriched pathways at the top/bottom of the ranking."),
              tags$ul(
                tags$li("Format: Gene\tScore (tab-separated, e.g., TP53\t2.5)"),
                tags$li("Scores can be log2 fold changes, t-statistics, etc."),
                tags$li("Higher scores = more upregulated genes")
              ),

              h3("Recluster Analysis"),
              p("Reduce redundancy in enrichment results by clustering related terms."),
              tags$ul(
                tags$li("Use results from ORA or GSEA analysis above"),
                tags$li("Adjust cluster number based on your data"),
                tags$li("Output shows clustered terms in an interactive table")
              ),

              h3("Available Species"),
              tags$ul(
                tags$li("Human: Homo sapiens (org.Hs.eg.db)"),
                tags$li("Mouse: Mus musculus (org.Mm.eg.db)")
              ),

              h3("Databases Available"),
              tags$ul(
                tags$li("GO Biological Process: Biological processes and pathways"),
                tags$li("GO Molecular Function: Molecular activities and functions"),
                tags$li("GO Cellular Component: Cellular locations and structures"),
                tags$li("KEGG Pathways: Kyoto Encyclopedia pathways (hsa for human, mmu for mouse)"),
                tags$li("Reactome: Reactome pathway database")
              ),

              h3("LLM Features"),
              if (!is.null(LLM)) {
                p("LLM summarization is enabled for Recluster Analysis. This feature provides AI-powered biological interpretation of clustered pathways and genes.")
              } else {
                p("LLM summarization is not available. To enable this feature, provide an ellmer Chat object when calling egt_web_interface(LLM = your_chat_object).")
              }
            )
          )
        )
      )
    )
  )

  # Server Logic
  server <- function(input, output, session) {

    # Reactive values to store results
    values <- reactiveValues(
      ora_result = NULL,
      gsea_result = NULL,
      recluster_result = NULL,
      llm_summary = NULL
    )

    # Helper function to get database
     get_database <- function(db_name, species = "org.Hs.eg.db") {
       # Convert species to OrgDB object
       eval(parse(text = paste0("library(", species, ")")))
       orgdb_obj <- eval(parse(text = species))

       # Get KEGG organism code
       kegg_org <- if(species == "org.Hs.eg.db") "hsa" else "mmu"

       tryCatch({
         switch(db_name,
           "GO_BP" = database_GO_BP(OrgDB = orgdb_obj),
           "GO_MF" = database_GO_MF(OrgDB = orgdb_obj),
           "GO_CC" = database_GO_CC(OrgDB = orgdb_obj),
           "KEGG" = database_KEGG(kegg_organism = kegg_org, OrgDB = orgdb_obj),
           "Reactome" = database_Reactome(OrgDB = orgdb_obj)
         )
       }, error = function(e) {
         showNotification(paste("Error loading database:", e$message), type = "error")
         return(NULL)
       })
     }

    # Parse gene input
    parse_gene_list <- function(input_text) {
      genes <- trimws(unlist(strsplit(input_text, "[\n,;\t]+")))
      genes <- genes[genes != ""]
      return(genes)
    }

    # Parse ranked gene input for GSEA
    parse_ranked_genes <- function(input_text) {
      lines <- trimws(unlist(strsplit(input_text, "\n")))
      lines <- lines[lines != ""]

      gene_scores <- list()
      for (line in lines) {
        parts <- trimws(unlist(strsplit(line, "\t")))
        if (length(parts) >= 2) {
          gene <- parts[1]
          score <- as.numeric(parts[2])
          if (!is.na(score)) {
            gene_scores[[gene]] <- score
          }
        }
      }

      if (length(gene_scores) == 0) {
        stop("No valid gene-score pairs found. Please use Gene\tScore format.")
      }

      return(unlist(gene_scores))
    }

    # ORA Analysis
    observeEvent(input$run_ora, {
      tryCatch({
        # Show progress
        showNotification("Running ORA analysis...", type = "message", duration = 3)

        # Parse genes
        genes <- parse_gene_list(input$ora_genes)
        if (length(genes) == 0) {
          showNotification("Please enter gene list", type = "error")
          return()
        }

        # Get database
        database <- get_database(input$ora_database, input$ora_species)

        # Run analysis
        result <- egt_enrichment_analysis(
          genes = genes,
          database = database,
          p_val_cut_off = input$ora_pval,
          min_geneset_size = input$ora_min_size,
          max_geneset_size = input$ora_max_size
        )

        # Check for significant results based on p.adjust values
        if ("p.adjust" %in% colnames(result)) {
          significant_results <- sum(result$p.adjust < input$ora_pval, na.rm = TRUE)
          if (significant_results == 0) {
            showNotification(paste("No significant results found. All adjusted p-values ≥", input$ora_pval, ". Try relaxing the p-value cutoff."), type = "warning", duration = 10)
            values$ora_result <- result
            return()
          }
        } else if (nrow(result) == 0) {
          showNotification("No significant results found. Try relaxing parameters.", type = "warning")
          values$ora_result <- result
          return()
        }

        values$ora_result <- result
        showNotification(paste("ORA analysis completed! Found", nrow(result), "significant pathways."), type = "message")

      }, error = function(e) {
        showNotification(paste("Error in ORA analysis:", e$message), type = "error", duration = 10)
      })
    })

    # GSEA Analysis
    observeEvent(input$run_gsea, {
      tryCatch({
        # Show progress
        showNotification("Running GSEA analysis...", type = "default", duration = 3)

        # Parse ranked genes
        ranked_genes <- parse_ranked_genes(input$gsea_genes)
        if (length(ranked_genes) == 0) {
          showNotification("Please enter ranked gene list", type = "error")
          return()
        }

        # Get database
        database <- get_database(input$gsea_database, input$gsea_species)

        # Run analysis
        result <- egt_gsea_analysis(
          genes = ranked_genes,
          database = database,
          p_val_cut_off = input$gsea_pval,
          min_geneset_size = input$gsea_min_size,
          max_geneset_size = input$gsea_max_size,
          gseaParam = input$gsea_param
        )

        # Check for significant results based on p.adjust values
        if ("p.adjust" %in% colnames(result)) {
          significant_results <- sum(result$p.adjust < input$gsea_pval, na.rm = TRUE)
          if (significant_results == 0) {
            showNotification(paste("No significant results found. All adjusted p-values ≥", input$gsea_pval, ". Try relaxing the p-value cutoff."), type = "warning", duration = 10)
            values$gsea_result <- NULL
            return()
          }
        } else if (nrow(result) == 0) {
          showNotification("No significant results found. Try relaxing parameters.", type = "warning")
          values$gsea_result <- NULL
          return()
        }

        values$gsea_result <- result
        showNotification(paste("GSEA analysis completed! Found", nrow(result), "significant pathways."), type = "message")

      }, error = function(e) {
        showNotification(paste("Error in GSEA analysis:", e$message), type = "error", duration = 10)
      })
    })

    # Recluster Analysis
    observeEvent(input$run_recluster, {
      tryCatch({
        # Get input data
        input_data <- NULL
        if (input$recluster_input_type == "ora") {
          input_data <- values$ora_result
          if (is.null(input_data)) {
            showNotification("Please run ORA analysis first", type = "error")
            return()
          }
        } else if (input$recluster_input_type == "gsea") {
          input_data <- values$gsea_result
          if (is.null(input_data)) {
            showNotification("Please run GSEA analysis first", type = "error")
            return()
          }
        }

        # Show progress
        showNotification("Running recluster analysis...", type = "message", duration = 3)

        # Run recluster analysis
        result <- egt_recluster_analysis(
          x = input_data,
          ClusterNum = input$recluster_num,
          P.adj = input$recluster_pval,
          nTop = input$recluster_ntop,
          method = input$recluster_method
        )

        values$recluster_result <- result

        # Generate LLM summary if requested and LLM is available
        if (!is.null(LLM) && isTRUE(input$recluster_llm_summary)) {
          tryCatch({
            showNotification("Generating LLM summary...", type = "default", duration = 3)
            showNotification("LLM summary is processing, it would be very slow...", type = "default", duration = 15)
            llm_summary <- egt_llm_summary(result, chat = LLM)
            values$llm_summary <- llm_summary
            showNotification("Recluster analysis completed with LLM summary!", type = "message")
          }, error = function(e) {
            values$llm_summary <- paste("Error generating LLM summary:", e$message)
            showNotification("Recluster analysis completed, but LLM summary failed.", type = "warning")
          })
        } else if (isTRUE(input$recluster_llm_summary) && is.null(LLM)) {
          values$llm_summary <- "LLM object not provided. To enable LLM summary, please provide an LLM object when calling egt_web_interface()."
          showNotification("Recluster analysis completed. LLM not available for summary.", type = "warning")
        } else {
          values$llm_summary <- NULL
          showNotification("Recluster analysis completed!", type = "message")
        }

      }, error = function(e) {
        showNotification(paste("Error in recluster analysis:", e$message), type = "error", duration = 10)
      })
    })

    # ORA Plot Output
    output$ora_plot <- renderPlot({
      if (is.null(values$ora_result)) {
        plot.new()
        text(0.5, 0.5, "Run ORA analysis to see results", cex = 1.5, col = "gray")
        return()
      }



      tryCatch({
        egt_plot_results(values$ora_result, ntop = 20, P.adj = 1)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error creating plot:", e$message), cex = 1, col = "red")
      })
    })

    # ORA Table Output
    output$ora_table <- DT::renderDataTable({
      if (is.null(values$ora_result)) {
        return(data.frame(Message = "Run ORA analysis to see results"))
      }

      # Format results for display
      display_data <- values$ora_result
      if ("p.adjust" %in% colnames(display_data)) {
        display_data$p.adjust <- format(display_data$p.adjust, scientific = TRUE, digits = 3)
      }
      if ("pvalue" %in% colnames(display_data)) {
        display_data$pvalue <- format(display_data$pvalue, scientific = TRUE, digits = 3)
      }

      display_data
    }, options = list(scrollX = TRUE, pageLength = 10))

    # GSEA Plot Output
    output$gsea_plot <- renderPlot({
      if (is.null(values$gsea_result)) {
        plot.new()
        text(0.5, 0.5, "Run GSEA analysis to see results", cex = 1.5, col = "gray")
        return()
      }



      tryCatch({
        egt_plot_results(values$gsea_result, ntop = 20, P.adj = 1)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error creating plot:", e$message), cex = 1, col = "red")
      })
    })

    # GSEA Table Output
    output$gsea_table <- DT::renderDataTable({
      if (is.null(values$gsea_result)) {
        return(data.frame(Message = "Run GSEA analysis to see results"))
      }

      # Format results for display
      display_data <- values$gsea_result
      if ("p.adjust" %in% colnames(display_data)) {
        display_data$p.adjust <- format(display_data$p.adjust, scientific = TRUE, digits = 3)
      }
      if ("pvalue" %in% colnames(display_data)) {
        display_data$pvalue <- format(display_data$pvalue, scientific = TRUE, digits = 3)
      }
      if ("NES" %in% colnames(display_data)) {
        display_data$NES <- round(display_data$NES, 3)
      }

      display_data
    }, options = list(scrollX = TRUE, pageLength = 10))

    # Recluster Table Output
    output$recluster_table <- gt::render_gt({
      if (is.null(values$recluster_result)) {
        # Return a simple gt table with message when no results
        data.frame(Message = "Run recluster analysis to see results") |>
          gt::gt() |>
          gt::tab_style(
            style = gt::cell_text(align = "center", size = "18px", color = "gray"),
            locations = gt::cells_body()
          ) |>
          gt::cols_label(Message = "")
      } else {
        tryCatch({
          # Extract gt table from EnrichGT object
          values$recluster_result@gt_object
        }, error = function(e) {
          # Return error message as gt table
          data.frame(Error = paste("Error displaying table:", e$message)) |>
            gt::gt() |>
            gt::tab_style(
              style = gt::cell_text(color = "red"),
              locations = gt::cells_body()
            ) |>
            gt::cols_label(Error = "")
        })
      }
    })

    # LLM Summary Output
     output$recluster_llm_summary <- renderUI({
       if (is.null(values$llm_summary)) {
         return(div(style = "color: #6c757d; font-style: italic; text-align: center; margin-top: 20px;", 
                    "No LLM summary available. Check 'Generate LLM Summary' option and ensure LLM object is provided."))
       }

       if (is.character(values$llm_summary)) {
         return(div(style = "color: #dc3545; margin: 15px;", values$llm_summary))
       }

       # If it's an EnrichGT object with LLM annotation, format the summary
       if (inherits(values$llm_summary, "EnrichGT_obj") && !is.null(values$llm_summary@LLM_Annotation)) {
         cluster_names <- values$llm_summary@LLM_Annotation@pathways$cluster_names
         if (length(cluster_names) == 0) {
           return(div(style = "color: #6c757d; font-style: italic; text-align: center; margin-top: 20px;", 
                      "No clusters found for LLM summary."))
         }

         # Create HTML content
         html_content <- div()
         
         for (i in seq_along(cluster_names)) {
           cluster_name <- cluster_names[i]
           cluster_div <- div(style = "margin-bottom: 25px; padding: 15px; background-color: white; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);")
           
           # Get cluster title
           title_idx <- which(values$llm_summary@LLM_Annotation@genes_and_title$clustersName == cluster_name)
           if (length(title_idx) > 0) {
             title <- values$llm_summary@LLM_Annotation@genes_and_title$resultsTitle[[title_idx]]
             cluster_div <- tagAppendChild(cluster_div, 
               h4(style = "color: #495057; margin-bottom: 15px; border-bottom: 2px solid #007bff; padding-bottom: 8px;", 
                  HTML(paste0(cluster_name, " - ", title))))
           } else {
             cluster_div <- tagAppendChild(cluster_div, 
               h4(style = "color: #495057; margin-bottom: 15px; border-bottom: 2px solid #007bff; padding-bottom: 8px;", 
                  cluster_name))
           }
           
           # Get pathway summary
           pathway_idx <- which(values$llm_summary@LLM_Annotation@pathways$cluster_names == cluster_name)
           if (length(pathway_idx) > 0) {
             pathway_summary <- values$llm_summary@LLM_Annotation@pathways$results[[pathway_idx]]
             cluster_div <- tagAppendChild(cluster_div,
               div(style = "margin-bottom: 15px;",
                 h5(style = "color: #28a745; margin-bottom: 8px;", "Pathway Summary:"),
                 div(style = "color: #495057; line-height: 1.6; white-space: pre-wrap; font-family: inherit;", 
                     HTML(gsub("\n", "<br>", pathway_summary)))
               ))
           }

           # Get gene summary
           gene_idx <- which(values$llm_summary@LLM_Annotation@genes_and_title$clustersName == cluster_name)
           if (length(gene_idx) > 0) {
             gene_summary <- values$llm_summary@LLM_Annotation@genes_and_title$results[[gene_idx]]
             cluster_div <- tagAppendChild(cluster_div,
               div(style = "margin-bottom: 10px;",
                 h5(style = "color: #17a2b8; margin-bottom: 8px;", "Gene Analysis:"),
                 div(style = "color: #495057; line-height: 1.6; white-space: pre-wrap; font-family: inherit;", 
                     HTML(gsub("\n", "<br>", gene_summary)))
               ))
           }
           
           html_content <- tagAppendChild(html_content, cluster_div)
         }
         
         return(html_content)
       }

       return(div(style = "color: #dc3545; margin: 15px;", "Unable to format LLM summary."))
     })


  }

  # Run the application
  shinyApp(ui = ui, server = server, options = list(port = port, host = host, launch.browser = launch.browser))
}
