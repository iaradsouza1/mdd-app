library(shiny)
library(shinythemes)
library(dplyr)
library(ggplot2)
library(DT)
library(ggrepel)

# Load tx differential expression table 
load("data/plot_tables.rda")

# Get the gene list
dte_genes <- unique(plot_table$hgnc_symbol)
g <- unique(isa_df$group)

# Define UI  --------------------------------------------------------------
ui <- navbarPage(
  theme = shinytheme("sandstone"),
  "Transcript expression in MDD",
  
  tabPanel("DTE",
   sidebarLayout(
      sidebarPanel(
        selectInput("gene", "Select a gene from the table", choices = dte_genes),
        tags$h5("Description:"),
        textOutput(outputId = "description_dte")
      ),
      mainPanel(
        plotOutput("plot")
      )
    ),
    dataTableOutput("table")
  ),
  
  tabPanel("DTU",
    sidebarLayout(
     
      sidebarPanel(
        selectInput("groups", "Select a group", choices = g, selected = "aINS_female", multiple = F),
        selectInput("id_dtu", "Select a gene (name + ensG)", choices = unique(isa_df$id)),
        tags$h5("Description:"),
        textOutput(outputId = "description_dtu")
      ),
     
      mainPanel(
        htmlOutput("file")
      )
    )
  )
)

# Define server -----------------------------------------------------------
server <- function(input, output, session) {
  
  # DTE server ----
  # Dataset reactive
  dataset <- reactive(plot_table)
  
  # Organize plot_table to display
  df_show <- reactive({
    dataset() %>% 
      mutate(group = paste(region, gender, sep = "_")) %>% 
      dplyr::select(group, gene, hgnc_symbol, entrezgene_id, tx, logFC_gene, 
                    logFC_tx, padj_gene, padj_tx, prop, canonical)
    
  })
  
  # Filter gene to plot
  plot <- reactive({
    dataset() %>%
      filter(hgnc_symbol %in% input$gene)
  })
  
  # Outputs
  output$description_dte <- renderText({
    desc_list[input$gene]
  })
  
  output$table <- renderDataTable(df_show(), 
                                  options = list(searchHighlight = TRUE,
                                                 scrollX = T), 
                                  filter = 'top')
  output$plot <- renderPlot({
    ggplot(plot()) +
      facet_wrap(region ~ gender, drop = T) +
      geom_hline(aes(yintercept = logFC_gene), lty = 2, color = "#0000b1ff") +
      geom_hline(yintercept = 0, color = "#00000a8a") +
      geom_text(aes(x = 0.2, y = logFC_gene+.1, label = "gene logFC"), color = "#0000b1ff") + 
      scale_x_continuous(limits = c(0, 1), breaks = seq(0,1, 0.1)) +
      geom_point(aes(x = 0.5, y = logFC_tx, col = signif, size = prop)) + 
      scale_y_continuous(name = "logFC") +
      scale_size_continuous(name = "Transcript proportion") + 
      labs(title = paste0("DTE analysis for ", input$gene, " gene"), x = "", y = "logFC") + 
      geom_label_repel(aes(x = 0.5, y = logFC_tx, label = tx, color = signif),
                       direction = "both",
                       segment.colour = "grey20",
                       nudge_x = 0.3,
                       segment.size = 0.3, box.padding = 0.5, show.legend = F) +
      scale_color_manual(name = "Transcript significance (padj <= 0.05)", values = c("S" = "red", "NS" = alpha("black", 0.7)),
                         labels = c("Not significant", "Significant")) +
      theme_bw() + 
      theme(axis.ticks.x = element_blank(), 
            axis.text.x = element_blank())
  })
  
  # DTU server ----
  
  # Create dynamic selectInput for genes in each group
  observeEvent(input$groups, {
    x <- input$groups
    y <- isa_df$id[isa_df$group == x]
    updateSelectInput(session, inputId = "id_dtu", label = "Select a gene (name + ensG)",
                      choices = y,
                      selected = tail(y, 1))
  })
  
  # List all files in www/
  files <- reactive({
    list.files("www", full.names = T)
  })
  
  # Get the gene 
  gene_ens <- reactive({
    gene <- as.character(unlist(strsplit(input$id_dtu, split = "_")))[1]
    ens <- as.character(unlist(strsplit(input$id_dtu, split = "_")))[2]
    return(c(gene, ens))
  })
  
  # Get the file corresponding to chosen gene
  graph <- reactive({
    x <- files()[( grepl(gene_ens()[1], files()) & grepl(gene_ens()[2], files()) ) & grepl(input$groups, files())]
    x <- gsub("www/", "", x)
    return(x)
  })
  
  # Outputs
  output$description_dtu <- renderText({
    desc_list[gene_ens()[1]]
  })
  
  observeEvent(input$id_dtu, {
    output$file <- renderUI({
      tags$iframe(src = graph(), style = "height:800px; width:100%;scrolling=yes")
    })
  })
  
  observeEvent(input$id_dtu, {
    print(isa_df$entrezgene_id[match(gene_ens()[1], isa_df$gene_name)])
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
