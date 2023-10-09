#link backend code
source('global.R')

ui <- dashboardPage(
  dashboardHeader(
    title = tags$span(
      style = "font-family: 'Roboto', sans-serif; font-size: 24px; font-weight: Thin;",
      "scRNA seq Analysis"
    )
  ), 
  dashboardSidebar( #this block creates tabs when the tab scRNAseq is selected then you have an option to upload an rds file and hit buttons, run or reset. 
    sidebarMenu(
      id = 'tab',
      useShinyjs(),
      menuItem("Home Page", tabName = "home", icon = icon("list")),
      menuItem("scRNAseq Analyzer", tabName = "input", icon = icon("edit")),
      conditionalPanel(
        condition = "input.tab == 'input'", 
        div(
          fileInput("file", "Upload File", multiple = FALSE, accept = c('.rds')),
          actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
          actionButton("run", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%")
        )
      )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(
        tabName = "input",
        tabsetPanel(
          id = 'main_tabs', 
          tabPanel("Instructions", includeMarkdown('./markdown/instructions.md')) #this links a Markdown file for the writing in the instructions tab
        )
      ),
      tabItem(
        tabName = "home", 
        tags$h1(
          HTML("<u>Welcome to The scRNAseq Seurat Analysis RShiny App</u>"),
          style = "font-family: 'Roboto Thin', serif; font-size: 30px; font-weight: Thin;"
        )
      )
    )
  )
)

server <- function(input, output, session) { 
  options(shiny.maxRequestSize = 300 * 1024^2) #this is the maximum size for the input file. 

  shinyjs::disable("run") #its disabled until a file is inputted

  observe({
    if (!is.null(input$file)) { #if a file is inputed, run is enabled. 
      shinyjs::enable('run')
      
    } else {
      shinyjs::disable("run")
    }
  })

  observeEvent(input$reset, {  #file is the upload file call recognized so that is reset and run is again disabled when the reset button is hit. 
    shinyjs::reset("file")
    shinyjs::disable("run")
    removeTab('main_tabs', 'UMAP') #the tabs UMAP and Gene Expression disappear. 
    removeTab('main_tabs', 'Gene Expression')
  })

  observeEvent(input$run, {  #when the run button is hit run is then disabled and a loading icon appears and we can load the input file's path
    shinyjs::disable("run")

    show_modal_spinner(text = "Prepping Plots ...")
    obj <- load_seurat_obj(input$file$datapath) #since the error comments are vectors anything that is not a vector will display the metadata UMAP plot and the plot of each gene.
    if (is.vector(obj)) {
      showModal(modalDialog(
        title = 'Error with file',
        HTML(paste("<h5>There is an error with the file you uploaded. See below for more details. </h5> <br>",
          paste(unlist(obj), collapse = "<br><br>"))
        )
      ))
      shinyjs::enable("run")
    } else {
      output$umap <- renderPlot({
        create_metadata_UMAP(obj, input$metadata_col)
      })

      output$featurePlot <- renderPlot({
        create_feature_plot(obj, input$gene)
      })

      output$download_umap <- downloadHandler( #this allows us to download the file with the following extensions and input
        filename = function(){
            paste0(input$metadata_col, '_UMAP.png')
        },
        content = function(file){
            plot <- create_metadata_UMAP(obj, input$metadata_col)
            ggsave(filename = file, width = 10, heigh =5, type = "cairo")
        }
      )

      output$downloadFeaturePlot <- downloadHandler(
        filename = function(){
            paste0(input$gene, '_feature_plot.png')
        },
        content = function(file){
            plot <- create_feature_plot(obj, input$gene)
            ggsave(filename = file, width = 10, heigh =5, type = "cairo")
        }
      )

      insertTab(  #these plots only get inputed if a seurat object is uploaded with no error messages.
        inputId = 'main_tabs',
        tabPanel(
          "UMAP",
          fluidRow(
            column(
              width = 8,
              plotOutput(outputId = 'umap'),
              downloadButton('download_umap', 'Download UMAP')
            ),
            column(
              width = 4,
              selectizeInput(
                "metadata_col",
                "Metadata Column",
                choices = colnames(obj@meta.data), #this displays the dropdown.

              )
            )
          )
        )
      )
      insertTab(
        inputId = 'main_tabs',
        tabPanel(
          "Gene Expression",
          fluidRow(
            column(
              width = 8,
              plotOutput(outputId = 'featurePlot'),
              downloadButton('downloadFeaturePlot', 'Download Feature Plot')
            ),
            column(
              width = 4,
              selectizeInput(
                "gene",
                "Genes",
                rownames(obj)
              )
            )
          )
        )
      )
      remove_modal_spinner() #if the above commands execute then the loading icon is removed. 
      shinyjs::enable("run")
    }
  })
}

shinyApp(ui, server)


