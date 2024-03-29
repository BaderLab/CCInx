#' Run the Shiny app to interactively explore cell-cell interactions.
#'
#' @param INX The list containing edgelist and node metadata as generated by
#'   \code{\link{BuildCCInx}}.
#' @param includeHeadHTML Default=NA. If you'd like an HTML script to be included
#'   the webpage <head> section (such as the Google Analytics tracking script,
#'   see https://shiny.rstudio.com/articles/google-analytics.html), pass the
#'   path to the script HTML file here.
#' @param ... Named options that should be passed to the
#'   \code{\link[shiny]{runApp}} call (these can be any of the following:
#'   "port", "launch.browser", "host", "quiet", "display.mode" and "test.mode").
#'
#' @export

ViewCCInx <- function(INX,
                      includeHeadHTML=NA,
                      ...) {

  if (!is.list(INX)) {
    stop("INX must be the output of the function CCInx::BuildCCInx()")
  }
  if (!all(names(INX) %in% c("nodes","edges"))) {
    stop("INX must be the output of the function CCInx::BuildCCInx()")
  }

  if (is.na(includeHeadHTML)) {
    includeHeadHTML <- system.file("blank.html",package="CCInx")
  }

  ui <- pageWithSidebar(
    titlePanel("CCInx Viewer"),
    sidebarPanel(
      tags$head(includeHTML(includeHeadHTML)),
      fluidRow(
        column(6,
               uiOutput("cellTypesA"),
               selectInput("proteinTypeA","Node Type A:",
                           choices=c("Receptor","Ligand","ECM"),selected="Ligand")
        ),
        column(6,
               uiOutput("cellTypesB"),
               selectInput("proteinTypeB","Node Type B:",
                           choices=c("Receptor","Ligand","ECM"),selected="Receptor")
        )
      ),
      hr(),
      uiOutput("FilterType"),
      uiOutput("Filter"),
      #        radioButtons("YSpacing","Y-axis:",inline=T,
      #                     choices=list(Relative="relative",
      #                                  Absolute="absolute"))

      hr(),
      fluidRow(
        column(8,
               downloadButton("PlotSave","Save Figure As:"),
               align="right"),
        column(4,
               selectInput("imageFileType",label=NULL,
                           choices=list(PDF="pdf",
                                        EPS="eps",
                                        TIFF="tiff",
                                        PNG="png"))
        )
      )
      #### TESTING ####
      # textOutput("TEST")
    ),
    mainPanel(
      plotOutput("CCInx",inline=T)
    )
  )

  server <- function(input,output,session) {
    output$cellTypesA <- renderUI({
      temp <- unique(INX$nodes$cellType)
      selectInput("cellTypeA","Cell Type A:",
                  choices=temp,selected=temp[1])
    })
    output$cellTypesB <- renderUI({
      temp <- unique(INX$nodes$cellType)
      selectInput("cellTypeB","Cell Type B:",
                  choices=temp,selected=temp[2])
    })

    output$FilterType <- renderUI({
      if (is.null(attr(INX,"GeneStatistic"))) {
        selectInput("WhichFilter","Filter network by:",
                    choices=list("All (could be slow!)"="All",
                                 "Expression magnitude"="Magn",
                                 "Top weighted edges"="Top",
                                 "Gene symbols"="Genes"),
                    selected="Top")
      } else {
        selectInput("WhichFilter","Filter network by:",
                    choices=list("All (could be slow!)"="All",
                                 "Differential expression significance"="Stat",
                                 "Differential expression magnitude"="Magn",
                                 "Top weighted edges"="Top",
                                 "Gene symbols"="Genes"),
                    selected="Top")
      }
    })
    output$Filter <- renderUI({
      switch(input$WhichFilter,
             "Stat"=numericInput("GeneStatistic",
                                 label=paste0("Maximum ",
                                              attr(temp_inx(),"GeneStatistic"),
                                              ":"),
                                 value=0.05),
             "Magn"=numericInput("GeneMagnitude",
                                 label=paste0("Minumum absolute ",
                                              attr(temp_inx(),"GeneMagnitude"),
                                              ":"),
                                 value=.1),
             "Top"=numericInput("TopN",
                                label=HTML(paste0(
                                  "Number of edges sorted by weight<br/>(absolute mean scaled ",
                                  attr(temp_inx(),"GeneMagnitude"),
                                  " of nodes):"
                                )),
                                value=20),
             "Genes"=selectInput("GeneNames",label="Search by gene symbols:",
                                 choices=temp_inx()$nodes$gene,multiple=T)
      )
    })

    temp_inx <- reactive({
      FilterInx_step1(INX,
                      cellTypeA=input$cellTypeA,
                      cellTypeB=input$cellTypeB,
                      proteinTypeA=input$proteinTypeA,
                      proteinTypeB=input$proteinTypeB)
    })

    temp_inx2 <- reactive({
      switch(input$WhichFilter,
             "All"=temp_inx(),
             "Stat"=FilterInx_GeneStatistic(temp_inx(),input$GeneStatistic),
             "Magn"=FilterInx_GeneMagnitude(temp_inx(),input$GeneMagnitude),
             "Top"=FilterInx_topN(temp_inx(),input$TopN),
             "Genes"=FilterInx_genenames(temp_inx(),
                                         unique(unlist(strsplit(input$GeneNames," |,"))))
      )
    })

    temp_figH <- function() {
      temp <- max(table(temp_inx2()$nodes$side)) *
        strheight("ABC",units="inches") *
        1.25 * 96
      if (temp <= 800) { temp <- 800 }
      return(temp)
    }

    output$CCInx <- renderPlot({
      DoPlotInx(temp_inx2()) #,input$YSpacing)
    },res=96,width=720,height=temp_figH)

    output$PlotSave <- downloadHandler(
      filename=reactive({
        paste0("CCInx_",
               input$cellTypeA,"_",input$proteinTypeA,"_",
               input$cellTypeB,"_",input$proteinTypeB,
               ".",input$imageFileType)
      }),
      content=function(file) {
        switch(input$imageFileType,
               "pdf"=grDevices::cairo_pdf(file,height=temp_figH()/96,width=7.5,fallback_resolution=300),
               "eps"=grDevices::cairo_ps(file,height=temp_figH()/96,width=7.5,fallback_resolution=300),
               "tiff"=grDevices::tiff(file,height=temp_figH()/96,width=7.5,units="in",res=300),
               "png"=grDevices::png(file,height=temp_figH()/96,width=7.5,units="in",res=300))
        DoPlotInx(temp_inx2())
        grDevices::dev.off()
      }
    )


    #### TESTING ####
    output$TEST <- renderPrint(paste0("CCInx_",
                                      input$cellTypeA,"_",input$proteinTypeA,"_",
                                      input$cellTypeB,"_",input$proteinTypeB,
                                      ".",input$imageFileType))
  }
  shinyApp(ui,server,options=list(...))
}
