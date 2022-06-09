#app.R
library(shiny)
source("functions/functions.R")
ui <- fluidPage(
  tabPanel("Transcription", fluid = TRUE,
           sidebarLayout(
             sidebarPanel(id = "sidebar",
                          selectInput(
                            inputId = "morph1",
                            label = "Compare...",
                            choices = c("Pachon","Molino","Tinaja","Rascon","Rio Choy")
                          ),
                          selectInput(
                            inputId = "morph2",
                            label = "to...",
                            choices = c("")
                          ),
                          radioButtons("trstat",
                                       label = textOutput("trstat_label"),
                                       choices = c("logFC", "p")
                          ),
                          radioButtons("direction",
                                       label = textOutput("dir_label"),
                                       choices = c("Above", "Below")
                          ),
                          
                          # If morph2 IS set to control, populate a drop-down of available conditions
                          conditionalPanel(
                            condition = "input.morph2 == 'Control'",
                            selectInput(
                              inputId = "condition",
                              label = "Condition?",
                              choices = c("Sleep deprivation")
                            )
                          ),
                          
                          conditionalPanel(
                            condition = "input.trstat == 'logFC'",
                            uiOutput("FCthresh_updater")
                          ),
                          conditionalPanel(
                            condition = "input.trstat == 'p'",
                            sliderInput(inputId = "pthresh",
                                        label = textOutput("pthresh_label"),
                                        min = 0, max = 1, value = 0.5, step = 0.000001)
                            
                          ),
                          actionButton("Transc_enter","Find Genes"),
                          # Change sorting of output table
                          radioButtons("TRwhich_sort",
                                       label = "Sort by...",
                                       choices = c("logFC", "p")
                          ),
                          conditionalPanel(
                            condition = "input.TRwhich_sort == 'p'",
                            radioButtons("TRsort_p",
                                         label = "...from...",
                                         choices = c("High to Low" = "TRp_hl", 
                                                     "Low to High" = "TRp_lh")
                            )
                          ),
                          conditionalPanel(
                            condition = "input.TRwhich_sort == 'logFC'",
                            radioButtons("TRsort_FC",
                                         label = "...from...",
                                         choices = c("High to Low" = "TRFC_hl", 
                                                     "Low to High" = "TRFC_lh")
                            )
                          ),
                          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                           tags$div("Crawling through the data...",id="loadmessage"))
             ),
             mainPanel(fluidRow(
               conditionalPanel(condition = "input.Transc_enter",
                                downloadButton("TranscDL", "Download", class = "download"),
                                tableOutput("transc_table_out"),
                                textOutput("transc_wrnings_out")
               )
             )
             )
           )
  ),
  
)

server <- function(input, output, session) {
  #empty data frame for control.
  condition_control <- data.frame()
  
  observe({
    
    transc_morph_choices <- c("Control", "Rio Choy")
    # Transcription Page: Update morph-selection widget to only enable
    # comparisons between current morph and morph which is NOT morph1
    updateSelectInput(session = getDefaultReactiveDomain(),
                      "morph2",
                      choices = transc_morph_choices[transc_morph_choices != input$morph1],
                      selected = tail(transc_morph_choices, 1)
    )
    # Transcription Page: Update min and max values for logFC based on data
    output$FCthresh_updater <- renderUI({
      if((input$trstat == "logFC") & (input$morph2 == "Control")){
        sec.table <- condition_control$logFC[
          (paste(c(input$morph1, "-", input$morph2), collapse = "")
           %in% condition_control$Comparison) & 
            (condition_control$Condition == input$condition)]
        sliderInput("FCthresh", "...the following threshold:", 
                    min = round(min(sec.table), 2), 
                    max = round(max(sec.table), 2), value = 0, step = 0.0005)
        
      }else if((input$trstat == "logFC") & (input$morph2 != "Control")){
        sec.table <- morph1.morph2$logFC[
          grepl(input$morph1, morph1.morph2$Comparison) &
            grepl(input$morph2, morph1.morph2$Comparison)
        ]
        sliderInput("FCthresh", "...the following threshold:", 
                    min = round(min(sec.table), 2), 
                    max = round(max(sec.table), 2), value = 0, step = 0.0005)
        
      }else{
        sliderInput(inputId = "FCthresh",
                    label = "Threshold selector will appear here once you have inputted enough data to determine min and max possible logFC values",
                    min = 0, max = 0, value = 0)
      }
    })
  })
  
  # Transcription Page labels
  output$trstat_label <- renderText("Display genes whose logFC/p-value...")
  output$dir_label <- renderText("... is above/below....")
  output$pthresh_label <- renderText("... the following threshold:")
  
  # Transcription Page: Output a table with specified transcription data
  transc_table <- eventReactive(input$Transc_enter, valueExpr = {
    if(input$morph2 == "Control"){
      condit <- input$condition
    }else{
      condit <- "Between morph"
    }
    if(input$trstat == "logFC"){
      TranscTable(morph1 = input$morph1,
                  morph2 = input$morph2,
                  condition = condit,
                  direction = input$direction,
                  tr.stat = input$trstat,
                  tr.thresh = input$FCthresh,
                  GOTable = GeneToGO)
      
    }else if(input$trstat == "p"){
      TranscTable(morph1 = input$morph1,
                  morph2 = input$morph2,
                  condition = condit,
                  direction = input$direction,
                  tr.stat = input$trstat,
                  tr.thresh = input$pthresh,
                  GOTable = GeneToGO)
    }
  })
  output$transc_table_out <- renderTable({
    # Change log fold changes to have 5 decimal points
    reformattedTranscT <- data.frame(
      transc_table()[[1]][,1:4],
      format(transc_table()[[1]][,5], digits = 5),
      format(transc_table()[[1]][,6], digits = 5),
      transc_table()[[1]][,7:10]
    )
    names(reformattedTranscT) <- names(transc_table()[[1]])
    # Sort table based on p or logFC value
    if(input$TRwhich_sort == "p"){
      if(input$TRsort_p == "TRp_hl"){
        reformattedTranscT <- reformattedTranscT[order(
          as.numeric(reformattedTranscT$`p-value`), decreasing = T),]
      }else{
        reformattedTranscT <- reformattedTranscT[order(
          as.numeric(reformattedTranscT$`p-value`)),]
      }
    }else if(input$TRwhich_sort == "logFC"){
      if(input$TRsort_FC == "TRFC_hl"){
        reformattedTranscT <- reformattedTranscT[order(reformattedTranscT$logFC, 
                                                       decreasing = T),]
      }else{
        reformattedTranscT <- reformattedTranscT[order(reformattedTranscT$logFC),]
      }
    }
    
    reformattedTranscT
  })
  output$transc_wrnings_out <- renderText({transc_table()[[2]]})
  
  # Transcription Page: Enable downloading of Transcription table
  output$TranscDL <- downloadHandler(
    filename = function() {
      paste("CaveCrawler-Transcription-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(transc_table(), file, row.names = F)
    }
  )
    
}

shinyApp(ui = ui, server = server)