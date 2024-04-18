library(shiny)
library(reshape2)
source("helper.R")
library(ggplot2)
rm(list=ls())

ui <- fluidPage(
  # Title of the APP: Sugarcane EST-DigitalNorthern
  titlePanel(title=div(img(src="LabBCES.png"), "RNA-Seq data Leaf Development again Genome SP80-3280 B2 - LabBCES")),
  fluidRow(
    p("Expression profiles from sugarcane using genome annotation of assembly B2 and RNA-Seq data from leaf development (", a(href='https://pubmed.ncbi.nlm.nih.gov/26714767/', 
    'Mattiello et al., 2015'),").",
      style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
    
    br(),
    p("We have two alternative annotations for the time being, one generated with GALBA and another with BRAKER3. The anotation with GALBA uses as extrisic evidence, for gene prediction, proteins sequences. While, BRAKER# uses both RNA-Seq data and protein sequences."
    ,style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
    
    br(),
    p("The goal of this tool is to aid you in visualizing the expression profile of a selected gene of interest.
    Users can control this with the slider ",strong('Percent of libraries with data:'),"."
      ,style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
    
    br(),
  ),
  sidebarLayout(
    # Sidebar for the first part of the analyses. Select the amount of samples to be taken into account to consider a gene as expressed.
    sidebarPanel(
      selectInput(
        inputId="selectedDataset",
        label = "Select dataset, either GALBA or BRAKER3",
        choices = c('GALBA','BRAKER3'),
        selected = NULL,
        multiple = FALSE,
        selectize = TRUE,
        width = NULL,
        size = NULL
      ),
    ),
    mainPanel(
      textOutput("selectedAnnotationReference"),
      plotOutput(outputId = "boxPlotExpressedGenes")
    )
  ),
  br(),
  p("The last figure can display the expression profile of each of the genes/transcripts kepth after the filters."
    ,style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
  
  br(),
  sidebarLayout(
    sidebarPanel(
      #htmlOutput("listSelectedGenesUI"),
      htmlOutput("listSelectedGenesUI2"),
    ),
    mainPanel(
      plotOutput(outputId = "expressionProfileSelectedGene"),
      tableOutput(outputId = "expressionTableSelectedGene")
    )
  ),
  p("Finaly, the user can download the whole dataset (all genes/transcripts), or the set after filtering.z"
    ,style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
  
  br(),
  sidebarLayout(
    sidebarPanel(
      selectInput("dataset", "Choose a dataset:",
                  choices = c("full", "filtered")),
      downloadButton("downloadData", "Download")
    ),
    mainPanel(
      tableOutput("table")
    )
  ),
  hr(),
  h1("Funding"),
  img(src="FAPESP.png"),
  img(src="CNPq.png"),
  img(src="RCG2I.png")
  
  
)

server <- function(input, output) {
  dataset <- reactive({
    fileNameTPM=paste("dataTPM_",input$selectedDataset,"_dups.rds",sep="")
    fileNameGeneIDs=paste("geneIDs_",input$selectedDataset,"_dups.rds",sep="")
    urlSourceTPM=paste("https://raw.githubusercontent.com/labbces/B2quant/main/",fileNameTPM,sep="")
    urlSourceGeneIDs=paste("https://raw.githubusercontent.com/labbces/B2quant/main/",fileNameGeneIDs,sep="")
    download.file(urlSourceTPM, fileNameTPM)
    download.file(urlSourceGeneIDs, fileNameGeneIDs)
    dataTPM<-readRDS(fileNameTPM)
    geneIDs<-readRDS(fileNameGeneIDs)
    return(list(dataTPM=dataTPM,geneIDs=geneIDs))
  })
  
  #Boxplot of genes expression values per EST library BRAKER3_dups.rds
  # dataset1<-reactive({
  #   input$selectedDataset
  # })
  # fileName=paste("dataTPM_",dataset1(),"_dups.rds")
  # urlSource=paste("https://raw.githubusercontent.com/labbces/B2quant/main/",fileName)
  # download.file(urlSource)
  # 
  output$boxPlotExpressedGenes <- renderPlot({
    # fileName=paste("dataTPM_",input$selectedDataset,"_dups.rds",sep="")
    # urlSource=paste("https://raw.githubusercontent.com/labbces/B2quant/main/",fileName,sep="")
    # download.file(urlSource, fileName)
    # dataTPM<-readRDS(fileName)
    dataTPM<-dataset()$dataTPM
    ggplot(dataTPM, aes(x=Sample,y=TPM,fill=DevStage)) +
      theme_bw()+
      theme(legend.position = 'none',axis.text.x = element_text(angle = 45))+
      geom_boxplot()+
      scale_y_log10()+
      ylab("TPM - Log10")
    
  })
  output$selectedAnnotationReference <- renderText({
    paste("You have selected: ",input$selectedDataset)
  })
  output$listSelectedGenesUI2 <-renderUI({
    #fileName=paste("dataTPM_",input$selectedDataset,"_dups.rds",sep="")
    #dataTPM<-readRDS(fileName)
    geneIDs<-dataset()$geneIDs
    selectInput(
    "selectedGene",p("Select the gene of interest, based on your previous filters:",
                     style="color:black; text-align:center"),
    choices=geneIDs[1:100]
  )
  })
  output$expressionProfileSelectedGene<-renderPlot({
    dataTPM<-dataset()$dataTPM
    ggplot(dataTPM[which(dataTPM$Gene == input$selectedGene),],aes(x=Sample, y=TPM, fill=DevStage))+
      theme_bw() +
      theme(legend.position = 'none',axis.text.x = element_text(angle = 45)) +
      geom_col() +
      ylab("TPM - Log10")+
      scale_y_log10() +
      ggtitle(input$selectedGene)
  })
  output$expressionTableSelectedGene<-renderTable({
    dataTPM<-dataset()$dataTPM
    dataTPM[which(dataTPM$Gene == input$selectedGene),]
  })
  datasetDownload <- reactive({
    switch(input$dataset,
           "full" = dataCPM,
           "filtered" = getExpressionProfile(dataCPM, input$minCPM,input$selectLib, input$percentLibs, libraries)
           )
  })
  output$table <- renderText({
    paste("Ready to download your selected dataset: ", input$dataset, "with ", nrow(datasetDownload()), "genes" )
  })
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$dataset, "_sugarcaneEST_CPM.csv", sep = "")
    },
    content = function(file) {
      write.csv(datasetDownload(), file, row.names = TRUE)
      }
  )
}

shinyApp(ui = ui, server = server)
