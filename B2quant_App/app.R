library(shiny)
library(ggplot2)
rm(list=ls())

ui <- fluidPage(
  # Title of the APP: Sugarcane EST-DigitalNorthern
  titlePanel(title=div(img(src="LabBCES.png"), "RNA-Seq data Leaf Development against Genome SP80-3280 B2 - LabBCES")),
  fluidRow(
    p("Expression profiles from sugarcane using genome annotation of assembly B2 and RNA-Seq data from leaf development (", a(href='https://pubmed.ncbi.nlm.nih.gov/26714767/', 
    'Mattiello et al., 2015'),"). In that publication we measure gene expression using RNA-Seq in four leaf segments of sugarcane SP80-3280. B0: Base Zero, B: Base, M: Middle part, P: Tip of leaf",
      style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
    
    p("The B2 genome assembly of cultivar SP80-3280 is a draft assembly generated with PacBio HiFi and Illumina HiC reads. It has a high N50.",
      style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
    
    br(),
    p("We have two alternative annotations for the time being using the B2 assembly. One generated with GALBA and another with BRAKER3. The anotation with GALBA uses as extrinsic evidence, for gene prediction, proteins sequences. BRAKER3 uses both RNA-Seq data and protein sequences."
    ,style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
    
    br(),
    p("The goal of this tool is to aid you in visualizing the expression profile of a selected (knonw) gene of interest. 
    You could select a gene of interest based on similarity searches using our ", a(href='http://labbces.cena.usp.br:4567/','Blast server'),".
      In the BLAST server the identifiers of the annotated transcriptional units and their deduced proteins have the following form XXXXX_gNNNNNN.tY, where XXXXX could either be BRAKER3 or GALBA,
      NNNNNN and Y are numbers. But please note that the expression values available here are for the gene, and not the transcriptional units. So in order to convert the identifiers from the transcript to the gene, you just need to remove the suffix .tN.
      for example for the transcript: BRAKER3_g152127.t2, the correspoding gene name that you can use here would be BRAKER3_g152127. Please also check that you select the proper dataset for the given gene identifier." 
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
        selected = 'BRAKER3',
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
  p("Type the identifier of your gene of interest."
    ,style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
  
  # br(),
  # sidebarLayout(
  #   sidebarPanel(
  #     #htmlOutput("listSelectedGenesUI"),
  #     htmlOutput("listSelectedGenesUI2"),
  #   ),
  #   mainPanel(
  #     plotOutput(outputId = "expressionProfileSelectedGene"),
  #     tableOutput(outputId = "expressionTableSelectedGene")
  #   )
  # ),
  br(),
  sidebarLayout(
    sidebarPanel(
      textInput(inputId='selectedGene',
                label='Gene of Interest',
                value = NULL,
                width = NULL,
                placeholder = "Type Gene ID")
    ),
    mainPanel(
      plotOutput(outputId = "expressionProfileSelectedGene"),
      tableOutput(outputId = "expressionTableSelectedGene")
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
    # fileNameGeneIDs=paste("geneIDs_",input$selectedDataset,"_dups.rds",sep="")
    urlSourceTPM=paste("https://raw.githubusercontent.com/labbces/B2quant/main/",fileNameTPM,sep="")
    # urlSourceGeneIDs=paste("https://raw.githubusercontent.com/labbces/B2quant/main/",fileNameGeneIDs,sep="")
    download.file(urlSourceTPM, fileNameTPM)
    # download.file(urlSourceGeneIDs, fileNameGeneIDs)
    dataTPM<-readRDS(fileNameTPM)
    # geneIDs<-readRDS(fileNameGeneIDs)
    # return(list(dataTPM=dataTPM,geneIDs=geneIDs))
    return(list(dataTPM=dataTPM))
  })

  output$boxPlotExpressedGenes <- renderPlot({
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
  output$expressionProfileSelectedGene<-renderPlot({
    dataTPM<-dataset()$dataTPM
    ggplot(dataTPM[which(dataTPM$Gene == input$selectedGene),],aes(x=Sample, y=TPM, fill=DevStage))+
      theme_bw() +
      theme(legend.position = 'none',axis.text.x = element_text(angle = 45)) +
      geom_col() +
      ylab("TPM - Log10")+
      ggtitle(input$selectedGene)
  })
  output$expressionTableSelectedGene<-renderTable({
    dataTPM<-dataset()$dataTPM
    dataTPM[which(dataTPM$Gene == input$selectedGene),]
  })
}

shinyApp(ui = ui, server = server)
