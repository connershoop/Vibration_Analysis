library(shiny)
library(ggplot2)
library(DT)
library(stats)
library(shinyWidgets)

server <- function(input, output, session){
  library(stats)
  library(pracma)
  
  ####################################################################################
  ## plots data of the class 1d matrix
  ## fs is sampling frequency
  ## freq Min and freq Max specify frequency bounds on plot
  ## log specifies whether or not to apply a log scale to the y axis
  ## numAv specifies the number of averages to give, when this increases, the distance between plot points in the frequency domain will increase
  ## machineFreq scales the frequency domain so as to make 1 the machine frequency
  ## grid puts a vertical line on every integer in the frequency axis
  ## displayPeaks shows peaks on the graph
  ## peaksQuantile sets the quantile threshhold for the peaks that would be excluded
  
  ## Function to plot spectrum
  plotSpectrum <- function(data, fs, freqMin = 0, freqMax = 100, log = FALSE, numAv = 1, machineFreq = 1,
                           grid = FALSE, displayPeaks = FALSE, peaksQuantile =0){
    #* This takes input data and depending on value for numAv, breaks the data into numAv segments, performs fft's on them,
    # and then averages the y values to a final value, essentialy performing TSA
    
    #initialize
    x.spec <- spectrum(data[0:floor((length(data)/numAv)*1-1),], plot=FALSE) 
    spy <-  2*x.spec$spec
    spx <- x.spec$freq*fs/machineFreq
    
    #run
    if (numAv ==1){}
    else{
      for (i in 2:numAv){
        x.spec <- spectrum(data[floor((length(data)/numAv)*(i-1)):floor((length(data)/numAv)*i-1),], plot=FALSE) 
        
        spy <-  spy + 2*x.spec$spec
      }
      spy = spy/numAv
    }
    
    
    if (log == TRUE){
      plot(spy~spx,xlab="frequency",ylab="spectral density",type="l", xlim = c(freqMin,freqMax), log = "y")
    }else{
      plot(spy~spx,xlab="frequency",ylab="spectral density",type="l", xlim = c(freqMin,freqMax))
    }
    #Add vertical grid
    if (grid == TRUE){
      abline(v = freqMin:freqMax,  lty = 2, col = "grey")
    }
    
    if (displayPeaks == TRUE){
      peaks <- findpeaks(spy, nups = 1, threshold = peaksQuantile*max(spy))
      points(x = peaks[,2]*fs/length(data)/machineFreq*numAv, y = peaks[, 1], pch = 23, col ="red")
    }
  }
  ####################################################################################
  ####################################################################################
  library(wavelets)
  library(reshape2)
  library(ggplot2)
  library(RColorBrewer)
  ## this function takes a time series of data and turns it into a wavelet plot
  ## data is a 1D matrix of time series data
  ## fs is the sampling frequency of the data
  ## time is the length of time wished to be evaluated
  ## wavelet is the type of wavelet to use during the analysis
  
  plotWaveletCoef <- function(data=data, fs=102400L, time=1L, wavelet="d8"){
    
    dwtData <- dwt(data[1:(2*time*fs)],filter = wavelet, boundary = "periodic")
    
    # New plot method for dwtData
    
    ## we are given a stacked sequence of arrays, with smaller and smaller values
    ## the below code fills a matrix, waveletMatrix, full to be the same number of rows
    ## as there are arrays, but each to have the length of the longest row, retaining
    ## much of the original values, evenly spaced.  This prepares for the data to be plotted to a 
    ## heat map
    
    waveletMatrix <- matrix(, ncol = length(dwtData@W[["W1"]]), nrow = length(dwtData@W))
    waveletMatrix[1,] <- dwtData@W[["W1"]]
    for (i in 2:length(dwtData@W)){
      length <- length(dwtData@W[[paste("W",i, sep = "")]])
      finalLength <-  length(dwtData@W[["W1"]])
      seg = floor(finalLength/length)
      for (j in 1:length){
        waveletMatrix[i, ((j-1)*seg+1):(j*seg)] <- dwtData@W[[paste("W",i, sep = "")]][j]
      }
      
    }
    
    
    meltedWaveletMatrix <- melt(waveletMatrix)
    meltedWaveletMatrix[is.na(meltedWaveletMatrix)] <- 0
    ## add negatives to equal 
    meltedWaveletMatrix$value <-  meltedWaveletMatrix$value + abs(min(meltedWaveletMatrix$value))
    meltedWaveletMatrix[meltedWaveletMatrix==0] <- min(meltedWaveletMatrix$value[meltedWaveletMatrix$value!=0])
    ## prepare for graphing
    meltedWaveletMatrix$Var2 <- meltedWaveletMatrix[,2]/fs
    
    my_breaks <-  seq(min(meltedWaveletMatrix$value), max(meltedWaveletMatrix$value),max(meltedWaveletMatrix$value)/5 )
    
    ggplot(data = meltedWaveletMatrix, aes(y=Var1, x=Var2, fill=value)) + 
      geom_tile() + scale_fill_gradient(name = "value", my_breaks, trans = "log") + 
      xlab("time(s)") + ylab("Wavelet size")+ scale_fill_gradientn(colors=c("black"
                                                                            ,"darkred","red","yellow","white"))
  }
  ####################################################################################
  
  
  
  myData <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(NULL)
    data <- read.csv(inFile$datapath, header = TRUE)
    data
  })
  
  output$contents <- DT::renderDataTable({
    DT::datatable(myData())       
  })
  
  
  fftPlot <- eventReactive(input$updateFft, {
    fftData <- as.matrix(myData())
    plotSpectrum(data=fftData, fs=input$fs, freqMin=input$freqMin, freqMax=input$freqMax, log = input$log,
                 numAv=input$numAv, machineFreq=input$machineFreq, grid = input$grid, displayPeaks = input$peaks,
                 peaksQuantile=input$peaksQuantile)
    
  })
  
  waveletPlot <- eventReactive(input$updateWavelet, {
    waveletData <- as.matrix(myData())
    plotWaveletCoef(data=waveletData, fs=input$fs, time=input$time,
                    wavelet=input$waveletType)
    
  })
    
  output$tsPlot <- renderPlot({
    timeData <- as.matrix(myData())
    plot.ts(y = timeData, x = seq(0, (length(timeData)-1))/input$fs, type = "l", xlab = "time(s)",
            ylab = "vibration amplitude")
  })
  
  output$fftPlot <- renderPlot({
    fftPlot()
  })
  
  output$waveletPlot <- renderPlot({
    waveletPlot()
  })

  
}

write.csv(data.frame(a = 1:10, b = letters[1:10]), 'test.csv')



ui<- shinyUI(fluidPage(
  titlePanel("Vibration Analysis"),
  sidebarLayout(
    sidebarPanel(
      tabsetPanel(
        tabPanel("data",
      fileInput('file1', 'Choose CSV File',
                accept=c('text/csv',
                         'text/comma-separated-values,text/plain',
                         '.csv'))
        ),
      tabPanel("spectrum",
               actionButton("updateFft", "Update Spectrum Graph"),
               switchInput("log", "Log Scale?"),
               switchInput("grid", "Vertical Grid?"),
               switchInput("peaks", "Display Peaks?"),
               numericInput("fs", "sampling frequency", 102400),
               numericInput("freqMin", "minimum display frequency", 0),
               numericInput("freqMax", "maximum display frequency", 100),
               numericInput("numAv", "number of averages",1),
               numericInput("machineFreq", "operating frequency of machine to scale by", 1),
               numericInput("peaksQuantile", "Quantile for peaks", 0.001)
    ),
    tabPanel("wavelet", 
             actionButton("updateWavelet", "Update Wavelet Spectrogram"),
             numericInput("fs", "sampling frequency", 102400),
             numericInput("time", "signal time to use (seconds)", 1),
             textInput("waveletType", "type of wavelet,
                        Daubechies
                        2,4,6,8,10,12,14,16,18,20.
                       
                        Least Asymetric
                        8,10,12,14,16,18,20.
                       
                        Best Localized
                        14,18,20.
                       
                        Coiflet
                        6,12,18,24,30.", "d8")
             ))),
    mainPanel(
      tabsetPanel(
      tabPanel("data", DT::dataTableOutput('contents')), 
      tabPanel("spectrum",plotOutput('fftPlot', height = "800px")),
      tabPanel("wavelet",plotOutput('waveletPlot', height = "800px"))
      ) 
    )
  ),fluidPage("Time Vibration Data", plotOutput('tsPlot'))
)
)

shinyApp(ui,server)