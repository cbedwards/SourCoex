library(shiny)
# library(devtools)
# install_github("cbedwards/SourCoex",
               # auth_token="133a2ab1fa7cf7a6d1d4e547e9204b27fa66c13d")
library(SourCoex)

#### Functions For Now: Definitely easier when my package is available to servers.




# Define UI for app that draws a histogram ----
ui <- fluidPage(

  # App title ----
  titlePanel("Visualizing microbe models"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      selectInput(
        inputId = "specid",
        label = "Species (to compare data to our curves)",
        choices = specMap$name.full

      ),
      selectInput(
        inputId = "modname",
        label = "Model",
        choices = c("Logistic Growth","One-Species Tilman")

      ),
      ## Logistic growth parameters
      conditionalPanel(
        condition = "input.modname=='Logistic Growth'",

        sliderInput(inputId = "r",
                    label = "r (growth rate)",
                    min = 0,
                    max = 2,
                    step=.01,
                    value = .2),

        # Input: Slider for the number of bins ----
        sliderInput(inputId = "k",
                    label = "k (carrying capacity in thousands of CFUs)",
                    min = 0,
                    max = 500,
                    value = 100)
      ),

      ## Tilman 1-species parameters
      conditionalPanel(
        condition = "input.modname=='One-Species Tilman'",

        sliderInput(inputId = "tilr",
                    label = "r (growth conversion factor)",
                    min = 0,
                    max = 1,
                    step=.0001,
                    value = .05),

        # Input: Slider for the number of bins ----
        sliderInput(inputId = "d",
                    label = "d (per capita growth rate)",
                    min = 0,
                    max = .5,
                    value = 0.01,
                    step=.0001),
        sliderInput(inputId = "R",
                    label = "R_0 (initial resource units per transfer)",
                    min = 1,
                    max = 200,
                    value = 10)
      ),
      # App description
      includeMarkdown("model-context.md"),
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: grah ----
      plotOutput(outputId = "trajecplot"),

      #Function descriptions
      conditionalPanel(
      condition = "input.modname=='Logistic Growth'",
      withMathJax(includeMarkdown("logistic.md"))
      ),
      conditionalPanel(
        condition = "input.modname=='One-Species Tilman'",
        withMathJax(includeMarkdown("tilman1.md"))
      )

      #text context




    )
  )
)


# Define server logic required to draw a histogram ----
server <- function(input, output) {

  # Histogram of the Old Faithful Geyser Data ----
  # with requested number of bins
  # This expression that generates a histogram is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$bins) change
  # 2. Its output type is a plot
  output$trajecplot <- renderPlot({
    #grab data of appropriate species
    specid=specMap$name.data[specMap$name.full==input$specid]
    dat.prepped=dataprep_landis_solo(specid)
    ## parameters etc for logistic, if appropriate
    if(input$modname=="Logistic Growth"){
      ode_cur=ode_log
      parms.cur=c(input$r, input$k)
      aug.ls=list()
      parmnames.cur=parmnamesLog
      units.cur=unitsLog
    }
    ## parameters etc for 1-species Tilman, if appropriate
    if(input$modname=="One-Species Tilman"){
      ode_cur=ode_Til_1spec
      parms.cur=c(input$tilr, input$d, input$R)
      aug.ls=list(Til=1)
      parmnames.cur=parmnamesTil1spec
      units.cur=unitsTil1spec
    }
    modfit=plotter_landis_solo(parms=parms.cur,
                               ode_fun=ode_cur,
                               parmnames = parmnames.cur,
                               parmunits = units.cur,
                               dat=dat.prepped$dat,
                               specid=specid,
                               noparms=TRUE,
                               aug.ls=aug.ls)
    modfit
  })

}

shinyApp(ui = ui, server = server)
