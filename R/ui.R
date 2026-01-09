#' Generating layout for the graphical interface.
#'
#' @param params input parameters
#' @return shiny ui
#' @keywords internal
ui <- function(params){
  shiny::fluidPage(
    theme = shinythemes::shinytheme("simplex"),
    shiny::navbarPage(
      "pandemonium",

  shiny::tabPanel("Input",
  shiny::fluidPage(

  # Sidebar with a slider input for number of bins
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::selectInput(
        "linkage",
        "Linkage",
        choices = c("complete", "single",
                    "average", "median",
                    "centroid", "mcquitty",
                    "ward.D", "ward.D2"),
        selected = "complete"
      ),
      shiny::checkboxInput("useCov", "Include covariance"),
      shiny::selectInput("coord", "Coordinates",
                         choices = c("Pull", "p-val")),
      shiny::selectInput("metric", "Distance measure",
                         choices = c("euclidean", "maximum", "manhattan", "canberra",
                                     "binary", "minkowski")),
      shiny::selectInput("kC", "Number of clusters", choices = 2:8),
      shiny::selectInput("px", "x", choices = params, selected = params[1]),
      shiny::selectInput("py", "y", choices = params, selected = params[2]),
      shiny::uiOutput("conditions1"),
      shiny::uiOutput("conditions2")

    ),

    # Show a plot of the generated distribution
    shiny::mainPanel(
      shiny::splitLayout(cellWidths = c("50%", "50%"),
                         cellArgs = list(style = "height: 325px"),
                  shiny::imageOutput("heatmap"),
                  shiny::imageOutput("chi2")),
      shiny::splitLayout(cellWidths = c("50%", "50%"),
                         cellArgs = list(style = "height: 325px"),
                  shiny::imageOutput("wc"),
                  shiny::imageOutput("sigbins"))
    )
  ))),
  shiny::tabPanel("Benchmarks",
                  shiny::fluidPage(
                    DT::DTOutput('benchmarks')
                  )),
  shiny::tabPanel("Distance breakdown",
                  shiny::fluidPage(
                    shiny::textOutput("distText"),
                    shiny::plotOutput("hist")
                    )),
  shiny::tabPanel("Coordinates",
                  shiny::fluidPage(
                    shiny::plotOutput("scalar"),
                    shiny::plotOutput("pc"),
                    shiny::plotOutput("pc2")
                  )),
  shiny::tabPanel("Dimension Reduction",
                  shiny::fluidPage(
                    shiny::splitLayout(cellWidths = c("30%", "30%", "30%"),
                    shiny::plotOutput("tsne"),
                    shiny::plotOutput("umap"),
                    shiny::plotOutput("lle"))
                    )

                  ),
  shiny::tabPanel("Tour display",
                  shiny::fluidPage(
                    shiny::fluidRow(shiny::column(6,
                    shiny::sliderInput(
                      "n", "", min = 1, max = 100, value = 1,
                      animate = shiny::animationOptions(interval = 100)
                    ),
                    shiny::imageOutput("tourImg", width = "100%")),
                    shiny::column(6,
                    shiny::sliderInput(
                      "npca", "", min = 1, max = 100, value = 1,
                      animate = shiny::animationOptions(interval = 100)
                    ),
                    shiny::imageOutput("tourImgPca", width = "100%")
                  )))),

  shiny::tabPanel("Comparison",
                  shiny::fluidPage(

                    # Sidebar with a slider input for number of bins
                    shiny::fluidRow(
                      shiny::column(2,
                        shiny::selectInput(
                          "linkageA",
                          "Linkage A",
                          choices = c("complete", "single",
                                      "average", "median",
                                      "centroid", "mcquitty",
                                      "ward.D", "ward.D2"),
                          selected = "complete"
                        ),
                        shiny::checkboxInput("useCovA", "Include covariance A"),

                        shiny::selectInput("coordA", "Coordinates A",
                                           choices = c("Pull", "p-val")),
                        shiny::selectInput("metricA", "Distance measure A",
                                           choices = c("euclidean", "maximum", "manhattan",
                                                       "canberra", "binary", "minkowski")),
                        shiny::selectInput("kA", "Number of clusters A", choices = 2:8)
                        ),
                      shiny::column(2,
                        shiny::selectInput(
                          "linkageB",
                          "Linkage B",
                          choices = c("complete", "single",
                                      "average", "median",
                                      "centroid", "mcquitty",
                                      "ward.D", "ward.D2"),
                          selected = "complete"
                        ),
                        shiny::checkboxInput("useCovB", "Include covariance B"),
                        shiny::selectInput("coordB", "Coordinates B",
                                           choices = c("Pull", "p-val")),
                        shiny::selectInput("metricB", "Distance measure B",
                                           choices = c("euclidean", "maximum", "manhattan",
                                                       "canberra", "binary", "minkowski")),
                        shiny::selectInput("kB", "Number of clusters B", choices = 2:8)
                      ),

                      # Show a plot of the generated distribution
                      shiny::column(8,
                        shiny::splitLayout(cellWidths = c("50%", "50%"),
                                           cellArgs = list(style = "height: 325px"),
                        shiny::plotOutput("heatmapA"),
                        shiny::plotOutput("heatmapB"))
                      )
                    ),
                    shiny::fluidRow(
                      shiny::column(4, shiny::plotOutput("tableAB")),
                      shiny::column(8, shiny::splitLayout(
                        cellWidths = c("50%", "50%"),
                        cellArgs = list(style = "height: 325px"),
                        shiny::plotOutput("wcA"),
                        shiny::plotOutput("wcB")))
                    ))),
  shiny::tabPanel("Statistics",
                  shiny::fluidPage(shiny::plotOutput("clusterstats"))),
  shiny::tabPanel("Detourr",
                  shiny::fluidPage(
                      shiny::fluidRow(shiny::column(6,
                                                  detourr::displayScatter2dOutput("detour1",width = "100%", height = "400px")),          
                                      shiny::column(6,
                                                  detourr::displayScatter2dOutput("detour2",width = "100%", height = "400px")))
                      ))
))}
