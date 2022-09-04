# setwd("~/Documents/GitHub/2022modeling/WebApp")
library('shiny')
library('shinyWidgets')

# for modeling
library('plotly')
library('pracma')
source('WLCModelFunctions.R')

# for testing
library(maps)
library(mapproj)
source("helpers.R")
counties <- readRDS("data/counties.rds")

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  setBackgroundColor("black"),
  tags$head(tags$style(
    HTML('
         #sidebar {
            background-color: black;
         }

         * { 
           font-family: "Arial";
           color: white;
         }
         p {
           font-size: 20px;
         }
         '
         )
  )),
  
  titlePanel(h1("WLC Polymer Length Optimization", align="center")),
  
  sidebarLayout(
    sidebarPanel(id="sidebar",
      img(src = "iGEM_logga_vit.png", width = "100%", align="center"),
      h1("What do you want to model?", align="center"),

      tags$style("#polymer_count {background-color: black; color: white;}"),
      selectInput(inputId="polymer_count",
                 label=h3("Model 1 or 2 polymers"),
                 choices = c("1 Polymer", "2 Polymers"),
                 selected = c("1 Polymer"),
                 selectize = FALSE),
      
      hr(),
      
      uiOutput("modeling_parameters"),
      ),
    
    mainPanel(
      hr(),
      
      h1("Usage instructions", align="center"),
      
      p("For a more detailed description, visit ",
        a("our webpage", 
          href = "https://2022.igem.wiki/chalmers-gothenburg/software")),
      p("If you have questions regarding usage or results from the web app, feel free to send an email to anestrandalvin@gmail.com.
        We would also be grateful for any feedback, especially if you have empirical data from studies where this tool has been used."),
      
      hr(),
      
      h1("Model results", align="center"),
      h3(textOutput("selected_polymer_count")),
      
      uiOutput("modeling_results"),
      
    )
  )
)

# Define server logic ----
server <- function(input, output, session) {
  output$selected_polymer_count <- renderText({
    paste("You have selected",input$polymer_count)
  })
  output$modeling_parameters <- renderUI({
    if (input$polymer_count == "1 Polymer") {
      tagList(
        textInput("some_text1", label = h3("Text input"), value = "Enter text..."),
      
        actionButton("run_one_polymer", label = "Run simulations with one polymer")
      )
    }
    else {
      tagList(
        textInput("some_text2", label = h3("Text input"), value = "Enter text..."),
        
        actionButton("run_two_polymers", label = "Run simulations with two polymers")
      )
    }
  })
  
  reactive_sim_value <- reactive({
    a = switch(input$polymer_count,
           "1 Polymer" = onePolymer(onePolymerTextParam),
           "2 Polymers" = twoPolymers(textOutput("some_text2")),
           )
  })
  # reactive_sim_value <- reactiveValues(a = reactive_simulation)
  
  doing_setting <- FALSE # not setting 
  observeEvent(input$run_one_polymer, {
 
    # Update modeling parameters
    output$some_text1 <- renderPrint({ input$some_text1 })
    onePolymerTextParam <- textOutput("some_text1")
    #reactiveTextParam <- reactiveValues(a = textOutput("some_text1"))
    
    output$modeling_results <- renderUI({
      h3(onePolymer(onePolymerTextParam)$a)
    })
    
    # try freezing
    # session$onFlushed(function() freezeReactiveValue(reactive_sim_value, "a"))
    # freezeReactiveValue(reactive_sim_value, "a")
    
  })
 
  
  observeEvent(input$run_two_polymers, {
    # Update modeling parameters
    output$some_text2 <- renderPrint({ input$some_text2 })
    
    output$modeling_results <- renderUI({
        h3(reactive_sim_value()$text)
    })
  })
}

shinyApp(ui = ui, server = server)


