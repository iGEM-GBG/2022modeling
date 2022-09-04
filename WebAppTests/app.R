# setwd("~/Documents/GitHub/2022modeling/WebAppTests")
library('shiny')
library('shinyWidgets')

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
  titlePanel(h1("WLC Chain Length Optimization", align="center")),
  
  sidebarLayout(
    sidebarPanel(id="sidebar",
                 img(src = "iGEM_logga_vit.png", width = "100%", align="center"),
                 h1("What do you want to model?", align="center"),
                 #checkboxGroupInput(inputId ="polymer_count", label="Model 1 or 2 polymers",
                 #  choices = list("1 polymer", "2 polymers")
                 #),
                 
                 tags$style("#polymer_count {background-color: black; color: white;}"),
                 selectInput(inputId="polymer_count",
                             label=h3("Model 1 or 2 polymers"),
                             choices = c("1 Polymer", "2 Polymers"),
                             selected = c("1 Polymer"),
                             selectize = FALSE),
                 sliderInput("slider1", label = h3("Slider"), min = 0, 
                             max = 100, value = 50),
                 column(4, verbatimTextOutput("value")),
                 
                 textInput("name", "What's your name?"),
                 textOutput("greeting"),
    ),
    
    mainPanel(
      
      h1("Usage instructions", align="center"),
      
      p("For a more detailed description, visit ",
        a("our webpage", 
          href = "https://2022.igem.wiki/chalmers-gothenburg/software")),
      p("If you have questions regarding usage or results from the web app, feel free to send an email to anestrandalvin@gmail.com.
        We would also be grateful for any feedback, especially if you have empirical data from studies where this tool has been used."),
      
      p(textOutput("selected_polymer_count")),
      
    )
  )
)

# Define server logic ----
server <- function(input, output, session) {
  output$greeting <- renderText({
    paste0("Hello ", input$name, "!")
  })
}

shinyApp(ui = ui, server = server)


