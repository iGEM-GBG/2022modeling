# setwd("~/Documents/GitHub/2022modeling/WebApp")
library('shiny')
library('shinyWidgets')

# for modeling
library('plotly')
library('pracma')
library('misc3d')
source('WLCModelFunctions.R')

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
  
  titlePanel(h1("PLOP - Polymer Length Optimization Program", align="center")),
  
  sidebarLayout(
    sidebarPanel(id="sidebar",
      img(src = "iGEM_logga_vit.png", width = "100%", align="center"),
      
      uiOutput("model_div"),
      
      uiOutput("rerun_button")
    ),
      
    
    mainPanel(
      hr(),
      
      actionButton("show_instructions", label = "Show/hide usage instructions", align="center"),
      uiOutput("usage_instructions"),
      
      
      hr(),
      
      h1("Model results", align="center"),
      uiOutput("selected_parameters"),
      
      uiOutput("modeling_results")
      
    )
  )
)



# Define server logic ----
server <- function(input, output, session) {
  # Initialize parameter values
  min_length_value <- reactiveVal(5)
  max_length_value <- reactiveVal(10)
  segment_count_value <- reactiveVal(30)
  length_count_value <- reactiveVal(2)
  radius1_value <- reactiveVal(1)
  radius2_value <- reactiveVal(1)
  persistence_length_value <- reactiveVal(1)
  aT_value <- reactiveVal(10)
  
  
  output$selected_parameters <- renderUI({
    #h2("You have selected the following parametres")
  })
  
  output$modeling_parameters <- renderUI({
    if (input$polymer_count == "1 Polymer") {
      tagList(
        # Choose polymer length range
        numericInput("min_length", min = 0, label = p("Enter the shortest polymer length in nm"), value = min_length_value()),
        numericInput("max_length", min = 0, label = p("Enter the longest polymer length in nm"), value = max_length_value()),

        hr(style = "border: 2px solid white;"),
        numericInput("segment_count", label = p("Enter the number of segments for the polymer"), value = segment_count_value()),
        
        hr(style = "border: 2px solid white;
           "),
        p("Below you chose how many lengths you want to model within the chosen range. 
          If 1 is chosen, the polymer will have length equal to the lowest bound of the range."),
        numericInput("length_count", min=1, label = p("Enter the number of lengths to model"), value = length_count_value()),
      
        hr(style = "border: 2px solid white;"),
        p("Below you chose how large the proteins are at the ends of the polymer. 
          The proteins are modelled as spheres. If radius 0 is chosen there are no proteins at the ends"),
        numericInput("radius1", min=0.01, label = p("Enter the radius of the protein at the start of the polymer in nm"), value = radius1_value()),
        numericInput("radius2", min=0.01, label = p("Enter the radius of the protein at the end of the polymer in nm"), value = radius2_value()),
        
        hr(style = "border: 2px solid white;"),
        p("Below you chose the persistence length of the polymer, which is a measure of polymer stiffness.
          It normally varies between 0.45 and 0.6 nm for protein linkers."),
        numericInput("persistence_length", min=0, label = p("Enter persistence length of the polymer in nm"), value = persistence_length_value()),
        
        hr(style = "border: 2px solid white;"),
        p("Below you chose how many simulations you want to run. 
          Depending on desired runtime and other parameters, the appropriate number of simulations could range from 10^4 to 10^6."),
        numericInput("aT", min=1, label = p("Enter the number of simulations for each polymer length"), value = aT_value()),
        
        hr(style = "border: 2px solid white;"),
        actionButton("run_one_polymer", label = "Run PLOP"),
      )
    }
    else {
      tagList(
        p("The option to model two polymers is under development"),
        
        #textInput("some_text2", label = h3("Text input"), value = "Enter text..."),
        
        actionButton("run_two_polymers", label = "Run PLOP")
      )
    }
  })
  
  observeEvent(input$run_one_polymer, {
    # TODO: validate that parameters allowed
    # Save parameter values
    min_length_value(input$min_length)
    max_length_value(input$max_length)
    segment_count_value(input$segment_count)
    length_count_value(input$length_count)
    radius1_value(input$radius1)
    radius2_value(input$radius2)
    persistence_length_value(input$persistence_length)
    aT_value(input$aT)
    
    onePolymerResult <- onePolymer(range = c(input$min_length, input$max_length), 
                                   N = input$segment_count,
                                   length_count = input$length_count, 
                                   radiuses = c(input$radius1, input$radius2),
                                   aT = input$aT,
                                   lp = input$persistence_length)
    
    
    
    lapply(1:input$length_count, function(i){ # Linker simulation outputs
      output[[paste("simulation",i)]] <- renderPlotly({ onePolymerResult[[i]]$example_chain })
      output[[paste("histogram",i)]] <- renderPlotly({ onePolymerResult[[i]]$distance_histogram })
    })
  
    output$modeling_results <- renderUI({
      tagList(
        
        h3("Contact probabilities"),
        lapply(1:input$length_count, function(i){
          p(paste("Linker length ", probs(input$length_count, onePolymerResult)$lengths[i], "nm: ", probs(input$length_count, onePolymerResult)$probs[i]))
        }),
       
        renderPlotly(
          plot_ly() %>%
            add_trace(x = probs(input$length_count, onePolymerResult)$lengths,
                      y = probs(input$length_count, onePolymerResult)$probs,
                      line = list(color = "blue"),
                      type='scatter', mode = 'lines')%>%
            layout(title = "Protein contact probability plot", paper_bgcolor='black',
                   plot_bgcolor = "black",
                   font = list(color = '#FFFFFF'),
                   xaxis=list(title=list(text ='Polymer length (nm)', font = list(color = "white")),gridcolor = 'white'),
                   yaxis= list(title=list(text ='Contact probability', font = list(color = "white")), gridcolor = 'white'))
        ),
        h3("Example simulations for all polymer lengths"),
        lapply(1:input$length_count, function(i){
            plotlyOutput(paste("simulation",i), width='80%')
        }),
        
        h3("Estimations of probability density functions for surface distances"),
        lapply(1:input$length_count, function(i){
          plotlyOutput(paste("histogram",i), width='80%')
          
        }),
        
        
      )
    })
    
    run()
  })
  
  
  observeEvent(input$run_two_polymers, {
    
    output$modeling_results <- renderUI({
      # TODO: show all modeling parameters
      
      h3(twoPolymers(input$some_text2)$text)
    })
    
    run()
  })
  probs <- function(n, result){ # Get contact probability plot
    contact_probs = 1:n
    for (i in 1:n){
      contact_probs[i] <- result[[i]]$contact_prob
    }
    lengths <- linspace(input$min_length, input$max_length, input$length_count)
    return (list("probs" = contact_probs, "lengths" = lengths))
  }
  run <- function(){ # Remove input choices and add button for running again
    output$model_div <- renderUI({  })
    output$rerun_button <- renderUI({actionBttn("run_again", label = "Delete results to run again") })
  }
  
  observeEvent(input$run_again, {
    output$model_div <- renderUI({ model_div })
    removeUI(selector = '#run_again')
    output$modeling_results <- renderUI({  })
  })
  
  shows_instructions <- reactiveVal(0)
  observeEvent(input$show_instructions, {
    if (shows_instructions() == 0){
      output$usage_instructions <- renderUI({ usage_instructions })
      shows_instructions(1)
    } else {
      output$usage_instructions <- renderUI({ })
      shows_instructions(0)
    }
  })
  
  output$model_div <- renderUI({ model_div })
  
  model_div <- div(id="model_div",
      h1("What do you want to model?", align="center"),
      
      tags$style("#polymer_count {background-color: black; color: white;}"),
      selectInput(inputId="polymer_count",
                  label=h3("Model 1 or 2 polymers"),
                  choices = c("1 Polymer", "2 Polymers"),
                  selected = c("1 Polymer"),
                  selectize = FALSE),
      
      hr(style = "border: 2px solid white;"),
      
      uiOutput("modeling_parameters"),
  )
  
  usage_instructions <- div(id="usage_instructions",
                            h1("Usage instructions", align="center"),
                            p("The easiest way to get a general understanding of what PLOP can do is 
        to press the Run PLOP button at the bottom left of this page. 
        PLOP calculate probabilities that proteins come in contact, 
        depending on polumer length of polymers connected to the proteins.
        This enables polymer length optimization. A polymer can for example be a protein linker 
        or a DNA strand."),
                            p("Note that depending on the chosen parameters, 
                              it can take a long time to run PLOP. You can choose to run PLOP
                              for a couple of minutes and get very questionable accuracy,
                              or run for hours and get very high accuracy."),
                            
                            h2("Research"),
                            p("Some things need to be known about a system system before running PLOP. 
        The following steps describe what research is needed."),
                            p("1. Research stiffness of polymers. The parameter for stiffness is called persistence length, which is needed to run PLOP.
        For example a protein linker’s persistence length is dependent on the choice of amino acids."),
                            p("2. Research size of proteins connected to your polymers."),
                            p("3. Optional: look at similar systems for a first guess of optimal linker length."),
                            
                            h2("Advice"),
                            p("There are also other things to take into consideration when running PLOP. 
        Here are some advice regarding each of the parameters PLOP needs."),
                            p("1. Range (shortest and longest polymer lengths): Choose a wide range of polymer lengths. After each run of PLOP, it is recommended to reduce the range according to the new information."),
                            p("2. Number of segments: Choose a low number of segments (20-30) for approximate results. Choose a high number of segments (50-100) for more exact results. The number of segments needed is dependent on the system."),
                            p("3. Polymer lengths: If possible, choose a low number (1-5) of polymer lengths to model. This is in order to reduce the time each run of PLOP takes. It is possible to run PLOP with many polymer lengths, but be prepared to wait for hours if it is high."),
                            p("4. Protein radiuses: Enter average radius for proteins. Approximations are often adequate. For modeling polymers not connected to proteins, set protein radiuses to zero."),
                            p("5. Persistence length: A lower persistence length results in the polymer bending more. PLOP is not made for optimizing persistence length, but it is possible to use it for testing how persistence length affects the system. This can be used for informing the decision of amino acid sequence in protein linkers."),
                            p("6. Number of simulations: Choose a low number of simulations (10^3-10^4) for approximate results. Choose a high number of simulations (10^5 - 10^6) for more exact results."),
                            
                            
                            p("For a more details and eplanations about PLOP, visit ",
                              a("our webpage", 
                                href = "https://2022.igem.wiki/chalmers-gothenburg/software")),
                            p("If you have questions regarding usage or results from the web app, feel free to send an email to anestrandalvin@gmail.com.
        We would also be grateful for any feedback, especially if you have empirical data from studies where this tool has been used."),
                            
  )
}

shinyApp(ui = ui, server = server)


