library(shiny)
library(ggplot2)
library(plotly)
# may want to add in addition to the pymol plot, a little blurb about what Nucleotide change caused the AA mut
# Define the user interface
ui <- fluidPage(
  
  # Define a slider input for selecting the range of n values to include
  sliderInput("n_range", "Select n range:",
              min = 1, max = max(finalDF$n), value = c(1, max(finalDF$n))),
  
  # Define a slider input for zooming based on the x-axis
  sliderInput("zoom_range", "Zoom",
              min = 1, max = max(as.numeric(finalDF$Position)), value = c(1, max(as.numeric(finalDF$Position)))),
  
  # Define a checkbox for turning on/off the labels
  checkboxInput("show_labels", "Show Labels", TRUE),
  
  checkboxInput("show_histogram", "Show Histogram", FALSE),
  
  # Define the select input for selecting the effect value
  selectInput("effect", "Select Effect:",
              choices = unique(finalDF$Effect),
              selected = unique(finalDF$Effect),
              multiple = TRUE),
  
  # Define the select input for selecting the region value
  selectInput("region", "Select Region:",
              choices = unique(finalDF$region),
              selected = unique(finalDF$region),
              multiple = TRUE),
  
  # Define the plot output
  plotlyOutput("plot")
)

# Define the server function
server <- function(input, output) {
  
  # Define the reactive subset of the data frame based on the selected n range
  subset_data <- reactive({
    finalDF %>% filter(n >= input$n_range[1] & n <= input$n_range[2]) %>%
      filter(Effect %in% input$effect & region %in% input$region)
  })
  
  # Define the plot to render based on the selected n range and checkbox input
  output$plot <- renderPlotly({
    
    # Subset the data frame based on the selected n range
    data <- subset_data()
    
    # Define the lollipop plot

    if (input$show_histogram) {
      
    p1 <- ggplot(data, aes(x = n)) + 
      geom_bar(stat = "count")
    plot(p1)
    
    } else {
    if (input$show_labels) {
      p <- ggplot(data, aes(x = as.numeric(Position), y = n)) +
        geom_segment(aes(xend = as.numeric(Position), yend = 0), color = "gray50") +
        geom_point(aes(color = region, shape = Effect), size = 3,
                   position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0, dodge.width = 0.5)) +
        scale_color_manual(values = dd.col) +
        geom_text(aes(label = Protein), nudge_y = 1.5, angle = 45, size = 4, check_overlap = TRUE) +
        scale_shape_manual(values = c(15, 17, 19)) +
        labs(x = "Position", y = "Count") +
        theme_classic() +
        theme(plot.margin = unit(c(1,1,2,1), "lines")) +
        scale_x_continuous(n.breaks = 12)

    } else {
      p <- ggplot(data, aes(x = as.numeric(Position), y = n)) +
        geom_segment(aes(xend = as.numeric(Position), yend = 0), color = "gray50") +
        geom_point(aes(color = region, shape = Effect), size = 3,
                   position = position_jitterdodge(jitter.width = 0.1, jitter.height = 0, dodge.width = 0.5)) +
        scale_color_manual(values = dd.col) +
        scale_shape_manual(values = c(15, 17, 19)) +
        labs(x = "Position", y = "Count") +
        theme_classic() +
        theme(plot.margin = unit(c(1,1,2,1), "lines")) +
        scale_x_continuous(n.breaks = 12)
    }

    p <- ggplotly(p, textangle = 45, check_overlap = TRUE)

    # Add ability to zoom and scale the x-axis using Plotly
    p <- p %>% layout(xaxis = list(range = input$zoom_range))
    

    # Return the plotly object
    return(p)
    return(p1)
    }
  })
  
}

# Run the app
shinyApp(ui = ui, server = server)    
