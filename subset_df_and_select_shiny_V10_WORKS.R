library(shiny)
library(gdata)
library(dplyr)
library(tidyr)
library(plotly)
library(pheatmap)
#install.packages("janitor")
library(janitor)

##### HEY HEY READ ME #####
#  IVE BEEN TESTING THIS WITH TEST_NUCLEAR, NEED TO CHANGE THIS

# Define UI for app
ui <- fluidPage(
  #subset by patient
  selectInput("patients", label = "Select patients to merge:",
              choices = list("Patient X" = "patientX", "Patient Y" = "patientY", "Patient Z" = "patientZ"),
              multiple = FALSE,
              selected = "patientx"),
  #subset by nuc/AA
  selectInput("nuc_aa", "Select nucleotide or amino acid data:",
              choices = list("Nucleotide" = "nuclear", "Amino Acid" = "amino_acid"),
              multiple = FALSE,
              selected = "nuclear"),
  
  selectInput("effect", "Select Effect:",
              choices = unique(finalDF$Effect),
              selected = unique(finalDF$Effect),
              multiple = TRUE),
  
  selectInput("Gene", "Select Gene:",
              choices = unique(test_nuclear_use$Gene),
              selected = unique(test_nuclear_use$Gene),
              multiple = TRUE),
  
  # Define the select input for selecting the region value
  selectInput("region", "Select Region:",
              choices = unique(test_nuclear_use$region),
              selected = unique(test_nuclear_use$region),
              multiple = TRUE),
  

    plotlyOutput("heatmap")
  
  
)

# Define server logic
server <- function(input, output) {
  
  
  # Define reactive function for merging selected patients
  merged_data <- reactive({
    
    # Filter data based on selected patients and data type
    selected_patients <- input$patients
    selected_data <- input$nuc_aa
    
    # CAN FILTER THIS, but not important
    filtered_data <- switch(selected_data,
                            "nuclear" = test_nuclear_use,
                            "amino_acid" = test_nuclear_use)
    
    filtered_data %>% filter(Gene %in% input$Gene & Effect %in% input$effect & region %in% input$region)
    
    filtered_data <- filtered_data[filtered_data$patient %in% selected_patients,]
    
    # just cleaning up a bit, can get rid of the select column later when figure out how to use metadata
    if (selected_data == "nuclear") {
      filtered_data <- filtered_data %>%
        select(Nucleotide, everything()) %>%
        select(-c(Gene, POS, REF, ALT, Effect, Protein, aa_ref,	Position,
                  aa_mut, region, nmax, percentage
        ))
    } else {
      filtered_data <- filtered_data %>%
        select(Protein, everything()) %>%
        select(-c(Gene, POS, REF, ALT, Effect, Nucleotide, aa_ref,	Position,
                  aa_mut, region, nmax, percentage
        ))
    }
    
    
    #%>%
    # select(-c(Gene, POS, REF, ALT, Effect, Protein,
    #                                patient, time,	Position,	region,
    #                                percentage))
    
    
    # Split data by time and patient
    split_data <- split(filtered_data, list(filtered_data$time, filtered_data$patient))
    
    
    # Use lapply to rename columns in each split dataframe
    split_data <- lapply(split_data, function(df) {
      # Get the unique values in the "time" column
      unique_times <- unique(df$time)
      # Get the maximum value in the "time" column
      max_time <- max(unique_times)
      # Get the index of the "time" column
      time_idx <- which(names(df) == "time")
      # Rename the "time" column to "time_max" where max is the maximum value found in the column
      names(df)[time_idx] <- paste0("time_", max_time)
      
      # Get the index of the "n" column
      n_idx <- which(names(df) == "n")
      # Rename the "n" column to "n_max" where max is the maximum value found in the "time" column
      names(df)[n_idx] <- paste0("n_", max_time)
      #subset Nucleotide and protein names to
      if (selected_data == "nuclear") {
        nucleotide_idx <- which(names(df) == "Nucleotide")
        names(df)[nucleotide_idx] <- paste0("Nucleotide_", max_time)
      } else {
        protein_idx <- which(names(df) == "Protein")
        names(df)[protein_idx] <- paste0("Protein_", max_time)  
      }
      
      
      # nucleotide_idx <- grep("^Nucleotide", names(df))
      # protein_idx <- grep("^Protein", names(df))
      # 
      # Return the renamed dataframe
      return(df)
    })
    
    # cbind back together with cbindx for diff rows
    merged_data <- do.call("cbindX", split_data)
    
    if (selected_data == "nuclear") {
      df_merged <- merged_data %>%
        unite("Nucleotide", starts_with("Nucleotide"), sep = "_", na.rm = TRUE) %>%
        mutate(Nucleotide = sub(".*_", "", Nucleotide))
      
      merged_data <- df_merged 
    } else {
      df_merged <- merged_data %>%
        unite("Protein", starts_with("Protein"), sep = "_", na.rm = TRUE) %>%
        mutate(Protein = sub(".*_", "", Protein))
    } 
    
    df_merged <- df_merged %>%
      select(-c(patient, starts_with("time")))
    
    merged_data <- df_merged %>% remove_rownames %>% column_to_rownames(var="Nucleotide")
    
    # df_merged_trans <- transpose(df_merged)
    # df_merged_trans <- df_merged_trans %>% row_to_names(row_number = 1)
    # rownames(df_merged_trans) <- colnames(df_merged)
    
    
    #merged_data <- df2
    
    
    
    
    # Return merged dataframe
    return(merged_data)
  })
  
  heatmap_data <- reactive({
    finalDF <- merged_data()
    
    # # Use lapply to rename columns in each split dataframe
    # split_data <- lapply(split_data, function(df) {
    #   # Get the unique values in the "time" column
    #   unique_times <- unique(df$time)
    #   # Get the maximum value in the "time" column
    #   max_time <- max(unique_times)
    #   # Get the index of the "
    #   finalDF <- merged_data %>% select(matches("_\\d+$"))
    
    # Replace NA with 0
    finalDF[is.na(finalDF)] <- 0
    
    # Convert data to matrix
    heatmap_matrix <- data.matrix(finalDF)
    
    # Create heatmap
    p <- plot_ly(x = colnames(heatmap_matrix), y = rownames(heatmap_matrix),
                 z = heatmap_matrix, type = "heatmap")
    
    # Convert plot to plotly object
    return(p)
  }) 
  
  output$heatmap <- renderPlotly({
    heatmap_data()
  })
  
}


# Run the app

shinyApp(ui = ui, server = server)
