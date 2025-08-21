# Load required packages
library(ggplot2)  
library(readxl)    
library(tidyr)     
library(scales)    

# Read data from Excel file
data <- read_excel("Cancer_sample_TAI.xlsx")

# Prepare plotting data: convert from wide format to long format
melted_data <- pivot_longer(data, cols = everything(), names_to = "stage", values_to = "sampleTAI")

# Define custom colors for each stage
colors <- c("#D693BE", "#8EC8ED", "#F5B3A5", "#AED594")

# Rename x-axis labels
stage_labels <- c("StageI","StageII","StageIII","StageIV" )

# Set Arial font for Windows system
windowsFonts(Arial = windowsFont("Arial"))

# Create the boxplot with jitter points
ggplot(melted_data, aes(x = stage, y = sampleTAI, color = stage)) +
        geom_boxplot(fill = NA) +  # Boxplot without fill color
        geom_jitter(alpha = 0.8, position = position_jitter(width = 0.2)) +  
        scale_color_manual(values = colors, guide = "none") +  
        scale_x_discrete(labels = stage_labels) +  
        scale_y_continuous(
                labels = scales::label_number(accuracy = 0.1)  
        ) +
        labs(x = "Stage", y = "TAI") +  
        theme_minimal() +  
        theme(
                text = element_text(family = "Arial"),  
                legend.position = "none",  
                axis.text = element_text(size = 12),  
                axis.title = element_text(size = 16)   
        )
