# Load required packages
library(readxl)
library(myTAI)
library(openxlsx)

# Read Excel file
data <- read_excel("Cancer\\XXX\\Merge_RNA_seq_StageX_FPKM.xlsx")

# Get sample names
sample_names <- colnames(data)[-c(1, 2)]

# Initialize an empty data frame to store results
result <- data.frame(Sample = character(), TAI = numeric())

# Loop to calculate TAI value for each sample
for (sample in sample_names) {
        # Extract gene age, gene, and expression level of current sample
        subset_data <- data[, c("AGE", "Tag", sample)]
        colnames(subset_data) <- c("AGE", "Tag", "Expression")
        
        # Calculate TAI value
        tai_value <- TAI(subset_data)
        
        # Add result to the result data frame
        new_row <- data.frame(Sample = sample, TAI = tai_value)
        result <- rbind(result, new_row)
}

# Save results to Excel file
write.xlsx(result, "Cancer\\XXX\\Merge_RNA_seq_StageX_sampleTAI.xlsx")