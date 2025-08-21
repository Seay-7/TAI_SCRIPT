# Load necessary packages
library(readxl)       
library(survival)     
library(survminer)    
library(dplyr)        
library(maxstat)      


excel_data <- read_excel("cancer_patient_TAI.xlsx")

# Read clinical data text file
txt_data <- read.delim("Cancer\\XXX\\Merge_clinical.txt", 
                       sep = "\t", 
                       stringsAsFactors = FALSE)

# Merge datasets by sample ID
merged_data <- merge(excel_data, txt_data, 
                     by.x = "Sample", 
                     by.y = "A0_Samples")

# Prepare survival data (time, status, and TAI values)
surv_data <- merged_data %>% 
        select(A1_OS, A2_Event, TAI) %>% 
        mutate(
                status = as.numeric(A2_Event == "Dead")  # Convert survival status to 0/1 (0=censored, 1=dead)
        )

# Use maximum selected rank statistics (MaxStat) to find optimal threshold
maxstat_result <- maxstat.test(
        Surv(time = A1_OS, event = status) ~ TAI,
        data = surv_data,
        smethod = "LogRank"
)

# Extract optimal threshold
optimal_threshold <- maxstat_result$estimate
cat("Optimal TAI threshold:", optimal_threshold, "\n")  

# Group samples by optimal threshold
merged_data <- merged_data %>%
        mutate(
                TAI_Group = case_when(
                        TAI <= optimal_threshold ~ "Low-TAI",
                        TAI > optimal_threshold ~ "High-TAI",
                        TRUE ~ NA_character_
                )
        ) %>%
        filter(!is.na(TAI_Group))  # Remove samples with NA TAI values

# Set factor levels for consistent plotting
merged_data$TAI_Group <- factor(merged_data$TAI_Group, levels = c("Low-TAI", "High-TAI"))

# Check sample sizes in each group
print(table(merged_data$TAI_Group))


# Survival analysis
final_data <- merged_data %>%
        filter(A1_OS <= xxx)  

# Create survival object
surv_obj <- Surv(time = final_data$A1_OS, event = final_data$A2_Event == "Dead")

# Fit Kaplan-Meier curves
km_fit <- survfit(surv_obj ~ TAI_Group, data = final_data)

# Custom theme to center plot title
custom_theme <- theme_survminer(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5))

# Plot survival curves
ggsurvplot(
        km_fit,
        data = final_data,
        pval = TRUE,
        pval.method = FALSE,  
        pval.coord = c(0.1, 0.2),
        conf.int = TRUE,       
        risk.table = TRUE,     
        palette = c("#377EB8", "#CA1512"),  
        title = "Survival by TAI Group",
        xlab = "Time (Days)",
        ylab = "Survival Probability",
        break.time.by = 500, 
        legend.labs = c("Low-TAI", "High-TAI"),
        legend.title = "Group",
        ggtheme = custom_theme
)

# COX proportional hazards analysis
# Select variables for COX model - adjust based on your actual data
cox_data <- final_data %>%
        select(A1_OS, A2_Event, TAI_Group, A17_Age, A18_Sex, A6_Stage) %>%
        mutate(
                status = as.numeric(A2_Event == "Dead")  # Convert to 0/1 status
        )
cox_data$A17_Age <- as.numeric(cox_data$A17_Age)  # Ensure age is numeric

# Create survival object for COX model
surv_obj_cox <- Surv(time = cox_data$A1_OS, event = cox_data$status)

# Fit COX proportional hazards model
cox_model <- coxph(surv_obj_cox ~ TAI_Group + A17_Age + A18_Sex, data = cox_data)

# Output model results
summary(cox_model)

# Extract hazard ratios (HR) and confidence intervals (CI)
hr <- exp(coef(cox_model))
ci <- exp(confint(cox_model))

# Print results
cat("Hazard Ratios:\n")
print(hr)
cat("Confidence Intervals:\n")
print(ci)
    