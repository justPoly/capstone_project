# Load required libraries
library(dplyr)
library(tidyr)
library(corrplot)
library(caret)  # For findCorrelation()
library(stats)  # For Chi-square test
library(ggplot2)
library(tinytex)
# Load data
data <- read.csv("brain-cancer-dataset.csv")

# Handle missing values
data <- data %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  mutate(across(where(is.character), ~ ifelse(is.na(.), "Unknown", .))) %>%
  mutate(across(where(is.factor), ~ ifelse(is.na(.), as.character(levels(.)[1]), .)))

# Check for any remaining missing values
if (sum(is.na(data)) > 0) {
  stop("Unresolved missing values in the dataset.")
}

# Convert categorical variables to factors
data <- data %>% mutate(across(where(is.character), as.factor))

# Calculate variance for numeric columns and filter
variances <- apply(data[, sapply(data, is.numeric)], 2, var, na.rm = TRUE)
data_filtered <- data[, variances > 0.01]

# Compute and visualize correlation matrix
cor_matrix <- cor(data_filtered[, sapply(data_filtered, is.numeric)], use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.7)

# Remove highly correlated features (threshold: 0.8)
high_corr <- findCorrelation(cor_matrix, cutoff = 0.8)
data_filtered <- data_filtered[, -high_corr]

# Trim and clean column names
colnames(data) <- trimws(colnames(data))
colnames(data) <- make.names(colnames(data), unique = TRUE)

# Ensure factors for specific columns
data$OS.Status <- as.factor(data$OS.Status)
data$DFS.Status <- as.factor(data$DFS.Status)

# Define categorical columns
categorical_columns <- c(
  "Study.ID", "Patient.ID", "Sample.ID", "Age.Class", "BRAF_RELA.Status", 
  "BRAF.Status", "BRAF.Status2", "Cancer.Predispositions", "Cancer.Type", 
  "Cancer.Type.Detailed", "Chemotherapy", "Chemotherapy.Agents", 
  "Chemotherapy.Type", "Clinical.Status.at.Collection.Event", 
  "CTNNB1.Status", "DFS.Status", "Ethnicity", "Extent.of.Tumor.Resection", 
  "External.Patient.ID", "Formulation", "H3F3A_CTNNB1.Status", 
  "Initial.CNS.Tumor.Diagnosis.Related.to.OS", "Initial.Diagnosis.Type", 
  "LGG_BRAF.Status", "Medical.Conditions", "Multiple.Cancer.Predispositions", 
  "Multiple.Medical.Conditions", "Multiple.Tumor.Locations", "Oncotree.Code", 
  "OS.Status", "Protocol.and.Treatment.Arm", "Race", "Radiation", 
  "Radiation.Site", "Radiation.Type", "Sample.Annotation", "Sample.Origin", 
  "Sex", "Surgery", "Treatment", "Treatment.Changed", "Treatment.Status", 
  "Tumor.Location.Condensed", "Tumor.Tissue.Site", "Tumor.Type", "Updated.Grade"
)

# Ensure all categorical columns are factors
data <- data %>%
  mutate(across(all_of(categorical_columns), as.factor))

# Chi-square test for OS.Status and DFS.Status
chi_sq_results_OS <- lapply(categorical_columns, function(col) {
  table_data <- table(data[[col]], data$OS.Status)
  if (any(dim(table_data) == 0)) {
    return(data.frame(Feature = col, p_value = NA))
  }
  test_result <- chisq.test(table_data)
  data.frame(Feature = col, p_value = test_result$p.value)
})

chi_sq_summary_OS <- do.call(rbind, chi_sq_results_OS)

chi_sq_results_DFS <- lapply(categorical_columns, function(col) {
  table_data <- table(data[[col]], data$DFS.Status)
  if (any(dim(table_data) == 0)) {
    return(data.frame(Feature = col, p_value = NA))
  }
  test_result <- chisq.test(table_data)
  data.frame(Feature = col, p_value = test_result$p.value)
})

chi_sq_summary_DFS <- do.call(rbind, chi_sq_results_DFS)

# Filter significant features
significant_features_OS <- chi_sq_summary_OS %>% filter(p_value < 0.05)
significant_features_DFS <- chi_sq_summary_DFS %>% filter(p_value < 0.05)

# Subset data for significant features
data_OS <- data %>% select(all_of(significant_features_OS$Feature), OS.Status)
data_DFS <- data %>% select(all_of(significant_features_DFS$Feature), DFS.Status)

# Extract selected features
selected_features_OS <- colnames(data_OS)[colnames(data_OS) != "OS.Status"]
selected_features_DFS <- colnames(data_DFS)[colnames(data_DFS) != "DFS.Status"]

# Visualize distributions
ggplot(data, aes(x = OS.Status, fill = OS.Status)) +
  geom_bar() +
  labs(title = "Distribution of OS.Status", x = "OS.Status", y = "Count") +
  theme_minimal()

ggplot(data, aes(x = DFS.Status, fill = DFS.Status)) +
  geom_bar() +
  labs(title = "Distribution of DFS.Status", x = "DFS.Status", y = "Count") +
  theme_minimal()

# Normalize numeric columns
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

numeric_cols <- names(data)[sapply(data, is.numeric)]
data[numeric_cols] <- lapply(data[numeric_cols], normalize)

# Visualize correlation matrix
cor_matrix <- cor(data[numeric_cols], use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.8, number.cex = 0.7, addCoef.col = "black")

# Create a summary of feature counts
feature_counts <- data.frame(
  Target = c("OS.Status", "DFS.Status"),
  Features = c(length(selected_features_OS), length(selected_features_DFS)),
  Common = length(intersect(selected_features_OS, selected_features_DFS))
)

# Create a bar plot
ggplot(feature_counts, aes(x = Target, y = Features, fill = Target)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  geom_text(aes(label = Features), vjust = -0.5, size = 4) +
  labs(
    title = "Number of Features Selected for OS.Status and DFS.Status",
    x = "Target Variable",
    y = "Number of Features"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")


# Categorize features
feature_categories <- data.frame(
  Feature = union(selected_features_OS, selected_features_DFS),
  Category = ifelse(
    union(selected_features_OS, selected_features_DFS) %in% intersect(selected_features_OS, selected_features_DFS), "Both",
    ifelse(union(selected_features_OS, selected_features_DFS) %in% selected_features_OS, "OS.Status", "DFS.Status")
  )
)

# Display the table
knitr::kable(
  feature_categories,
  caption = "Feature Categories by Target Variable",
  col.names = c("Feature Name", "Category")
)


