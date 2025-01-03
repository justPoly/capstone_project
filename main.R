# Load required libraries
library(dplyr)
library(tidyr)
library(corrplot)
library(caret)  # For findCorrelation()
library(stats)  # For Chi-square test
library(ggplot2)
library(tinytex)
library(quarto)
library(knitr)
library(tidyverse)
library(pROC)
library(glmnet)
library(reshape2)# Melt data for easier plotting

# Load data
my_data <- read.csv("brain-cancer-dataset.csv")

# Handle missing values
my_data <- my_data %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), median(., na.rm = TRUE), .))) %>%
  mutate(across(where(is.character), ~ ifelse(is.na(.), "Unknown", .))) %>%
  mutate(across(where(is.factor), ~ ifelse(is.na(.), as.character(levels(.)[1]), .)))

# Check for any remaining missing values
if (sum(is.na(my_data)) > 0) {
  stop("Unresolved missing values in the dataset.")
}

# Convert categorical variables to factors
my_data <- my_data %>% mutate(across(where(is.character), as.factor))

# Calculate variance for numeric columns and filter
variances <- apply(my_data[, sapply(my_data, is.numeric)], 2, var, na.rm = TRUE)
data_filtered <- my_data[, variances > 0.01]

# Compute and visualize correlation matrix
cor_matrix <- cor(data_filtered[, sapply(data_filtered, is.numeric)], use = "complete.obs")
corrplot(cor_matrix, method = "color", tl.cex = 0.7)

# Remove highly correlated features (threshold: 0.8)
high_corr <- findCorrelation(cor_matrix, cutoff = 0.8)
data_filtered <- data_filtered[, -high_corr]

# Trim and clean column names
colnames(my_data) <- trimws(colnames(my_data))
colnames(my_data) <- make.names(colnames(my_data), unique = TRUE)

# Ensure factors for specific columns
my_data$OS.Status <- as.factor(my_data$OS.Status)
my_data$DFS.Status <- as.factor(my_data$DFS.Status)

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
my_data <- my_data %>%
  mutate(across(all_of(categorical_columns), as.factor))

# Chi-square test for OS.Status and DFS.Status
chi_sq_results_OS <- lapply(categorical_columns, function(col) {
  table_data <- table(my_data[[col]], my_data$OS.Status)
  if (any(dim(table_data) == 0)) {
    return(data.frame(Feature = col, p_value = NA))
  }
  test_result <- chisq.test(table_data)
  data.frame(Feature = col, p_value = test_result$p.value)
})

chi_sq_summary_OS <- do.call(rbind, chi_sq_results_OS)

chi_sq_results_DFS <- lapply(categorical_columns, function(col) {
  table_data <- table(my_data[[col]], my_data$DFS.Status)
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
data_OS <- my_data %>% select(all_of(significant_features_OS$Feature), OS.Status)
data_DFS <- my_data %>% select(all_of(significant_features_DFS$Feature), DFS.Status)

# Extract selected features
selected_features_OS <- colnames(data_OS)[colnames(data_OS) != "OS.Status"]
selected_features_DFS <- colnames(data_DFS)[colnames(data_DFS) != "DFS.Status"]

# Visualize distributions
ggplot(my_data, aes(x = OS.Status, fill = OS.Status)) +
  geom_bar() +
  labs(title = "Distribution of OS.Status", x = "OS.Status", y = "Count") +
  theme_minimal()

ggplot(my_data, aes(x = DFS.Status, fill = DFS.Status)) +
  geom_bar() +
  labs(title = "Distribution of DFS.Status", x = "DFS.Status", y = "Count") +
  theme_minimal()

# Normalize numeric columns
normalize <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

numeric_cols <- names(my_data)[sapply(my_data, is.numeric)]
my_data[numeric_cols] <- lapply(my_data[numeric_cols], normalize)

# Visualize correlation matrix
cor_matrix <- cor(my_data[numeric_cols], use = "complete.obs")
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

# Set seed for reproducibility
set.seed(123)

# Split the data for OS.Status (80% for training, 20% for testing)
train_OS <- createDataPartition(data_OS$OS.Status, p = 0.8, list = FALSE)
train_OS_data <- data_OS[train_OS, ]
test_OS_data <- data_OS[-train_OS, ]

# Identify columns with only one level (factor columns)
single_level_factors <- sapply(train_OS_data, function(x) is.factor(x) && length(unique(x)) == 1)
single_level_columns <- names(single_level_factors[single_level_factors])
print(single_level_columns)

# Remove columns with only one level
train_OS_data <- train_OS_data[, !colnames(train_OS_data) %in% single_level_columns]

# Set up the training control
train_control <- trainControl(method = "cv", number = 10)

# Train the logistic regression model
model_logistic_OS <- train(OS.Status ~ ., data = train_OS_data, method = "glm", family = "binomial", trControl = train_control)

# Check model summary
summary(model_logistic_OS)

# Predict on the training data
predictions_train <- predict(model_logistic_OS, newdata = train_OS_data)

# Confusion matrix to evaluate performance
conf_matrix_train <- confusionMatrix(predictions_train, train_OS_data$OS.Status)
print(conf_matrix_train)

# Predict on the test data
predictions_test <- predict(model_logistic_OS, newdata = test_OS_data)

# Confusion matrix to evaluate performance on the test data
conf_matrix_test <- confusionMatrix(predictions_test, test_OS_data$OS.Status)
print(conf_matrix_test)

# Feature importance
feature_importance <- varImp(model_logistic_OS, scale = FALSE)

# Print feature importance
print(feature_importance)

# Plot feature importance
plot(feature_importance)

# Cross-validation results
cv_results <- model_logistic_OS$resample
print(cv_results)

# Plot cross-validation results
ggplot(cv_results, aes(x = Resample, y = Accuracy)) +
  geom_boxplot() +
  labs(title = "Cross-validation Results - Accuracy", y = "Accuracy", x = "Fold") +
  theme_minimal()

# Get predicted probabilities (needed for ROC curve)
probabilities <- predict(model_logistic_OS, newdata = train_OS_data, type = "prob")[, 2]

# ROC curve and AUC
roc_curve <- roc(train_OS_data$OS.Status, probabilities)
plot(roc_curve, main = "ROC Curve for OS.Status", col = "blue")
auc(roc_curve)  # This will give you the AUC value

# Example of hyperparameter tuning for logistic regression using `caret`
tune_grid <- expand.grid(alpha = 0:1, lambda = seq(0, 1, by = 0.1))

model_tuned <- train(OS.Status ~ ., data = train_OS_data, method = "glmnet", 
                     trControl = train_control, tuneGrid = tune_grid, family = "binomial")

# Save the final model for later use
saveRDS(model_logistic_OS, "model_logistic_OS.rds")

# Load the model (when needed)
# model_logistic_OS <- readRDS("model_logistic_OS.rds")

# Split the data for DFS.Status (80% for training, 20% for testing)
train_DFS <- createDataPartition(data_DFS$DFS.Status, p = 0.8, list = FALSE)
train_DFS_data <- data_DFS[train_DFS, ]
test_DFS_data <- data_DFS[-train_DFS, ]

# Identify columns with only one level (factor columns)
single_level_factors_DFS <- sapply(train_DFS_data, function(x) is.factor(x) && length(unique(x)) == 1)
single_level_columns_DFS <- names(single_level_factors_DFS[single_level_factors_DFS])
print(single_level_columns_DFS)

# Remove columns with only one level
train_DFS_data <- train_DFS_data[, !colnames(train_DFS_data) %in% single_level_columns_DFS]

# Set up the training control for DFS.Status
train_control_DFS <- trainControl(method = "cv", number = 10)

# Train the logistic regression model for DFS.Status
model_logistic_DFS <- train(DFS.Status ~ ., data = train_DFS_data, method = "glm", family = "binomial", trControl = train_control_DFS)

# Check model summary for DFS.Status
summary(model_logistic_DFS)

# Predict on the training data for DFS.Status
predictions_train_DFS <- predict(model_logistic_DFS, newdata = train_DFS_data)

# Confusion matrix to evaluate performance on training data for DFS.Status
conf_matrix_train_DFS <- confusionMatrix(predictions_train_DFS, train_DFS_data$DFS.Status)
print(conf_matrix_train_DFS)

# Predict on the test data for DFS.Status
predictions_test_DFS <- predict(model_logistic_DFS, newdata = test_DFS_data)

# Confusion matrix to evaluate performance on test data for DFS.Status
conf_matrix_test_DFS <- confusionMatrix(predictions_test_DFS, test_DFS_data$DFS.Status)
print(conf_matrix_test_DFS)

# Feature importance for DFS.Status
feature_importance_DFS <- varImp(model_logistic_DFS, scale = FALSE)

# Print feature importance for DFS.Status
print(feature_importance_DFS)

# Plot feature importance for DFS.Status
plot(feature_importance_DFS)

# Cross-validation results for DFS.Status
cv_results_DFS <- model_logistic_DFS$resample
print(cv_results_DFS)

# Plot cross-validation results for DFS.Status
ggplot(cv_results_DFS, aes(x = Resample, y = Accuracy)) +
  geom_boxplot() +
  labs(title = "Cross-validation Results - Accuracy (DFS.Status)", y = "Accuracy", x = "Fold") +
  theme_minimal()

# Get predicted probabilities for DFS.Status (needed for ROC curve)
probabilities_DFS <- predict(model_logistic_DFS, newdata = train_DFS_data, type = "prob")[, 2]

# ROC curve and AUC for DFS.Status
roc_curve_DFS <- roc(train_DFS_data$DFS.Status, probabilities_DFS)
plot(roc_curve_DFS, main = "ROC Curve for DFS.Status", col = "blue")
auc(roc_curve_DFS)  # This will give you the AUC value

# Hyperparameter tuning for DFS.Status using glmnet
tune_grid_DFS <- expand.grid(alpha = 0:1, lambda = seq(0, 1, by = 0.1))

model_tuned_DFS <- train(DFS.Status ~ ., data = train_DFS_data, method = "glmnet", 
                         trControl = train_control_DFS, tuneGrid = tune_grid_DFS, family = "binomial")

#--Save Model here
# Save the final DFS.Status model
saveRDS(model_logistic_DFS, "model_logistic_DFS.rds")

# Load the model (when needed)
# model_logistic_DFS <- readRDS("model_logistic_DFS.rds")

# Train Random Forest model for OS.Status
model_rf_OS <- train(OS.Status ~ ., data = train_OS_data, method = "rf", trControl = train_control)

# Check model summary for Random Forest
summary(model_rf_OS)

# Get model evaluation results for Random Forest
print(model_rf_OS)

# Predict on the training data (or test data)
predictions_rf_OS <- predict(model_rf_OS, newdata = train_OS_data)

# Confusion matrix to evaluate Random Forest performance
conf_matrix_rf_OS <- confusionMatrix(predictions_rf_OS, train_OS_data$OS.Status)
print(conf_matrix_rf_OS)

# Plot feature importance for Random Forest
feature_importance_rf <- varImp(model_rf_OS, scale = FALSE)
print(feature_importance_rf)
plot(feature_importance_rf)

# Cross-validation results for Random Forest
cv_results_rf <- model_rf_OS$resample
print(cv_results_rf)

# Plot cross-validation results for Random Forest
ggplot(cv_results_rf, aes(x = Resample, y = Accuracy)) +
  geom_boxplot() +
  labs(title = "Cross-validation Results - Accuracy (Random Forest)", y = "Accuracy", x = "Fold") +
  theme_minimal()

# ROC curve and AUC for Random Forest
probabilities_rf <- predict(model_rf_OS, newdata = train_OS_data, type = "prob")[, 2]
roc_curve_rf <- roc(train_OS_data$OS.Status, probabilities_rf)
plot(roc_curve_rf, main = "ROC Curve for OS.Status (Random Forest)", col = "blue")
auc(roc_curve_rf)  # This will give you the AUC value

# Train Random Forest model for DFS.Status
model_rf_DFS <- train(DFS.Status ~ ., data = train_DFS_data, method = "rf", trControl = train_control)

# Check model summary for Random Forest
summary(model_rf_DFS)

# Get model evaluation results for Random Forest
print(model_rf_DFS)

# Predict on the training data (or test data)
predictions_rf_DFS <- predict(model_rf_DFS, newdata = train_DFS_data)

# Confusion matrix to evaluate Random Forest performance
conf_matrix_rf_DFS <- confusionMatrix(predictions_rf_DFS, train_DFS_data$DFS.Status)
print(conf_matrix_rf_DFS)

# Plot feature importance for Random Forest
feature_importance_rf_DFS <- varImp(model_rf_DFS, scale = FALSE)
print(feature_importance_rf_DFS)
plot(feature_importance_rf_DFS)

# Cross-validation results for Random Forest
cv_results_rf_DFS <- model_rf_DFS$resample
print(cv_results_rf_DFS)

# Plot cross-validation results for Random Forest
ggplot(cv_results_rf_DFS, aes(x = Resample, y = Accuracy)) +
  geom_boxplot() +
  labs(title = "Cross-validation Results - Accuracy (Random Forest for DFS)", y = "Accuracy", x = "Fold") +
  theme_minimal()

# ROC curve and AUC for Random Forest
probabilities_rf_DFS <- predict(model_rf_DFS, newdata = train_DFS_data, type = "prob")[, 2]
roc_curve_rf_DFS <- roc(train_DFS_data$DFS.Status, probabilities_rf_DFS)
plot(roc_curve_rf_DFS, main = "ROC Curve for DFS.Status (Random Forest)", col = "blue")
auc(roc_curve_rf_DFS)  # This will give you the AUC value

# Compare Logistic Regression and Random Forest for OS.Status
model_comparison_OS <- data.frame(
  Model = c("Logistic Regression", "Random Forest"),
  Accuracy = c(conf_matrix$overall['Accuracy'], conf_matrix_rf_OS$overall['Accuracy']),
  AUC = c(auc(roc_curve), auc(roc_curve_rf))
)

print(model_comparison_OS)

# Compare Logistic Regression and Random Forest for DFS.Status
model_comparison_DFS <- data.frame(
  Model = c("Logistic Regression", "Random Forest"),
  Accuracy = c(conf_matrix_train_DFS$overall['Accuracy'], conf_matrix_rf_DFS$overall['Accuracy']),
  AUC = c(auc(roc_curve_DFS), auc(roc_curve_rf_DFS))
)

print(model_comparison_DFS)

# Predictions on the test data for OS.Status (Logistic Regression)
predictions_test_OS_logistic <- predict(model_logistic_OS, newdata = test_OS_data)
conf_matrix_test_OS_logistic <- confusionMatrix(predictions_test_OS_logistic, test_OS_data$OS.Status)
print(conf_matrix_test_OS_logistic)

# Predictions on the test data for OS.Status (Random Forest)
predictions_test_OS_rf <- predict(model_rf_OS, newdata = test_OS_data)
conf_matrix_test_OS_rf <- confusionMatrix(predictions_test_OS_rf, test_OS_data$OS.Status)
print(conf_matrix_test_OS_rf)

# Predictions on the test data for DFS.Status (Logistic Regression)
predictions_test_DFS_logistic <- predict(model_logistic_DFS, newdata = test_DFS_data)
conf_matrix_test_DFS_logistic <- confusionMatrix(predictions_test_DFS_logistic, test_DFS_data$DFS.Status)
print(conf_matrix_test_DFS_logistic)

# Predictions on the test data for DFS.Status (Random Forest)
predictions_test_DFS_rf <- predict(model_rf_DFS, newdata = test_DFS_data)
conf_matrix_test_DFS_rf <- confusionMatrix(predictions_test_DFS_rf, test_DFS_data$DFS.Status)
print(conf_matrix_test_DFS_rf)

# Compare performance for OS.Status models (Logistic Regression vs Random Forest)
model_comparison_OS_test <- data.frame(
  Model = c("Logistic Regression", "Random Forest"),
  Accuracy = c(conf_matrix_test_OS_logistic$overall['Accuracy'], conf_matrix_test_OS_rf$overall['Accuracy']),
  AUC = c(roc(test_OS_data$OS.Status, predict(model_logistic_OS, newdata = test_OS_data, type = "prob")[, 2])$auc,
          roc(test_OS_data$OS.Status, predict(model_rf_OS, newdata = test_OS_data, type = "prob")[, 2])$auc)
)
print(model_comparison_OS_test)

# Compare performance for DFS.Status models (Logistic Regression vs Random Forest)
model_comparison_DFS_test <- data.frame(
  Model = c("Logistic Regression", "Random Forest"),
  Accuracy = c(conf_matrix_test_DFS_logistic$overall['Accuracy'], conf_matrix_test_DFS_rf$overall['Accuracy']),
  AUC = c(roc(test_DFS_data$DFS.Status, predict(model_logistic_DFS, newdata = test_DFS_data, type = "prob")[, 2])$auc,
          roc(test_DFS_data$DFS.Status, predict(model_rf_DFS, newdata = test_DFS_data, type = "prob")[, 2])$auc)
)
print(model_comparison_DFS_test)

# OS.Status performance
ggplot(model_comparison_OS_test, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Model Comparison - OS.Status", x = "Model", y = "Accuracy")

# DFS.Status performance
ggplot(model_comparison_DFS_test, aes(x = Model, y = Accuracy, fill = Model)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Model Comparison - DFS.Status", x = "Model", y = "Accuracy")

# Model Tuning to Improve Performance
# Hyperparameter tuning for Random Forest (DFS.Status)
tune_rf_DFS <- train(
  DFS.Status ~ ., 
  data = train_DFS_data,
  method = "rf",
  trControl = train_control,
  tuneGrid = expand.grid(mtry = seq(2, sqrt(ncol(train_DFS_data)), by = 1))
)

# Hyperparameter tuning for Random Forest (OS.Status)
tune_rf_OS <- train(
  OS.Status ~ ., 
  data = train_OS_data,
  method = "rf",
  trControl = train_control,
  tuneGrid = expand.grid(mtry = seq(2, sqrt(ncol(train_OS_data)), by = 1))
)

# Check the best parameters
print(tune_rf_DFS$bestTune)
print(tune_rf_OS$bestTune)

# Predictions for DFS.Status
predictions_rf_DFS_test <- predict(tune_rf_DFS, newdata = test_DFS_data)

# Predictions for OS.Status
predictions_rf_OS_test <- predict(tune_rf_OS, newdata = test_OS_data)

# DFS.Status Evaluation
conf_matrix_rf_DFS_test <- confusionMatrix(predictions_rf_DFS_test, test_DFS_data$DFS.Status)
print(conf_matrix_rf_DFS_test)

# OS.Status Evaluation
conf_matrix_rf_OS_test <- confusionMatrix(predictions_rf_OS_test, test_OS_data$OS.Status)
print(conf_matrix_rf_OS_test)

# DFS.Status AUC
probabilities_rf_DFS_test <- predict(tune_rf_DFS, newdata = test_DFS_data, type = "prob")[, 2]
roc_rf_DFS_test <- roc(test_DFS_data$DFS.Status, probabilities_rf_DFS_test)
auc_rf_DFS_test <- auc(roc_rf_DFS_test)
print(auc_rf_DFS_test)

# OS.Status AUC
probabilities_rf_OS_test <- predict(tune_rf_OS, newdata = test_OS_data, type = "prob")[, 2]
roc_rf_OS_test <- roc(test_OS_data$OS.Status, probabilities_rf_OS_test)
auc_rf_OS_test <- auc(roc_rf_OS_test)
print(auc_rf_OS_test)

# Compare performance for OS.Status models (Logistic Regression vs Tuned Random Forest)
model_comparison_OS_test <- data.frame(
  Model = c("Logistic Regression", "Tuned Random Forest"),
  Accuracy = c(
    conf_matrix_test_OS_logistic$overall['Accuracy'], 
    confusionMatrix(predict(tune_rf_OS, newdata = test_OS_data), test_OS_data$OS.Status)$overall['Accuracy']
  ),
  AUC = c(
    roc(test_OS_data$OS.Status, predict(model_logistic_OS, newdata = test_OS_data, type = "prob")[, 2])$auc,
    roc(test_OS_data$OS.Status, predict(tune_rf_OS, newdata = test_OS_data, type = "prob")[, 2])$auc
  )
)
print(model_comparison_OS_test)

# Compare performance for DFS.Status models (Logistic Regression vs Tuned Random Forest)
model_comparison_DFS_test <- data.frame(
  Model = c("Logistic Regression", "Tuned Random Forest"),
  Accuracy = c(
    conf_matrix_test_DFS_logistic$overall['Accuracy'], 
    confusionMatrix(predict(tune_rf_DFS, newdata = test_DFS_data), test_DFS_data$DFS.Status)$overall['Accuracy']
  ),
  AUC = c(
    roc(test_DFS_data$DFS.Status, predict(model_logistic_DFS, newdata = test_DFS_data, type = "prob")[, 2])$auc,
    roc(test_DFS_data$DFS.Status, predict(tune_rf_DFS, newdata = test_DFS_data, type = "prob")[, 2])$auc
  )
)
print(model_comparison_DFS_test)

# Prepare data for visualization (OS.Status)
visualization_data_OS <- data.frame(
  Model = c("Logistic Regression", "Tuned Random Forest"),
  Accuracy = c(
    conf_matrix_test_OS_logistic$overall['Accuracy'], 
    confusionMatrix(predict(tune_rf_OS, newdata = test_OS_data), test_OS_data$OS.Status)$overall['Accuracy']
  ),
  AUC = c(
    roc(test_OS_data$OS.Status, predict(model_logistic_OS, newdata = test_OS_data, type = "prob")[, 2])$auc,
    roc(test_OS_data$OS.Status, predict(tune_rf_OS, newdata = test_OS_data, type = "prob")[, 2])$auc
  )
)

visualization_data_OS_melted <- melt(visualization_data_OS, id.vars = "Model", variable.name = "Metric", value.name = "Value")

# Plot Accuracy and AUC for OS.Status
ggplot(visualization_data_OS_melted, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Performance Comparison for OS.Status Models",
       x = "Model",
       y = "Performance Metric") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")


# Prepare data for visualization (DFS.Status)
visualization_data_DFS <- data.frame(
  Model = c("Logistic Regression", "Tuned Random Forest"),
  Accuracy = c(
    conf_matrix_test_DFS_logistic$overall['Accuracy'], 
    confusionMatrix(predict(tune_rf_DFS, newdata = test_DFS_data), test_DFS_data$DFS.Status)$overall['Accuracy']
  ),
  AUC = c(
    roc(test_DFS_data$DFS.Status, predict(model_logistic_DFS, newdata = test_DFS_data, type = "prob")[, 2])$auc,
    roc(test_DFS_data$DFS.Status, predict(tune_rf_DFS, newdata = test_DFS_data, type = "prob")[, 2])$auc
  )
)

# Melt data for easier plotting
visualization_data_DFS_melted <- melt(visualization_data_DFS, id.vars = "Model", variable.name = "Metric", value.name = "Value")

# Plot Accuracy and AUC for DFS.Status
ggplot(visualization_data_DFS_melted, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Performance Comparison for DFS.Status Models",
       x = "Model",
       y = "Performance Metric") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set1")

