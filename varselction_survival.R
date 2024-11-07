##### Variable Importance Selection #####
library(Boruta)
# Set random seed for reproducibility
set.seed(123)
boruta_output <- Boruta(status ~ ., data=msigma, pValue = 0.01, mcAdj = TRUE, doTrace = 2)
summary(boruta_output)

# Extract variable importance values from Boruta results
importance <- boruta_output$ImpHistory 
importance

# Save Boruta plot as a PDF
pdf(paste0("./boruta", name, ".pdf"), height = 10, width = 10)
plot(boruta_output, las = 2, xlab = '', main = 'Variable Importance')
dev.off()

# Save the importance data as a CSV file
write.csv(importance, file = paste0("./borutavip", name, ".csv"), row.names = FALSE)

# Refine tentative features using the TentativeRoughFix function
bor <- TentativeRoughFix(boruta_output)
print(bor)

# Display the final results of feature selection
attStats(bor)
write.csv(attStats(bor), file = paste0('./borutavip', name, "-attStats.csv"))

# Extract all important features
getNonRejectedFormula(bor)

# Function to prepare Boruta importance data for further analysis
boruta.imp <- function(x) {
  imp <- reshape2::melt(x$ImpHistory, na.rm = TRUE)[, -1]
  colnames(imp) <- c("Variable", "Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable = names(x$finalDecision), 
                            finalDecision = x$finalDecision)
  
  # Add shadow feature groups for reference
  showGrp <- data.frame(Variable = c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision = c("shadowMax", "shadowMean", "shadowMin"))
  variableGrp <- rbind(variableGrp, showGrp)
  
  # Merge data to show importance scores with decisions
  boruta.variable.imp <- merge(imp, variableGrp, all.x = TRUE)
  
  # Sort variables by median importance
  sortedVariable <- boruta.variable.imp %>%
    group_by(Variable) %>%
    summarise(median = median(Importance)) %>%
    arrange(median) %>%
    .$Variable %>%
    as.vector()
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels = sortedVariable)
  invisible(boruta.variable.imp)
}

# Generate Boruta importance data and display a sample
boruta.variable.imp <- boruta.imp(boruta_output)
head(boruta.variable.imp)

# Prepare important features for survival analysis
expr_imp <- msigma[, rownames(attStats(bor)[attStats(bor)$decision == "Confirmed",])]
expr_imp[] <- lapply(expr_imp, as.numeric)
expr_imp$status <- as.numeric(msigma$status)
expr_imp$time <- unclass(round(msigma$time))

#### Survival Analysis ####
library(mlr3verse)
library(mlr3extralearners)
library(mlr3proba)
library(mlr3measures)
library(pseudo)
library(survivalmodels)

# Adjust time variable and set up survival analysis task
expr_imp$time <- expr_imp$time - floor(min(expr_imp$time) / 10) * 10
set.seed(111)
task <- as_task_surv(expr_imp, time = "time", event = "status")
split <- partition(task) # Split data into training and testing

# Train multiple survival models
lrn_cox <- lrn("surv.coxph")
lrn_cox$train(task, row_ids = split$train)

lrn_ranger <- lrn("surv.ranger")
lrn_ranger$train(task, row_ids = split$train)

lrn_deep <- lrn("surv.deepsurv")
lrn_deep$train(task, row_ids = split$train)

lrn_deephit <- lrn("surv.deephit")
lrn_deephit$train(task, row_ids = split$train)

lrn_xgboost <- lrn("surv.xgboost.cox")
lrn_xgboost$train(task, row_ids = split$train)

## Model Interpretation ##
library(survival)
library(survex)

# Prepare data for explanations using test data
credit_x <- task$data(rows = split$test, cols = task$feature_names)
credit_y <- task$data(rows = split$test, cols = task$target_names)

# Create model explainers
ranger_explainer <- explain(lrn_ranger, data = credit_x, y = Surv(credit_y$time, credit_y$status), label = "ranger model")
deepsurve_explainer <- explain(lrn_deep, data = credit_x, y = Surv(credit_y$time, credit_y$status), label = "deepsurve model")
cox_explainer <- explain(lrn_cox, data = credit_x, y = Surv(credit_y$time, credit_y$status), label = "cox model")
xgboost_explainer <- explain(lrn_xgboost, data = credit_x, y = Surv(credit_y$time, credit_y$status), label = "xgboost model")
deephit_explainer <- explain(lrn_deephit, data = credit_x, y = Surv(credit_y$time, credit_y$status), label = "deephit model")

## Model Performance Evaluation ##
perf_credit_ranger <- model_performance(ranger_explainer)
perf_credit_deepsurve <- model_performance(deepsurve_explainer)
perf_credit_deephit <- model_performance(deephit_explainer)
perf_credit_cox <- model_performance(cox_explainer)
perf_credit_xgboost <- model_performance(xgboost_explainer)

# Compile AUC results and sort
auc_RES <- list(
  "xgboost_explainer" = as.numeric(perf_credit_xgboost$result$`Integrated C/D AUC`),
  "ranger_explainer" = as.numeric(perf_credit_ranger$result$`Integrated C/D AUC`),
  "cox_explainer" = as.numeric(perf_credit_cox$result$`Integrated C/D AUC`),
  "deepsurve_explainer"= as.numeric(perf_credit_deepsurve$result$`Integrated C/D AUC`),
  "deephit_explainer" = as.numeric(perf_credit_deephit$result$`Integrated C/D AUC`)
)

# Sort AUC values and display the best-performing models
my_vector <- unlist(auc_RES)
sorted_vector <- sort(my_vector, decreasing = TRUE)
sorted_list <- as.list(sorted_vector)
names(sorted_list)

# Plot model performance over time and as a scalar metric
setwd("~/AD/survival")
pdf(paste0("modelperformance", name, ".pdf"), height = 7, width = 15)
plot(perf_credit_ranger, perf_credit_deepsurve, perf_credit_deephit, perf_credit_cox, perf_credit_xgboost)
dev.off()

pdf(paste0("modelperformance", name, "-scalar.pdf"), height = 7, width = 15)
plot(perf_credit_ranger, perf_credit_deepsurve, perf_credit_deephit, perf_credit_cox, perf_credit_xgboost, metrics_type = "scalar")
dev.off()

# Feature Importance Analysis
effect1 <- model_parts(ranger_explainer)
p1 <- plot(effect1)

effect2 <- model_parts(deepsurve_explainer)
p2 <- plot(effect2)

effect3 <- model_parts(cox_explainer)
p3 <- plot(effect3)

# Combine importance plots and save
library(cowplot)
pdf(paste0("model_parts", name, ".pdf"), height = 7, width = 20)
p1 + p2 + p3
dev.off()

# Rank features based on weighted importance
column_sums1 <- colSums(abs(effect1$result[, 3:length(effect1$result) - 3]))
column_sums2 <- colSums(abs(effect2$result[, 3:length(effect2$result) - 3]))
column_sums3 <- colSums(abs(effect3$result[, 3:length(effect3$result) - 3]))

rank_sum <- list()
rank_sum[names(sorted_list)[1]] <- list(column_sums1)
rank_sum[names(sorted_list)[2]] <- list(column_sums2)
rank_sum[names(sorted_list)[3]] <- list(column_sums3)
rank_sum <- data.frame(rank_sum)

# Compute weighted sum and sort for top features
rank_sum$sum <- column_sums1 * 1.0 + column_sums2 * 0.8 + column_sums3 * 0.6
rank_sum <- rank_sum[order(rank_sum$sum, decreasing = TRUE), ]

# Save ranked features
write.csv(rank_sum, file = paste0("surtop3-ranksum", name, ".csv"))

# Select top 20% features (or top 5 if fewer)
rank_sum <- rank_sum[!rownames(rank_sum) %in% c("_times_", "_full_model_"),]
top_20 <- rownames(rank_sum)[1:round(length(rank_sum$sum) * 0.2)]

# Display the top features
print(top_20)
