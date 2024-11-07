# Load necessary libraries for model setup
library(tidymodels)
library(DALEXtra)
library(iml)
library(forcats)

# Define various machine learning models with hyperparameter tuning
{# XGBoost
  xgb_mod <- boost_tree(tree_depth = tune(),
                        learn_rate = tune(), 
                        loss_reduction = tune(), 
                        min_n = tune(), 
                        sample_size = tune(), 
                        trees = 200L
  ) %>%           
    set_engine("xgboost") %>%           
    set_mode("classification")   
  
  # Decision Tree
  dt_mod <- decision_tree(cost_complexity = tune(), min_n = tune()) %>%           
    set_engine("rpart") %>%           
    set_mode("classification")
  
  # Logistic Regression
  logistic_mod <- logistic_reg(penalty = tune(), mixture = tune()) %>%
    set_engine('glmnet')
  
  # Neural Network (nnet)
  nnet_mod <- mlp(hidden_units = tune(), 
                  penalty = tune(), 
                  epochs = tune()) %>%          
    set_engine('nnet') %>%          
    set_mode('classification')
  
  # Naive Bayes
  naive_mod <- naive_Bayes(smoothness = tune(), Laplace = tune()) %>% 
    set_engine("naivebayes") %>% 
    translate()
  
  # Generalized Additive Model
  gen_mod <- gen_additive_mod(adjust_deg_free = tune(), select_features = tune()) %>% 
    set_engine("mgcv") %>% 
    set_mode("classification") %>% 
    translate()
  
  # Linear Discriminant Analysis
  dis_mod <- discrim_linear(penalty = tune()) %>% 
    set_engine("mda") %>% 
    translate()
  
  # K-Nearest Neighbors
  kknn_mod <- nearest_neighbor(neighbors = tune(), 
                               dist_power = tune(), 
                               weight_func = tune()) %>%          
    set_engine('kknn') %>%          
    set_mode('classification')
  
  # Random Forest
  rf_mod <- rand_forest(mtry = tune(), min_n = tune(), trees = tune()) %>%          
    set_engine('ranger', importance = "impurity") %>%          
    set_mode('classification')
  
  # Support Vector Machine (SVM)
  svm_mod <- svm_rbf(cost = tune(), rbf_sigma = tune()) %>%
    set_engine('kernlab') %>%
    set_mode('classification')
  
  # Partial Least Squares (PLS)
  pls_mod <- pls(num_comp = tune()) %>%
    set_engine('mixOmics') %>%
    set_mode('classification')
  
  # Lasso Regression
  lasso_mod <- logistic_reg(penalty = tune(), mixture = 1) %>% 
    set_engine("glmnet") %>% 
    set_mode("classification")
  
  # LightGBM
  lgb_mod <- boost_tree(trees = tune(),
                        min_n = tune(),
                        tree_depth = tune()) %>% 
    set_engine("lightgbm") %>%
    set_mode("classification")
  
  # Bagged Trees
  bag_mod <- bag_tree(tree_depth = tune(), min_n = tune(), cost_complexity = tune()) %>% 
    set_engine("rpart") %>% 
    set_mode("classification") %>% 
    translate()
  
  # Flexible Discriminant Analysis
  mars_disc_spec <- discrim_flexible(prod_degree = tune()) %>% 
    set_engine("earth")
}

# Data Preparation and Splitting
set.seed(123)
split_pbp <- initial_split(expr, 0.7, strata = Group)
train_data <- training(split_pbp) # Training set
test_data <- testing(split_pbp) # Test set

# Recipe for data preprocessing
pbp_rec <- recipe(formula = formula_string, data = train_data)  %>%
  step_corr(all_numeric(), threshold = 0.7) %>% # Remove highly correlated variables
  step_center(all_numeric()) %>%                # Center numeric values
  step_zv(all_predictors())                     # Remove zero-variance predictors

summary(pbp_rec)
print(tail(pbp_rec[["var_info"]][,]))

# Cross-validation setup
set.seed(20220520)
df_folds <- vfold_cv(train_data, strata = Group, repeats = 10)
keep_pred <- control_resamples(save_pred = TRUE, verbose = TRUE)

# Workflow setup for models
wf <- workflow_set(preproc = list(pbp_rec),          
                   models = list(
                     pls = pls_mod,
                     xgb = xgb_mod,          
                     dt = dt_mod,          
                     bag = bag_mod,     
                     svm = svm_mod,          
                     nnet = nnet_mod,          
                     knn = kknn_mod,          
                     rf = rf_mod,
                     lasso = lasso_mod,
                     da = mars_disc_spec))

# Set up parallel processing
library(parallel)
library(finetune)

cl <- makeCluster(4)
doParallel::registerDoParallel(cl)
set.seed(123)

# Simulated Annealing Tuning
grid_ctrl <- control_sim_anneal(save_pred = TRUE, save_workflow = TRUE, parallel_over = "everything")

system.time({
  grid_res <- workflow_map(object = wf,
                           fn = "tune_sim_anneal",
                           resamples = df_folds,
                           iter = 30, 
                           control = grid_ctrl,
                           verbose = TRUE)
})
stopCluster(cl)

#### Model Importance Interpretation ####
# Function to visualize variable importance
ggplot_imp <- function(...) {
  obj <- list(...)
  metric_name <- attr(obj[[1]], "loss_name")
  metric_lab <- paste(metric_name, "after permutations\n(higher indicates more important)")
  
  full_vip <- bind_rows(obj) %>%
    filter(variable != "_baseline_")
  
  perm_vals <- full_vip %>% 
    filter(variable == "_full_model_") %>% 
    group_by(label) %>% 
    aggregate(dropout_loss ~ label, FUN = mean)
  
  p <- full_vip %>%
    filter(variable != "_full_model_") %>% 
    mutate(variable = fct_reorder(variable, dropout_loss)) %>%
    ggplot(aes(dropout_loss, variable)) 
  
  if(length(obj) > 1) {
    p <- p + 
      facet_wrap(vars(label)) +
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss, color = label),
                 linewidth = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(aes(color = label, fill = label), alpha = 0.2)
  } else {
    p <- p + 
      geom_vline(data = perm_vals, aes(xintercept = dropout_loss),
                 linewidth = 1.4, lty = 2, alpha = 0.7) +
      geom_boxplot(fill = "#91CBD765", alpha = 0.4)
  }
  p +
    theme(legend.position = "none") +
    labs(x = metric_lab, y = NULL, fill = NULL, color = NULL)
}

# Function for model ranking and comparison
expmod <- function(rankdata, name1) {
  set.seed(123)
  split_pbp <- initial_split(expr, 0.7, strata = Group)
  train_data <- training(split_pbp) 
  test_data <- testing(split_pbp)
  
  pbp_rec <- recipe(formula = formula_string, data = train_data)  %>%
    step_corr(all_numeric(), threshold = 0.7) %>% 
    step_zv(all_predictors())  
  
  # Select and finalize top three models
  wf_1 <- grid_res %>% extract_workflow_set_result(rankdata$Row.names[1]) %>% select_best(metric = "roc_auc")  
  final_wf1 <- grid_res %>% extract_workflow(rankdata$Row.names[1]) %>% finalize_workflow(wf_1) %>% update_recipe(pbp_rec) %>% fit(data = test_data)
  
  wf_2 <- grid_res %>% extract_workflow_set_result(rankdata$Row.names[2]) %>% select_best(metric = "roc_auc")
  final_wf2 <- grid_res %>% extract_workflow(rankdata$Row.names[2]) %>% finalize_workflow(wf_2) %>% update_recipe(pbp_rec) %>% fit(data = test_data)
  
  wf_3 <- grid_res %>% extract_workflow_set_result(rankdata$Row.names[3]) %>% select_best(metric = "roc_auc")
  final_wf3 <- grid_res %>% extract_workflow(rankdata$Row.names[3]) %>% finalize_workflow(wf_3) %>% update_recipe(pbp_rec) %>% fit(data = test_data)
  
  # Generate explainers for each model
  explainer_1 <- explain_tidymodels(final_wf1, data = test_data[, colnames(test_data) != "Group"], label = mod[rankdata$Row.names[1]], y = as.numeric(test_data$Group))
  explainer_2 <- explain_tidymodels(final_wf2, data = test_data[, colnames(test_data) != "Group"], label = mod[rankdata$Row.names[2]], y = as.numeric(test_data$Group))
  explainer_3 <- explain_tidymodels(final_wf3, data = test_data[, colnames(test_data) != "Group"], label = mod[rankdata$Row.names[3]], y = as.numeric(test_data$Group))
  
  # Calculate variable importance for each model
  set.seed(1803)
  vip_1 <- model_parts(explainer_1, loss_function = loss_root_mean_square)
  set.seed(1804)
  vip_2 <- model_parts(explainer_2, loss_function = loss_root_mean_square)
  set.seed(1805)
  vip_3 <- model_parts(explainer_3, loss_function = loss_root_mean_square)
  
  # Plot and save the importance plot
  setwd("./survival/")
  p <- ggplot_imp(vip_1, vip_2, vip_3) +  theme_bw()
  ggsave(paste("./vip3-tune", name1, ".pdf", sep = "-"), p, height = 15, width = 15)
  
  # Summarize and rank the variable importance
  summary_1 <- aggregate(dropout_loss ~ variable, data = vip_1, FUN = sum)
  summary_2 <- aggregate(dropout_loss ~ variable, data = vip_2, FUN = sum)
  summary_3 <- aggregate(dropout_loss ~ variable, data = vip_3, FUN = sum)
  
  merged_df_vip <- merge(merge(summary_1, summary_2, by = "variable"), summary_3, by = "variable")
  merged_df_vip$sum <- merged_df_vip[, 2] * 1.0 + merged_df_vip[, 3] * 0.8 + merged_df_vip[, 4] * 0.6
  sorted_vip <- merged_df_vip[order(merged_df_vip$sum, decreasing = TRUE), ]
  write.csv(sorted_vip, file = paste0(paste("./vip3-tune", name1, sep = "-"), ".csv"))
  
  return(sorted_vip)
}

# Run the model ranking function
sorted_vip <- expmod(rank, name)
