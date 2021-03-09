## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Bioconductor, eval = FALSE------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("https://bioconductor.org/biocLite.R")
#  biocLite()

## ----installation, eval = FALSE------------------------------------------
#  # Set one level up methylearning folder.
#  setwd("path/to/folder/containing/methylearning/folder")
#  devtools::install("methylearning")

## ----loading-------------------------------------------------------------
library(methylearning)

## ----loading_dataset, eval = FALSE---------------------------------------
#  if (require(GEOquery) &
#      require(Biobase)) {
#    # Load methylation array experiment from GEO.
#    gse <- getGEO("GSE69229", getGPL = FALSE)
#  }
#  # Create an ml_data object. The column with the sample labels to be used as
#  # classification labels must be specified. In this particular dataset, this
#  # column is named 'outcome:ch1' in 'phenoData' AnnotatedDataFrame.
#  set.seed(1234)
#  ml <- get_GEO_methylarray(gse, target_name = "outcome:ch1",
#                            random_filter = 1000)
#  summary(ml)

## ----load_data_real------------------------------------------------------
# Load the object.
load("demo_data/demo_1.RData")
ml <- ml_data(df, labels_column = dim(df)[2])
summary(ml)

## ----data_split----------------------------------------------------------
# Splitting data to 2/3 for training and 1/3 for test (this is by default).
set.seed(1234)
ml_s <- split_data(ml, splitting = 2/3)
summary(ml_s)

## ----list_fs_methods-----------------------------------------------------
# Print all the feature selection methods currently implemented. To
# return a list with the method names, set 'verbose = FALSE' (default).
available_fs_methods(verbose = TRUE)

## ----feature_selection---------------------------------------------------
set.seed(1234)
fs_methods_to_use <- c("random_fs", "anova_fs", "limma_fs", 
                       "information_gain_fs", "boruta_fs")
ml_f <- ml_fs(ml_s, fs_methods = fs_methods_to_use, selection_size = 20, 
              cores = 3)

## ----fs_exploration------------------------------------------------------
# Summary of the ml_f object.
summary(ml_f)
# Print selection methods used.
selection_methods(ml_f)
# Print all the feature selections obtained.
selection_results(ml_f)
# Plot a Venn diagram for up to 5 feature selection methods to see the overlap.
meths <- c("anova_fs", "limma_fs")
fill <- c("red", "blue")
plot_venn_fs(ml_f, fs_methods = meths, fill = fill, category = meths)
meths <- c("anova_fs", "limma_fs", "boruta_fs")
fill <- c("red", "blue", "orange")
plot_venn_fs(ml_f, fs_methods = meths, fill = fill, category = meths)
# Create a data.frame with feature selection results.
fs_results <- selection_df(ml_f)
fs_results
# Plot the computation time elapsed in each feature selection.
plot_fs_time(ml_f)

## ----classification------------------------------------------------------
# By default, all implemented classification methods will be applied to all 
# feature sets. We can select specific methods using cls_methods and fs_methods
# parameters. Also by default, 10-fold 3x repeated cross-validation is 
# implemented to select best models based on Kappa statistic. As data was
# previously partitioned, test evaluation is also performed when 
# test_eval = TRUE. Argument 'cores' is set to 1 to keep memory usage low.
set.seed(1234)
cls_methods_to_use <- c("knn", "svmRadial", "rf", "lda")
ml_c <- ml_cls(ml_f, cls_methods = cls_methods_to_use, test_eval = TRUE, 
               cores = 1)

## ----cls_evaluation------------------------------------------------------
# Classification summary.
summary(ml_c)
# Get training results.
print(get_training_results(ml_c))
# Get test results.
print(get_test_results(ml_c))
# Get best model based on training data validation.
get_best_model(ml_c)
# Get best feature set based on training data validation.
get_best_feature_set(ml_c)
# Get best model based on test data evaluation.
get_best_model_test(ml_c)
# Get best feature set based on test data evaluation.
get_best_feature_set_test(ml_c)
# Get the caret::confusionMatrix for selected prediction on test data.
get_confusionMatrix_test(ml_c, fs_method = "boruta_fs", cls_method = "knn")
# Plot mean accuracy (with sd) for the best tune from training data validation.
plot(ml_c)
# Plot mean Kappa (with sd) for the best tune from training data validation.
plot(ml_c, mode = "Kappa")
# Plot ROC curves of selected best models from training data validation.
plot(ml_c, mode = "ROC", fs_ROC = "random_fs", cls_ROC = "rf")
plot(ml_c, mode = "ROC", fs_ROC = "anova_fs", cls_ROC = "knn")
plot(ml_c, mode = "ROC", fs_ROC = "boruta_fs", cls_ROC = "lda")
plot(ml_c, mode = "ROC", fs_ROC = "information_gain_fs", cls_ROC = "svmRadial")
# Plot best models accuracy on test data evaluation.
plot_test(ml_c, mode = "Accuracy")
# Plot best models Kappa on test data evaluation.
plot_test(ml_c, mode = "Kappa")
# Plot computation time (classification step only).
plot(ml_c, mode = "Time", time_mode = "cls_only")
# Plot total computation time.
plot(ml_c, mode = "Time")

## ----sessionInfo---------------------------------------------------------
sessionInfo()

