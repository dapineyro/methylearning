#===============================================================================
# Implementation of classification methods using caret package.
# Author: David Piñeyro
# License: GPL-3
# Date: 2018-05-06
#===============================================================================
#' Wrapper to apply caret::train function
#'
#' This is a wrapper of the \link[caret]{train} function. It is configured to
#' perform "repeatedcv" method exclusively. To use with
#' \link[methylearning]{ml_cls} constructor function.
#'
#' @param cls_method [character] a \code{character} vector of length 1 with the
#'   classification method to be applied. To use in "method" argument of
#'   \link[caret]{train}.
#' @param feature_set [character] a \code{character} vector with the feature
#'   selection to be used.
#' @param training_data [data.frame] training data set (only features, not
#'   labels). Samples as rows and features as columns.
#' @param training_labels [factor] training class labels.
#' @param cv_folds [numeric] an \code{integer} indicating the k folds to use in
#'   k-fold cross-validation method. To use in "number" argument of
#'   \link[caret]{trainControl}. Default: 10.
#' @param rep [numeric] an \code{integer} indicating the times to repeat the
#'   cross-validation method. To use in "repeats" argument of
#'   \link[caret]{trainControl}. Default: 3.
#' @param tune_length [numeric] an \code{integer} to use in "tuneLength"
#'   argument of \link[caret]{train}. Default: 10.
#' @param metric [character] a \code{character} vector of length 1 that can be
#'   either "Kappa" or "Accuracy", to indicate the metric use to select the
#'   best tune for the classification model. To use in "metric" argument of
#'   \link[caret]{train}. Default: "Kappa".
#'
#' @return an object of \code{train} class.
#'
#' @examples
#' set.seed(1234)
#' tr_dat <- data.frame(a = runif(10),
#'                      b = runif(10),
#'                      c = runif(10),
#'                      d = runif(10))
#' tr_lab <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' model_cv <- apply_cls_method("knn", c("a", "b"), tr_dat, tr_lab)
#' print(model_cv)
#'
#' @importFrom Matrix Matrix
#' @importFrom glmnet glmnet
#' @importFrom MASS lda
#' @importFrom kernlab ksvm
#' @importFrom nnet nnet
#' @importFrom plyr aaply
#' @importFrom randomForest randomForest
#' @importFrom C50 C5.0
#'
#' @export
apply_cls_method <- function(cls_method, feature_set, training_data,
                             training_labels, cv_folds = 10, rep = 3,
                             tune_length = 10, metric = "Kappa") {
  # Reduce data to the selected features.
  t_data <- training_data[feature_set]
  # Merge training_labels as a column.
  t_data$training_labels <- training_labels
  # Set control paramenters.
  ctrl <- caret::trainControl(method = "repeatedcv",
                              number = cv_folds,
                              repeats = rep,
                              classProbs = TRUE,
                              savePredictions = TRUE)
  model_cv <- caret::train(training_labels ~ .,
                           data = t_data,
                           method = cls_method,
                           trControl = ctrl,
                           tuneLength = tune_length,
                           metric = "Kappa")
  return(model_cv)
}

#' Wrapper for predict + confusion matrix
#'
#' This function is a wrapper that performs \link[stats]{predict} to predict
#' test labels from a test \code{data.frame} using a classification model. Then,
#' a \code{confusionMatrix} object is created using
#' \link[caret]{confusionMatrix}.
#'
#' @param train_obj [train] an object of the \code{train} class, returned by
#'   \link[caret]{train} function, or any other model object to be used with
#'   \link[stats]{predict} function.
#' @param test_data [data.frame] data to be used as test data. It must contain
#'   the same rows and columns as the data used to train the model.
#' @param test_labels [factor] sample observed (real) labels.
#'
#' @examples
#' set.seed(1234)
#' tr_dat <- data.frame(a = runif(10),
#'                      b = runif(10),
#'                      c = runif(10),
#'                      d = runif(10))
#' tr_lab <- factor(rbinom(10, 1, 0.5), labels = c("control", "treatment"))
#' # Train the model with this random data, using only 2 features.
#' model_cv <- apply_cls_method("knn", c("a", "b"), tr_dat, tr_lab)
#' # Evaluating in the same training data.
#' cm <- evaluation(model_cv, tr_dat[c("a", "b")], tr_lab)
#' print(cm)
#'
#' @export
evaluation <- function(train_obj, test_data, test_labels) {
  # Predict test labels.
  p <- stats::predict(train_obj, test_data)
  # Create confusion matrix.
  cm <- caret::confusionMatrix(p, test_labels)
  return(cm)
}

#' Available Classification Methods
#'
#' This function generates a vector with the classification methods that have
#' been tested to use with \link[methylearning]{apply_cls_method}.
#'
#' @details While, in principle, all \link[caret]{train} "methods" would be
#'   suitable to use, only tested methods will be allowed, to ensure
#'   compatibility with \link[methylearning]{ml_cls} methods.
#'
#' @export
available_cls_methods <- function() {
  return(c("knn", "C5.0", "rf", "svmLinear", "svmRadial", "nnet", "glmnet",
           "lda"))
}
