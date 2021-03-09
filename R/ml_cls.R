#===============================================================================
# Definition of 'ml_cls' class and its methods.
# Author: David Piñeyro
# License: GPL-3
# Date: 2018-05-06
#===============================================================================
# Golbal variables
globalVariables(c("obs", "p_c"))

#' Create the list of requested features from an ml_fs object
#'
#' This function creates a list of features requested, from the performed
#' feature selection of an ml_fs object. This function is meant to be used by
#' ml_cls constructor.
#'
#' @param ml_fs_obj [ml_fs] an object of the \code{ml_fs} class.
#' @param requested [character] a \code{character} vector with the feature
#'   selection methods to be selected. Possible values: c("none", "all",
#'   "filters", "wrappers", methylearning::selection_methods(ml_fs_obj)).
#'   Values "all", "filters" and "wrappers" can only be indicated alone, or in
#'   combination with "None". Value "none" means to use all the features,
#'   regardless any feature selection method applied.
#'
#' @return A named \code{list} whith each feature selection.
#'
create_feature_list <- function(ml_fs_obj, requested) {
  ### Create a list with all the features to be used. Each element of the list
  # will contain the features selected by a requested method.
  feature_list <- list()
  i <- 1  # feature_list position iterator.
  if ("none" %in% requested) {
    # Select all the features.
    feature_list[[i]] <- ml_fs_obj$feature_names
    names(feature_list)[i] <- "none"
    i <- i + 1
    # Remove "none" from the requested_methods.
    requested <- requested[-which(requested == "none")]
  }
  if ("all" %in% requested) {
    if (length(requested) > 1) {
      stop("Error: not additional fs_methods can be specified when using 'all'")
    } else {
      # Store all feature selections performed.
      s_res <- selection_results(ml_fs_obj)
      j <- 1  # Iterator for s_res.
      while (j <= length(s_res)) {
        feature_list[[i]] <- s_res[[j]]
        names(feature_list)[i] <- names(s_res[j])
        i <- i + 1
        j <- j + 1
      }
      return(feature_list)
    }
  }
  # Check whether requested were available.
  available <- selection_methods(ml_fs_obj)
  availability <- vapply(requested, function(x) x %in% available, logical(1))
  if (length(requested) != sum(availability)) {
    stop(paste("Error: the following methods are not available:",
               names(availability)[availability == FALSE]))
  }
  # Given that all methods were performed in ml_fs object, collect and return
  # the features.
  s_res <- selection_results(ml_fs_obj)  # All feature selections.
  for (r in requested) {
    feature_list[[i]] <- s_res[[r]]
    names(feature_list)[i] <- r
    i <- i + 1
  }
  return(feature_list)
}

#' Creates \code{ml_cls} objects
#'
#' Constructor for \code{ml_cls} objects. This object is created from a
#' \code{ml_fs} object, using the feature selections that it contains to
#' train and validate the indicated classification methods.
#'
#' Object members:
#' \describe{
#' \item{fs_methods [character]}{a character vector with the feature selection
#'   methods applied.}
#' \item{cls_methods [character]}{a character vector with the classification
#'   methods applied.}
#' \item{validation_method [character]}{an string describing the validation
#'   method applied.}
#' \item{metric [character]}{metric used for select the best model during
#'   validation. One of c("Accuracy", "Kappa").}
#' \item{feature_list [list]}{a named list whith each feature selection.}
#' \item{fs_lengths [list]}{a named list with the length of each feature
#'   selection.}
#' \item{original_features [character]}{a character vector with the complete set
#'   of original features.}
#' \item{sample_names [character]}{a character vector with the sample names.}
#' \item{training_labels [factor]}{a factor with the training labels for
#'   classification.}
#' \item{data_partitioned [logical]}{a logical (TRUE or FALSE) indicating
#'   whether data was previously partitioned in training and test datasets.}
#' \item{cls_results [list]}{a list of lists. Each element of the list
#'   corresponds to a feature selection method used that in turn is a list in
#'   which each element is a classification method used. Is the return of
#'   \link[methylearning]{apply_cls_method} on each element.}
#' \item{test_evaluation [list]}{a list of lists. Each element of the list
#'   corresponds to a feature selection method used that in turn is a list in
#'   which each element is a classification method used. Is the return of
#'   \link[methylearning]{evaluation} on each element.}
#' \item{fs_computation_time [list]}{List of numeric vectors, with the
#'   computation time elapsed in each feature selection computation.}
#' \item{cls_computation_time [list]}{a list of lists. Each element of the list
#'   corresponds to a feature selection method used that in turn is a list in
#'   which each element is a classification method used. The values are
#'   numeric vectors with the elapsed computation time for each classification.}
#' }
#'
#' Object methods:
#' [summary]
#' [plot]
#' [plot_test]
#' [get_training_results]
#' [get_test_results]
#' [get_best_model]
#' [get_best_model_test]
#' [get_confusionMatrix_test]
#' [get_best_feature_set]
#' [get_best_feature_set_test]
#'
#' @param ml_fs_obj [ml_fs] input \code{data.frame} with samples as rows and
#'   features as columns. There must be an extra column with class labels.
#' @param cls_methods [character] a \code{character} vector with the
#'   classification methods to be applied. Possible values: c("all",
#'   methylearning::available_cls_methods()). Default: "all".
#' @param fs_methods [character] a \code{character} vector with the feature
#'   selection methods to be selected. Possible values: c("none", "all",
#'   "filters", "wrappers", methylearning::selection_methods(ml_fs_obj)).
#'   Values "all", "filters" and "wrappers" can only be indicated alone, or in
#'   combination with "None". Value "none" means to use all the features,
#'   regardless any feature selection method applied. Default: "all".
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
#' @param test_eval [logical] whether to predict using test data and store the
#'   performance evaluation. Only applicable if data was previously partitioned
#'   in training and test data. Default: FALSE.
#' @param cores [numeric] number of cores to be used in parallel parts.
#'   Default: 1.
#'
#' @return An object of [ml_cls] class.
#'
#' @examples
#' set.seed(1234)
#' \dontrun{
#' if (require(GEOquery) &
#'     require(Biobase)) {
#'   gse <- getGEO("GSE69229")
#'   # Generate ml_data object.
#'   ml <- get_GEO_methylarray(gse = gse, target_name = "outcome:ch1")
#'   # Generate ml_fs object.
#'   ml_f <- ml_fs(ml, fs_methods = "anova_fs", selection_size = 50)
#'   # Generate ml_cls object.
#'   ml_c <- ml_cls(ml_f, cls_methods = "knn")
#'   summary(ml_c)
#' }
#' }
#'
#' @include cls_methods.R
#'
#' @export
ml_cls <- function(ml_fs_obj, cls_methods = "all", fs_methods = "all",
                   cv_folds = 10, rep = 3, tune_length = 10, metric = "Kappa",
                   test_eval = FALSE, cores = 1) {
  # Check if input object is from ml_fs class.
  if (!("ml_fs" %in% class(ml_fs_obj))) {
    stop("ERROR: Input object is not from 'ml_fs' class.")
  }
  # Collect validation method.
  validation_method <- paste(cv_folds, "fold x", rep,
                             "repeated cross-validation, testing", tune_length,
                             "values for each tuning parameter and",
                             "selecting models by", metric, "\n")
  # Collect input data and labels. Training will be on training data set,
  # regardless data was previously partitioned.
  training_data <- ml_fs_obj$training_data
  training_labels <- ml_fs_obj$training_labels
  # Create the list of features to be used, based on fs_methods argument and
  # checking their validity.
  feature_list <- create_feature_list(ml_fs_obj, fs_methods)
  # Check cls_methods for its validity.
  available_cls <- available_cls_methods()
  if ("all" %in% cls_methods) {
    cls_methods <- available_cls
  }
  valid_cls <- vapply(cls_methods, function(x) {
    ifelse(x %in% available_cls, TRUE, FALSE)
    }, logical(1))
  if (length(cls_methods) > sum(valid_cls)) {
    stop(paste("ERROR: Classification method not implemented:",
               cls_methods[!valid_cls]))
  }
  # glmnet method does not work with only 1 feature. To prevent errors, check
  # whether any feature selection is of length 1 and report an error.
  f_lengths <- unlist(lapply(feature_list, length))
  if (sum(f_lengths < 2) & "glmnet" %in% cls_methods) {
    stop(paste("Error: there was detected a feature selection with length < 2",
               "and 'glmnet' classification method was selected. This method",
               "requires at least 2 features to work. Please, remove this",
               "feature selection or do not execute 'glmnet' classifier."))
  }
  # Name cls_methods to get named lists afterwards.
  names(cls_methods) <- cls_methods
  ### Apply classification methods.
  # Check if core number is allowed.
  if (cores < 1) {
    stop("ERROR: execution with less than 1 core is not allowed.")
  } else if (cores > parallel::detectCores()) {
    warning(paste0("The specified number of cores is higher than the system ,",
                   "cores, setting parallel processing to ",
                   parallel::detectCores(), " cores."))
    cores <- parallel::detectCores()
  }
  # Then, apply all the selected classification methods to all the selected
  # feature sets. cls_results will be a list of lists, each element containing
  # a list of classification results for a particular feature set. At the same
  # time, computation time is calculated and returned as a list of the same
  # structure as cls_results.
  cls_results <- list()
  cls_computation_time <- list()
  for (i in 1:length(feature_list)) {
    fs <- feature_list[[i]]
    fs_name <- names(feature_list)[i]
    r <- parallel::mclapply(cls_methods, function(x) {
      c_t <- system.time(
        c_r <- apply_cls_method(x, fs, training_data,
                                           training_labels, cv_folds, rep,
                                           tune_length, metric)
      )
      return(list(c_results = c_r, c_time = as.numeric(c_t["elapsed"])))
      }, mc.cores = cores)
    only_results <- lapply(r, function(x) x[["c_results"]])
    only_time <- lapply(r, function(x) x[["c_time"]])
    cls_results[[fs_name]] <- only_results
    cls_computation_time[[fs_name]] <- only_time
  }
  # Performance evaluation if data was partitioned.
  if (ml_fs_obj$data_partitioned & test_eval) {
    # Create an empty list to store results.
    test_evaluation <- list()
    for (f in names(feature_list)) {
      features <- feature_list[[f]]
      test_evaluation[[f]] <- list()
      for (c in cls_methods) {
        cm <- evaluation(cls_results[[f]][[c]], ml_fs_obj$test_data[features],
                        ml_fs_obj$test_labels)
        test_evaluation[[f]][[c]] <- cm
      }
    }
  } else {
    test_evaluation <- FALSE
  }
  # Constructor.
  structure(list(fs_methods = names(feature_list),
                 cls_methods = cls_methods,
                 validation_method = validation_method,
                 metric = metric,
                 feature_list = feature_list,
                 fs_lengths = lapply(feature_list, length),
                 original_features = ml_fs_obj$feature_names,
                 sample_names = ml_fs_obj$sample_names,
                 training_labels = ml_fs_obj$training_labels,
                 data_partitioned = ml_fs_obj$data_partitioned,
                 cls_results = cls_results,
                 test_evaluation = test_evaluation,
                 fs_computation_time = ml_fs_obj$computation_time,
                 cls_computation_time = cls_computation_time),
            class = "ml_cls")
}

#===============================================================================
# Methods for 'ml_cls' class objects.
#===============================================================================
#' Prints a summary for a \code{ml_cls} object
#'
#' \code{summary} method for objects of the \code{ml_cls} class.
#'
#' @param object [ml_cls_obj] an object of the \code{ml_cls} class.
#' @param ... additional arguments to pass to \code{summary.default}. See
#'   \link[base]{summary}.
#'
#' @return a printed summary.
#'
#' @examples
#' set.seed(1234)
#' \dontrun{
#' if (require(GEOquery) &
#'     require(Biobase)) {
#'   gse <- getGEO("GSE69229")
#'   # Generate ml_data object.
#'   ml <- get_GEO_methylarray(gse = gse, target_name = "outcome:ch1")
#'   # Generate ml_fs object.
#'   ml_f <- ml_fs(ml, fs_methods = "anova_fs", selection_size = 50)
#'   # Generate ml_cls object.
#'   ml_c <- ml_cls(ml_f, cls_methods = "knn")
#'   summary(ml_c)
#' }
#' }
#'
#' @export
summary.ml_cls <- function(object, ...) {
  if (!("ml_cls" %in% class(object))) {
    stop("Error: method 'summary.ml_cls' is only available for 'ml_cls' objects.")
  }
  cat(paste0("Summary for object '", deparse(substitute(object)),
             "' of the 'ml_cls' class\n\n"))

  cat("Feature selection methods performed:\n")
  cat(paste("\t", paste(object$fs_methods, collapse = ", "), "\n\n"))

  cat("Classification methods performed:\n")
  cat(paste("\t", paste(object$cls_methods, collapse = ", "), "\n\n"))

  cat("Full data set dimensions:\n")
  cat(paste("\tNumber of features:", length(object$original_features),
            "\n"))
  cat(paste("\tNumber of training samples:",
            length(object$training_labels), "\n"))
  cat(paste("\tNumber of test samples:",
            (length(object$sample_names) -
             length(object$training_labels)), "\n"))
  cat(paste("\tClass labels:", paste(levels(object$training_labels),
                                         collapse = ", "), "\n\n"))

  cat("Validation method used in training data classification:\n")
  cat(paste("\t", object$validation_method, "\n"))

  cat(paste0("Best model selected by ", object$metric, ":\n"))
  best <- get_best_model(object)
  cat(paste("\tBest feature selection method:", best$best_fs, "\n"))
  cat(paste("\tBest classification method:", best$best_cls, "\n"))
  cat(paste("\tBest accuracy:", best$best_acc, "\n"))
  cat(paste("\tBest Kappa:", best$best_kpp, "\n\n"))

  cat("Computation time:\n")
  fs <- names(object$cls_computation_time)
  for (fs_n in fs) {
    cl <- names(object$cls_computation_time[[fs_n]])
    for (cl_n in cl) {
      cat(paste0("\t", fs_n, " + ", cl_n, ":\n"))
      fs_t <- round(object$fs_computation_time[[fs_n]], 3)
      cl_t <- round(object$cls_computation_time[[fs_n]][[cl_n]], 3)
      cat(paste0("\t\tFeature Selection only: ", fs_t , " seconds.\n"))
      cat(paste0("\t\tClassification only: ", cl_t, " seconds.\n"))
      cat(paste0("\t\tTotal: ", fs_t + cl_t, " seconds.\n"))
    }
  }
  cat("\n")

  cat("Summary of classification results using training data:\n")
  print(get_training_results(object))

  cat(paste("Data partitioned:", object$data_partitioned, "\n"))
  if (object$data_partitioned & (class(object$test_evaluation) ==
                                     "list")) {
    cat("Summary of classification results using test data:\n")
    print(get_test_results(object))
  }
}

#===============================================================================
#' Get Training Results for Accuracy and Kappa
#'
#' Get a list with two \code{data.frame}, one for accuracy and the other for
#' kappa, where rows are feature selection methods and columns are
#' classification methods.
#'
#' @param ml_cls_obj [ml_cls] an object of 'ml_cls' class.
#'
#' @export
get_training_results <- function(ml_cls_obj) {
  UseMethod("get_training_results", ml_cls_obj)
}
#' @export
get_training_results.default <- function(ml_cls_obj) {
  message("get_training_results method is only defined for objects of ml_cls class.\n")
  message("Nothing was done.\n")
  return(ml_cls_obj)
}
#' @export
get_training_results.ml_cls <- function(ml_cls_obj) {
  if (!("ml_cls" %in% class(ml_cls_obj))) {
    stop("Error: method 'get_training_results.ml_cls' is only available for 'ml_cls' objects.")
  }
  # In each data.frame to return, feature selection methods will be rows and
  # classification methods will be columns.
  c_names <- ml_cls_obj$cls_methods
  r_names <- ml_cls_obj$fs_methods
  # Initialize and empty data.frame for each measure to return: "Accuracy" and
  # "Kappa".
  acc <- data.frame(matrix(NA, nrow = length(r_names), ncol = length(c_names),
                           dimnames = list(r_names, c_names)))
  kpp <- acc
  # Filling up accuracy and kappa data.frames.
  for (f in r_names) {
    for (c in c_names) {
      acc[f, c] <- max(ml_cls_obj$cls_results[[f]][[c]]$results$Accuracy)
      kpp[f, c] <- max(ml_cls_obj$cls_results[[f]][[c]]$results$Kappa)
    }
  }
  return(list(Accuracy = acc, Kappa = kpp))
}

#===============================================================================
#' Get Test Results for Accuracy and Kappa
#'
#' Get a list with two \code{data.frame}, one for accuracy and the other for
#' kappa, where rows are feature selection methods and columns are
#' classification methods.
#'
#' @param ml_cls_obj [ml_cls] an object of 'ml_cls' class.
#'
#' @return a \code{list} with two elements named "Accuracy" and "Kappa", each
#'   containing a \code{data.frame} with the obtained values, where rows are
#'   feature selection methods and columns are classification methods.
#'
#' @export
get_test_results <- function(ml_cls_obj) {
  UseMethod("get_test_results", ml_cls_obj)
}
#' @export
get_test_results.default <- function(ml_cls_obj) {
  message("get_test_results method is only defined for objects of ml_cls class.\n")
  message("Nothing was done.\n")
  return(ml_cls_obj)
}
#' @export
get_test_results.ml_cls <- function(ml_cls_obj) {
  if (!("ml_cls" %in% class(ml_cls_obj))) {
    stop("Error: method 'get_test_results.ml_cls' is only available for 'ml_cls' objects.")
  }
  # Check whether data was partitioned.
  if (!(ml_cls_obj$data_partitioned)) {
    stop("Error: method 'get_test_results.ml_cls' is only available for 'ml_cls' objects with partitioned data.")
  }
  # In each data.frame to return, feature selection methods will be rows and
  # classification methods will be columns.
  c_names <- ml_cls_obj$cls_methods
  r_names <- ml_cls_obj$fs_methods
  # Initialize and empty data.frame for each measure to return: "Accuracy" and
  # "Kappa".
  acc <- data.frame(matrix(NA, nrow = length(r_names), ncol = length(c_names),
                           dimnames = list(r_names, c_names)))
  kpp <- acc
  # Filling up accuracy and kappa data.frames.
  for (f in r_names) {
    for (c in c_names) {
      acc[f, c] <- max(ml_cls_obj$test_evaluation[[f]][[c]]$overall["Accuracy"])
      kpp[f, c] <- max(ml_cls_obj$test_evaluation[[f]][[c]]$overall["Kappa"])
    }
  }
  return(list(Accuracy = acc, Kappa = kpp))
}

#===============================================================================
#' Get Best Model
#'
#' This function selects the best model based on training data validation and
#' metric used in \link[methylearning]{ml_cls} function and return a list with
#' some objects (see Value).
#'
#' @param ml_cls_obj [ml_cls] an object of 'ml_cls' class.
#'
#' @return a \code{list} with the following elements:
#' \describe{
#'   \item{best_model [train]}{a \code{train} object returned by
#'     \link[caret]{train}.}
#'   \item{best_fs [character]}{best feature selection method.}
#'   \item{best_cls [character]}{best classification method.}
#'   \item{best_acc [numeric]}{best accuracy obtained.}
#'   \item{best_kpp [numeric]}{best Kappa obtained.}
#' }
#'
#' @export
get_best_model <- function(ml_cls_obj) {
  UseMethod("get_best_model", ml_cls_obj)
}
#' @export
get_best_model.default <- function(ml_cls_obj) {
  message("get_best_model method is only defined for objects of ml_cls class.\n")
  message("Nothing was done.\n")
  return(ml_cls_obj)
}
#' @export
get_best_model.ml_cls <- function(ml_cls_obj) {
  if (!("ml_cls" %in% class(ml_cls_obj))) {
    stop("Error: method 'get_best_model.ml_cls' is only available for 'ml_cls' objects.")
  }
  # First, decide which is the best model. ml_cls_obj$metric will be used as
  # metric.
  # list of data.frame with accuracy and kappa training results.
  res <- get_training_results(ml_cls_obj)
  # names of the best. Only firs result selected.
  best_first <- which(res[[ml_cls_obj$metric]] == max(res[[ml_cls_obj$metric]]),
                      arr.ind = TRUE)[1, ]
  best_fs <- rownames(res[[ml_cls_obj$metric]])[best_first[1]]
  best_cls <- colnames(res[[ml_cls_obj$metric]])[best_first[2]]
  # Get best model (train object).
  best_model <- ml_cls_obj$cls_results[[best_fs]][[best_cls]]
  best_model_acc <- res[["Accuracy"]][best_fs, best_cls]
  best_model_kpp <- res[["Kappa"]][best_fs, best_cls]
  return(list(best_model = best_model,
              best_fs = best_fs,
              best_cls = best_cls,
              best_acc = best_model_acc,
              best_kpp = best_model_kpp))
}

#===============================================================================
#' Get Best Model From Test Evaluation
#'
#' This function selects the best model based on test data evaluation, if
#' available. It returns a list with some objects (see Value).
#'
#' @param ml_cls_obj [ml_cls] an object of 'ml_cls' class.
#'
#' @return a \code{list} with the following elements:
#' \describe{
#'   \item{best_model [train]}{a \code{train} object returned by
#'     \link[caret]{train}.}
#'   \item{best_fs [character]}{best feature selection method.}
#'   \item{best_cls [character]}{best classification method.}
#'   \item{best_acc [numeric]}{best accuracy obtained.}
#'   \item{best_kpp [numeric]}{best Kappa obtained.}
#' }
#'
#' @export
get_best_model_test <- function(ml_cls_obj) {
  UseMethod("get_best_model_test", ml_cls_obj)
}
#' @export
get_best_model_test.default <- function(ml_cls_obj) {
  message("get_best_model_test method is only defined for objects of ml_cls class.\n")
  message("Nothing was done.\n")
  return(ml_cls_obj)
}
#' @export
get_best_model_test.ml_cls <- function(ml_cls_obj) {
  if (!("ml_cls" %in% class(ml_cls_obj))) {
    stop("Error: method 'get_best_model_test.ml_cls' is only available for 'ml_cls' objects.")
  }
  # Check whether data was partitioned.
  if (!(ml_cls_obj$data_partitioned)) {
    stop("Error: method 'get_best_model_test.ml_cls' is only available for 'ml_cls' objects with partitioned data.")
  }
  # First, decide which is the best model. ml_cls_obj$metric will be used as
  # metric.
  # list of data.frame with accuracy and kappa test results.
  res <- get_test_results(ml_cls_obj)
  # names of the best. Only firs result selected.
  best_first <- which(res[[ml_cls_obj$metric]] == max(res[[ml_cls_obj$metric]]),
                      arr.ind = TRUE)[1, ]
  best_fs <- rownames(res[[ml_cls_obj$metric]])[best_first[1]]
  best_cls <- colnames(res[[ml_cls_obj$metric]])[best_first[2]]
  # Get best model (train object).
  best_model <- ml_cls_obj$cls_results[[best_fs]][[best_cls]]
  best_model_acc <- res[["Accuracy"]][best_fs, best_cls]
  best_model_kpp <- res[["Kappa"]][best_fs, best_cls]
  return(list(best_model = best_model,
              best_fs = best_fs,
              best_cls = best_cls,
              best_acc = best_model_acc,
              best_kpp = best_model_kpp))
}

#===============================================================================
#' Get Confusion Matrix From Test Evaluation
#'
#' This function returns a \code{confusionMatrix} object from the selected
#' combination of \code{fs_method} and \code{cls_method}.
#'
#' @param ml_cls_obj [ml_cls] an object of 'ml_cls' class.
#' @param ... additional arguments for \link[methylearning]{get_confusionMatrix_test.ml_cls}
#'   function.
#'
#' @return a \code{confusionMatrix} object, generated by
#'   \link[caret]{confusionMatrix}
#'
#' @export
get_confusionMatrix_test <- function(ml_cls_obj, ...) {
  UseMethod("get_confusionMatrix_test")
}
#' @export
get_confusionMatrix_test.default <- function(ml_cls_obj, ...) {
  message("get_confusionMatrix_test method is only defined for objects of ml_cls class.\n")
  message("Nothing was done.\n")
  return(ml_cls_obj)
}
#' Get Confusion Matrix From Test Evaluation
#'
#' This function returns a \code{confusionMatrix} object from the selected
#' combination of \code{fs_method} and \code{cls_method}.
#'
#' @param ml_cls_obj [ml_cls] an object of 'ml_cls' class.
#' @param fs_method [character] an applied feature selection method.
#' @param cls_method [character] an applied classification method.
#' @param ... additional arguments. None available yet...
#'
#' @export
get_confusionMatrix_test.ml_cls <- function(ml_cls_obj, fs_method, cls_method,
                                            ...) {
  if (!("ml_cls" %in% class(ml_cls_obj))) {
    stop("Error: method 'get_confusionMatrix_test.ml_cls' is only available for 'ml_cls' objects.")
  }
  # Check whether data was partitioned.
  if (!(ml_cls_obj$data_partitioned)) {
    stop("Error: method 'get_confusionMatrix_test.ml_cls' is only available for 'ml_cls' objects with partitioned data.")
  }
  # Check whether fs_method and cls_method are available.
  if (!(fs_method %in% ml_cls_obj$fs_methods)){
    stop(paste("Error:", fs_method, "was not applied to this data set."))
  }
  if (!(cls_method %in% ml_cls_obj$cls_methods)){
    stop(paste("Error:", cls_method, "was not applied to this data set."))
  }
  return(ml_cls_obj$test_evaluation[[fs_method]][[cls_method]])
}

#===============================================================================
#' Get Best Performant Feature Set Based on Training Validation
#'
#' This function returns the best performant feature set based on training data
#' validation and metric used in \link[methylearning]{ml_cls} function
#'
#' @param ml_cls_obj [ml_cls] an object of 'ml_cls' class.
#'
#' @return a \code{character} vector with the best performant feature set.
#'
#' @export
get_best_feature_set <- function(ml_cls_obj) {
  UseMethod("get_best_feature_set", ml_cls_obj)
}
#' @export
get_best_feature_set.default <- function(ml_cls_obj) {
  message("get_best_feature_set method is only defined for objects of ml_cls class.\n")
  message("Nothing was done.\n")
  return(ml_cls_obj)
}
#' @export
get_best_feature_set.ml_cls <- function(ml_cls_obj) {
  if (!("ml_cls" %in% class(ml_cls_obj))) {
    stop("Error: method 'get_best_feature_set.ml_cls' is only available for 'ml_cls' objects.")
  }
  best_fs <- get_best_model(ml_cls_obj)$best_fs
  return(ml_cls_obj$feature_list[[best_fs]])
}

#===============================================================================
#' Get Best Performant Feature Set Based on Test Evaluation
#'
#' This function returns the best performant feature set based on test data
#' evaluation and metric used in \link[methylearning]{ml_cls} function
#'
#' @param ml_cls_obj [ml_cls] an object of 'ml_cls' class.
#'
#' @return a \code{character} vector with the best performant feature set.
#'
#' @export
get_best_feature_set_test <- function(ml_cls_obj) {
  UseMethod("get_best_feature_set_test", ml_cls_obj)
}
#' @export
get_best_feature_set_test.default <- function(ml_cls_obj) {
  message("get_best_feature_set_test method is only defined for objects of ml_cls class.\n")
  message("Nothing was done.\n")
  return(ml_cls_obj)
}
#' @export
get_best_feature_set_test.ml_cls <- function(ml_cls_obj) {
  if (!("ml_cls" %in% class(ml_cls_obj))) {
    stop("Error: method 'get_best_feature_set_test.ml_cls' is only available for 'ml_cls' objects.")
  }
  # Check whether data was partitioned.
  if (!(ml_cls_obj$data_partitioned)) {
    stop("Error: method 'get_best_feature_set_test.ml_cls' is only available for 'ml_cls' objects with partitioned data.")
  }
  best_fs <- get_best_model_test(ml_cls_obj)$best_fs
  return(ml_cls_obj$feature_list[[best_fs]])
}

#===============================================================================
#' Plot function for 'ml_cls' objects
#'
#' This function plot some graphics from a 'ml_cls' object, based on training
#' validation.
#'
#' @details This function can draw four type of graphics all of them based
#'   on training validation:
#'   \describe{
#'     \item{Accuracy barplot}{Barplot with accuracy values and accuracy sd.}
#'     \item{Kappa barplot}{Barplot with Kappa values and Kappa sd.}
#'     \item{ROC curves}{ROC curves and AUC values. For two-classes
#'       classification problems.}
#'     \item{Computation time barplot}{Barplot with the computation time took
#'       by classification or by feature selection + classification.}
#'   }
#'
#' @usage \method{plot}{ml_cls}(x, mode = "Accuracy", fs_ROC = "", cls_ROC = "", time_mode = "sum", ...)
#'
#' @param x [ml_cls] an object of 'ml_cls' class.
#' @param mode [character] one of the following strings: c("Accuracy", "Kappa",
#'   "ROC", "Time") to indicate the type of graphic to draw. Default: "Accuracy".
#' @param fs_ROC [character] a feature selection method previously performed, to
#'   construct ROC curve. When mode = "ROC", this should be a valid method.
#' @param cls_ROC [character] a classification method previously performed, to
#'   construct ROC curve. When mode = "ROC", this should be a valid method.
#' @param time_mode [character] whether to plot classification times only
#'   (time_mode = 'cls_only') or to plot the sum of feature selection and
#'   classification times (time_mode = 'sum'). Default: 'sum'.
#' @param ... additional arguments to pass to \link[graphics]{plot}.
#'
#' @return a plot.
#'
#' @export
plot.ml_cls <- function(x, mode = "Accuracy", fs_ROC = "",
                        cls_ROC = "", time_mode = "sum", ...) {
  if (!("ml_cls" %in% class(x))) {
    stop("Error: method 'plot.ml_cls' is only available for 'ml_cls' objects.")
  }
  # Check that is two-class classification for "ROC" mode.
  if (mode == "ROC") {
    if (nlevels(x$training_labels) != 2) {
      stop("Error: ROC curves are only available for two-class datasets.")
    }
    # Check correct fs_ROC and cls_ROC.
    if (!(fs_ROC %in% x$fs_methods)) {
      stop("Error: fs_ROC feature selection method not performed.")
    }
    if (!(cls_ROC %in% x$cls_methods)) {
      stop("Error: cls_ROC classification method not performed.")
    }
    if (length(fs_ROC) != 1) {
      stop("Error: only one method can be indicated in fs_ROC argument.")
    }
    if (length(cls_ROC) != 1) {
      stop("Error: only one method can be indicated in cls_ROC argument.")
    }
  }
  # Check correct mode input.
  mode_available <- c("Accuracy", "Kappa", "ROC", "Time")
  if (length(mode) > 1 | !(mode %in% mode_available)) {
    stop(paste("Error: incorrect plot mode indicated. Use one of",
               paste(mode_available, collapse = ", ")))
  }
  if (mode == "Accuracy") {
    # Best accuracy values data.frame.
    acc <- get_training_results(x)$Accuracy
    # Get the same data.frame but for the sd.
    # In each data.frame to return, feature selection methods will be rows and
    # classification methods will be columns.
    c_names <- x$cls_methods
    r_names <- x$fs_methods
    # Initialize and empty data.frame.
    acc_sd <- data.frame(matrix(NA, nrow = length(r_names), ncol = length(c_names),
                                dimnames = list(r_names, c_names)))
    # Filling up acc_sd.
    for (f in r_names) {
      for (c in c_names) {
        idx <- which(x$cls_results[[f]][[c]]$results$Accuracy == max(
          x$cls_results[[f]][[c]]$results$Accuracy))[1]
        acc_sd[f, c] <- x$cls_results[[f]][[c]]$results$AccuracySD[idx]
      }
    }
    ### Plot a grouped barchart with ggplot2.
    # Transform each data.frame from n_fs x n_cls to all_values x
    # (fs + cls + acc).
    # All values as a vectors.
    accuracy <- unlist(acc)
    accuracy_sd <- unlist(acc_sd)
    fs_method <- as.factor(rep(rownames(acc), times = ncol(acc)))
    cls_method <- as.factor(rep(colnames(acc), each = nrow(acc)))
    df2plot <- data.frame(fs_method, cls_method, accuracy, accuracy_sd)
    # Plot
    ggplot2::ggplot(df2plot,
                    ggplot2::aes(x=cls_method, y=accuracy, fill=fs_method)) +
      ggplot2::geom_bar(position=ggplot2::position_dodge(),
                        stat="identity",
                        colour="black", # Use black outlines,
                        size=.3) +      # Thinner lines
      ggplot2::geom_errorbar(ggplot2::aes(ymin=accuracy-accuracy_sd,
                                          ymax=accuracy+accuracy_sd),
                             size=.3,    # Thinner lines
                             width=.2,
                             position=ggplot2::position_dodge(.9)) +
      ggplot2::xlab("Classification Method") +
      ggplot2::ylab(mode) +
      ggplot2::scale_fill_hue(name="Feature selection", # Legend label, use darker colors
                              breaks=rownames(acc),
                              labels=rownames(acc)) +
      ggplot2::ggtitle("Classification Performance on Training Data") +
      ggplot2::scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1)) +
      ggplot2::theme_bw()
  } else if (mode == "Kappa") {
    # Best kappa values data.frame.
    kpp <- get_training_results(x)$Kappa
    # Get the same data.frame but for the sd.
    # In each data.frame to return, feature selection methods will be rows and
    # classification methods will be columns.
    c_names <-x$cls_methods
    r_names <-x$fs_methods
    # Initialize and empty data.frame.
    kpp_sd <- data.frame(matrix(NA, nrow = length(r_names), ncol = length(c_names),
                                dimnames = list(r_names, c_names)))
    # Filling up kpp_sd.
    for (f in r_names) {
      for (c in c_names) {
        idx <- which(x$cls_results[[f]][[c]]$results$Kappa == max(
          x$cls_results[[f]][[c]]$results$Kappa))[1]
        kpp_sd[f, c] <- x$cls_results[[f]][[c]]$results$KappaSD[idx]
      }
    }
    ### Plot a grouped barchart with ggplot2.
    # Transform each data.frame from n_fs x n_cls to all_values x
    # (fs + cls + kpp).
    # All values as a vectors.
    kappa <- unlist(kpp)
    kappa_sd <- unlist(kpp_sd)
    fs_method <- as.factor(rep(rownames(kpp), times = ncol(kpp)))
    cls_method <- as.factor(rep(colnames(kpp), each = nrow(kpp)))
    df2plot <- data.frame(fs_method, cls_method, kappa, kappa_sd)
    # Plot
    ggplot2::ggplot(df2plot,
                    ggplot2::aes(x=cls_method, y=kappa, fill=fs_method)) +
      ggplot2::geom_bar(position=ggplot2::position_dodge(),
                        stat="identity",
                        colour="black", # Use black outlines,
                        size=.3) +      # Thinner lines
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin=kappa-kappa_sd, ymax=kappa+kappa_sd),
        size=.3,    # Thinner lines
        width=.2,
        position=ggplot2::position_dodge(.9)) +
      ggplot2::xlab("Classification Method") +
      ggplot2::ylab(mode) +
      ggplot2::scale_fill_hue(name="Feature selection", # Legend label, use darker colors
                              breaks=rownames(kpp),
                              labels=rownames(kpp)) +
      ggplot2::ggtitle("Classification Performance on Training Data") +
      ggplot2::scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1)) +
      ggplot2::theme_bw()
  } else if (mode == "Time") {
    # Check correct time_mode.
    time_mode_available <- c("cls_only", "sum")
    if ((length(time_mode) > 1) | !(time_mode %in% time_mode_available)) {
      stop(paste("Error: incorrect time_mode indicated. Use one of",
                 paste(time_mode_available, collapse = ", ")))
    }
    if (time_mode == "cls_only") {
      feature_selection <- character(0)
      classification <- character(0)
      elapsed_time <- numeric(0)
      fs <- names(x$cls_computation_time)
      for (fs_n in fs) {
        cl <- names(x$cls_computation_time[[fs_n]])
        for (cl_n in cl) {
          elapsed_time <- c(elapsed_time,
                     round(x$cls_computation_time[[fs_n]][[cl_n]], 3))
          feature_selection <- c(feature_selection, fs_n)
          classification <- c(classification, cl_n)
        }
      }
      df2plot <- data.frame(as.factor(feature_selection),
                            as.factor(classification), elapsed_time)
      # Plot
      ggplot2::ggplot(df2plot,
                      ggplot2::aes(x=classification, y=elapsed_time,
                                   fill=feature_selection)) +
        ggplot2::geom_bar(position=ggplot2::position_dodge(),
                          stat="identity",
                          colour="black", # Use black outlines,
                          size=.3) +      # Thinner lines
        ggplot2::xlab("Classification methods") +
        ggplot2::ylab("Elapsed time (seconds)") +
        ggplot2::scale_fill_hue(name="FS methods") +
        ggplot2::ggtitle("Computation Time Classification Only") +
        ggplot2::theme_bw()
    } else if (time_mode == "sum") {
      feature_selection <- character(0)
      classification <- character(0)
      elapsed_time <- numeric(0)
      fs <- names(x$cls_computation_time)
      for (fs_n in fs) {
        cl <- names(x$cls_computation_time[[fs_n]])
        for (cl_n in cl) {
          sum_f_c <- (x$cls_computation_time[[fs_n]][[cl_n]] +
                        x$fs_computation_time[[fs_n]])
          elapsed_time <- c(elapsed_time, round(sum_f_c, 3))
          feature_selection <- c(feature_selection, fs_n)
          classification <- c(classification, cl_n)
        }
      }
      df2plot <- data.frame(as.factor(feature_selection),
                            as.factor(classification), elapsed_time)
      # Plot
      ggplot2::ggplot(df2plot,
                      ggplot2::aes(x=classification, y=elapsed_time,
                                   fill=feature_selection)) +
        ggplot2::geom_bar(position=ggplot2::position_dodge(),
                          stat="identity",
                          colour="black", # Use black outlines,
                          size=.3) +      # Thinner lines
        ggplot2::xlab("Classification methods") +
        ggplot2::ylab("Elapsed time (seconds)") +
        ggplot2::scale_fill_hue(name="FS methods") +
        ggplot2::ggtitle("Feature Selection + Classification Total Computation Time") +
        ggplot2::theme_bw()
    }
  }  else if (mode == "ROC") {
    # Set positive class as the first class.
    pos_cls <- levels(x$training_labels)[2]
    # Plot only one ROC curve, for each selected combination of fs_method +
    # cls_method.
    f <- fs_ROC
    c <- cls_ROC
     # Use predicted values of the best tune.
    bt <- x$cls_results[[f]][[c]]$bestTune
    names_bt <- names(bt)
    pred <- x$cls_results[[f]][[c]]$pred
    # For many tune parameters.
    for (i in 1:length(bt)) {
      idx <- pred[,names_bt[i]] %in% bt[1, i]
      pred <- pred[idx, ]
    }
    # Change pos_cls column name to p_c, to always use the same in plots.
    names(pred)[names(pred) == pos_cls] <- "p_c"
    # Plot ROC curve and AUC value.
    roc2plot <- ggplot2::ggplot(pred, ggplot2::aes(m = p_c, d = obs)) +
                  plotROC::geom_roc() +
                  ggplot2::coord_equal() +
                  plotROC::style_roc() +
                  ggplot2::ggtitle(paste("ROC curve for", f, "selection +", c,
                                   "classification")) +
                  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    roc2plot <- roc2plot + ggplot2::annotate("text", x = .75, y = .25,
              label = paste("AUC =", round(plotROC::calc_auc(roc2plot)$AUC, 3)))
    return(roc2plot)
  }
}

#===============================================================================
#' Specialized plot function for test data evaluation
#'
#' This function plot 'Accuracy' or 'Kappa' barcharts from an 'ml_cls' object,
#' based on test data evaluation.
#'
#' @details This function can draw two types of graphics all of them based
#'   on test data evaluation:
#'   \describe{
#'     \item{Accuracy barplot}{Barplot with accuracy values.}
#'     \item{Kappa barplot}{Barplot with Kappa values.}
#'   }
#'
#' @param ml_cls_obj [ml_cls] an object of 'ml_cls' class.#'
#' @param ... additional arguments for \link[methylearning]{plot_test.ml_cls}.
#'
#' @return a plot.
#'
#' @export
plot_test <- function(ml_cls_obj, ...) {
UseMethod("plot_test")
}
#' @export
plot_test.default <- function(ml_cls_obj, ...) {
  message("plot_test is only defined for objects of ml_cls class.\n")
  message("Nothing was done.\n")
  return(ml_cls_obj)
}
#' Specialized plot function for test data evaluation
#'
#' This function plot 'Accuracy' or 'Kappa' barcharts from an 'ml_cls' object,
#' based on test data evaluation.
#'
#' @details This function can draw two types of graphics all of them based
#'   on test data evaluation:
#'   \describe{
#'     \item{Accuracy barplot}{Barplot with accuracy values.}
#'     \item{Kappa barplot}{Barplot with Kappa values.}
#'   }
#'
#' @param ml_cls_obj [ml_cls] an object of 'ml_cls' class.
#' @param mode [character] one of the following strings: c("Accuracy", "Kappa")
#'   to indicate the type of graphic to draw.
#' @param ... additional arguments. None available yet.
#'
#' @export
plot_test.ml_cls <- function(ml_cls_obj, mode = "Accuracy", ...) {
  if (!("ml_cls" %in% class(ml_cls_obj))) {
    stop("Error: method 'plot_test.ml_cls' is only available for 'ml_cls' objects.")
  }
  # Check whether data was partitioned.
  if (!(ml_cls_obj$data_partitioned)) {
    stop("Error: method 'plot_test.ml_cls' is only available for 'ml_cls' objects with partitioned data.")
  }
  # Check correct mode input.
  mode_available <- c("Accuracy", "Kappa")
  if (length(mode) > 1 | !(mode %in% mode_available)) {
    stop(paste("Error: incorrect plot mode indicated. Use one of",
               paste(mode_available, collapse = ", ")))
  }
  if (mode == "Accuracy") {
    # Best accuracy values data.frame.
    acc <- get_test_results(ml_cls_obj)$Accuracy
    ### Plot a grouped barchart with ggplot2.
    # Transform acc data.frame from n_fs x n_cls to all_values x
    # (fs + cls + acc).
    # All values as a vectors.
    accuracy <- unlist(acc)
    fs_method <- as.factor(rep(rownames(acc), times = ncol(acc)))
    cls_method <- as.factor(rep(colnames(acc), each = nrow(acc)))
    df2plot <- data.frame(fs_method, cls_method, accuracy)
    # Plot
    ggplot2::ggplot(df2plot,
                    ggplot2::aes(x=cls_method, y=accuracy, fill=fs_method)) +
      ggplot2::geom_bar(position=ggplot2::position_dodge(),
                        stat="identity",
                        colour="black", # Use black outlines,
                        size=.3) +      # Thinner lines
      ggplot2::xlab("Classification Method") +
      ggplot2::ylab(mode) +
      ggplot2::scale_fill_hue(name="Feature selection", # Legend label, use darker colors
                              breaks=rownames(acc),
                              labels=rownames(acc)) +
      ggplot2::ggtitle("Classification Performance on Test Data") +
      ggplot2::scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1)) +
      ggplot2::theme_bw()
  } else if (mode == "Kappa") {
    # Best kappa values data.frame.
    kpp <- get_test_results(ml_cls_obj)$Kappa
    ### Plot a grouped barchart with ggplot2.
    # Transform each data.frame from n_fs x n_cls to all_values x
    # (fs + cls + kpp).
    # All values as a vectors.
    kappa <- unlist(kpp)
    fs_method <- as.factor(rep(rownames(kpp), times = ncol(kpp)))
    cls_method <- as.factor(rep(colnames(kpp), each = nrow(kpp)))
    df2plot <- data.frame(fs_method, cls_method, kappa)
    # Plot
    ggplot2::ggplot(df2plot,
                    ggplot2::aes(x=cls_method, y=kappa, fill=fs_method)) +
      ggplot2::geom_bar(position=ggplot2::position_dodge(),
                        stat="identity",
                        colour="black", # Use black outlines,
                        size=.3) +      # Thinner lines
      ggplot2::xlab("Classification Method") +
      ggplot2::ylab(mode) +
      ggplot2::scale_fill_hue(name="Feature selection", # Legend label, use darker colors
                              breaks=rownames(kpp),
                              labels=rownames(kpp)) +
      ggplot2::ggtitle("Classification Performance on Test Data") +
      ggplot2::scale_y_continuous(breaks=c(0.0, 0.2, 0.4, 0.6, 0.8, 1)) +
      ggplot2::theme_bw()
  }
}
