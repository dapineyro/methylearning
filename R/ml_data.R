#===============================================================================
# Definition of 'ml_data' class and its methods.
# Author: David Piñeyro
# License: GPL-3
# Date: 2018-04-21
#===============================================================================

#' Creates \code{ml_data} objects
#'
#' Constructor for \code{ml_data} objects. This object encapsulates input data
#' to be used with rest of package functions.
#'
#' Object members:
#' \describe{
#' \item{classes [character]}{class labels present in the data set.}
#' \item{feature_types [character]}{object classes of feature variables.}
#' \item{sample_names [character]}{all samples names. They are always
#'   \code{input data.frame} \code{rownames}.}
#' \item{training_data [data.frame]}{training data set (only features, not
#'   labels). Samples as rows and features as columns.}
#' \item{training_labels [factor]}{training class labels.}
#' \item{training_sample_names [character]}{training sample names.}
#' \item{test_data [data.frame]}{test data set (only features). It is only
#'   present when \code{data.partitioned} is TRUE, otherwise NULL.}
#' \item{test_labels [factor]}{test class labels. It is only present when
#'   \code{data.partitioned} is TRUE, otherwise NULL.}
#' \item{test_sample_labels [character]}{test sample names. It is only
#'   present when \code{data.partitioned} is TRUE, otherwise NULL.}
#' \item{data_partitioned [logic]}{TRUE when \code{test_samples} are
#'   provided.}
#' }
#'
#' Object methods:
#' [summary]
#' [split_data]
#'
#' @param input [data.frame] input \code{data.frame} with samples as rows and
#'   features as columns. There must be an extra column with class labels.
#' @param labels_column [numeric] a number indicating the column index of the
#'   class label variable.
#'
#' @return An object of [ml_data] class.
#'
#' @examples
#' set.seed(1234)
#' df <- data.frame(a = runif(10, 0, 1),
#'                  b = runif(10, 0, 1),
#'                  c = as.factor(as.integer(runif(10, 0, 2))))
#' data <- ml_data(df, 3)
#' summary(data)
#'
#' @export
ml_data <- function(input, labels_column = 1) {
  # Preparing input data.
  if (!is.data.frame(input)) {
    stop("input must be a data.frame.")
  }
  if (dim(input)[1] < 2 | dim(input)[2] < 2) {
    stop("input must have 2 or more rows/columns.")
  }
  # Row names are taken as sample names.
  sample_names <- rownames(input)
  # Get feature types.
  feature_types <- as.character(unique(lapply(input[, -labels_column], class)))
  # Preparing labels.
  if (labels_column > dim(input)[2] | labels_column < 1) {
    stop("Wrong labels_column. Out of input range.")
  }
  data_labels <- input[, labels_column]

  # Convert data_labels into a factor and create valid level names.
  data_labels <- as.factor(data_labels)
  levels(data_labels) <- make.names(levels(data_labels))
  if (nlevels(data_labels) < 2) {
    stop("labels_column must have at least 2 levels.")
  }
  classes <- levels(data_labels)
  # Preparing training/test data.
  test_data <- NULL
  test_labels <- NULL
  test_sample_names <- NULL
  training_data <- input[, -labels_column]
  training_labels <- data_labels
  training_sample_names <- sample_names
  data_partitioned <- FALSE
  # Constructor.
  structure(list(classes = classes,
                 feature_types = feature_types,
                 sample_names = sample_names,
                 feature_names = colnames(input[-labels_column]),
                 training_data = training_data,
                 training_labels = training_labels,
                 training_sample_names = training_sample_names,
                 test_data = test_data,
                 test_labels = test_labels,
                 test_sample_names = test_sample_names,
                 data_partitioned = data_partitioned), class = "ml_data")
}

#===============================================================================
# Methods for 'ml_data' class objects.
#===============================================================================
#' Prints a summary for a \code{ml_data} object
#'
#' \code{summary} method for objects of the class \code{ml_data}.
#'
#' @param object [ml_data] an object of the class \code{ml_data}.
#' @param ... additional arguments to pass to \code{summary.default}. See
#'   \link[base]{summary}.
#'
#' @return a printed summary.
#'
#' @examples
#' set.seed(1234)
#' df <- data.frame(a = runif(10, 0, 1),
#'                  b = runif(10, 0, 1),
#'                  c = as.factor(as.integer(runif(10, 0, 2))))
#' data <- ml_data(df, 3)
#' summary(data)
#'
#' @export
summary.ml_data <- function(object, ...) {
  cat("Summary of 'ml_data' class object.\n")
  cat(paste("Data partitioned:", object$data_partitioned, "\n"))
  if (object$data_partitioned) {
    cat(paste("Training samples:", length(object$training_labels), "\n"))
    cat(paste("Test samples:", length(object$test_labels), "\n"))
    cat("Training data classes:\n")
    for (i in 1:length(table(object$training_labels))) {
      cat("\tClass label:", names(table(object$training_labels)[i]),
          "/ Number of cases:", table(object$training_labels)[i], "\n")
    }
    cat("Test data classes:\n")
    for (i in 1:length(table(object$test_labels))) {
      cat("\tClass label:", names(table(object$test_labels)[i]),
          "/ Number of cases:", table(object$test_labels)[i], "\n")
    }
    cat(paste("Number of features:", length(object$feature_names)), "\n")
    cat("Feature types:\n")
    cat("\t", paste(object$feature_types, collapse = " "), "\n")
  } else {
    cat(paste("Training samples:", length(object$training_labels), "\n"))
    cat("Test data not present.\n")
    cat("Training data classes:\n")
    for (i in 1:length(table(object$training_labels))) {
      cat("\tClass label:", names(table(object$training_labels)[i]),
          "/ Number of cases:", table(object$training_labels)[i], "\n")
    }
    cat(paste("Number of features:", length(object$feature_names)), "\n")
    cat("Feature types:\n")
    cat("\t", paste(object$feature_types, collapse = " "), "\n")
  }
}

#===============================================================================
#' Splits data into training and data set.
#'
#' Method only available for objects of the \code{ml_data} class. It splits
#' the actual training data of an unpartitioned \code{ml_data} object into
#' training data and test data, using caret::createDataPartition().
#'
#' @param ml_obj [ml_data] an object of the class \code{ml_data}.
#' @param splitting [numeric] an optional numeric to express the proportion
#'   of cases to be used in training data. Default: 2/3.
#'
#' @return the original ml_obj, but with data partitioned.
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' if (require(GEOquery) &
#'     require(Biobase)) {
#'   gse <- getGEO("GSE69229")
#'   ml <- get_GEO_methylarray(gse = gse, target_name = "outcome:ch1")
#'   summary(ml)
#'   ml_partitioned <- split_data(ml)
#'   summary(ml_partitioned)
#' }
#' }
#'
#' @export
split_data <- function(ml_obj, splitting) {
  cat(paste("Splitting data of", deparse(substitute(ml_obj)), "object.\n"))
  UseMethod("split_data", ml_obj)
}
#' @export
split_data.default <- function(ml_obj, splitting) {
  message("split_data method is only defined for objects of ml_data class.\n")
  message("Nothing was done.\n")
  return(ml_obj)
}
#' @export
split_data.ml_data <- function(ml_obj, splitting = 2/3) {
  cat(paste("Training data:", round(splitting, 4) * 100, "% of the data.\n"))
  cat(paste("Test data:", round(1 - splitting, 4) * 100, "% of the data.\n"))
  # Checking if data was already partitioned. If it was, nothing is done.
  if (ml_obj$data_partitioned) {
    message("The object was already partitioned. No changes were made.\n")
    return(ml_obj)
  }
  all_labels <- ml_obj$training_labels
  # Split samples using caret::createDataPartition().
  partition <- caret::createDataPartition(all_labels, p = splitting)[[1]]
  training_data <- ml_obj$training_data[partition, ]
  training_labels <- all_labels[partition]
  test_data <- ml_obj$training_data[-partition, ]
  test_labels <- all_labels[-partition]
  # Training data should contain all class labels.
  if (nlevels(training_labels) != nlevels(all_labels)) {
    warning("Training data must contain examples of all classes.\n")
    warning("This is usually a consequence of small sample size.\n")
    warning("As splitting is a random process, another call to split_data() may
            solve this problem.\n")
    warning("No changes were made.\n")
    return(ml_obj)
  }
  # Checking partition scheme and warn about data descompensation.
  if (length(test_labels) > length(training_labels)) {
    warning("Training data is smaller than test data. This may be correct but a
            double check of splitting argument is suggested.\n")
  }
  # Modify ml_data object and return it.
  ml_obj$training_data <- training_data
  ml_obj$training_labels <- training_labels
  ml_obj$test_data <- test_data
  ml_obj$test_labels <- test_labels
  ml_obj$data_partitioned <- TRUE

  return(ml_obj)
}
