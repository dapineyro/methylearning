#===============================================================================
# Definition of 'ml_fs' class and its methods. This is an 'ml_data' derived
# class.
# Author: David Piñeyro
# License: GPL-3
# Date: 2018-04-28
#===============================================================================

#' Create \code{ml_fs} objects
#'
#' Constructor for \code{ml_fs} objects. This object takes an
#' \code{ml_data} object, performs one or more feature selection procedures and
#' returns an \code{ml_fs} object, which encapsulates all
#' \code{ml_data} and results of the feature selection procedures.
#'
#' Object members:
#' All the \code{ml_data} members, plus:
#' \describe{
#' \item{selection_methods [character]}{Selections methods applied.}
#' \item{selection_results [list]}{List of character vectors, with the feature
#'   selection of each method applied.}
#' \item{computation_time [list]}{List of numeric vectors, with the computation
#'   time elapsed in each feature selection computation.}
#' }
#'
#' Object methods:
#' all the \code{ml_data} objects plus:
#' \describe{
#' [selection_methods]
#' [selection_results]
#' [plot_venn_fs]
#' [selection_df]
#' [plot_fs_time]
#' }
#'
#' @param ml_data_obj [ml_data] an object of \code{ml_data} class.
#' @param fs_methods [character] a character vector containing one or more of
#'   the implemented methods. Use \code{available_fs_methods} to get all the
#'   methods currently implemented. Alternatively, one can specify "all", to
#'   select all the implemented methods, "filters", to select all the filter
#'   methods and "wrappers", to select all the wrappers.
#' @param selection_size [numeric] the maximum number of features to collect
#'   from each method.
#' @param use_all_data [logical] TRUE makes to use all the data (training +
#'   test) to feature selection. FALSE makes to use only the training data.
#'   Default: TRUE.
#' @param cores [numeric] number of processor cores to use. Default: 1.
#'
#' @return An object of \code{ml_fs} class.
#'
#' @examples
#' \dontrun{
#' if (require(GEOquery) &
#'     require(Biobase)) {
#'   gse <- getGEO("GSE69229")
#'   ml <- get_GEO_methylarray(gse = gse, target_name = "outcome:ch1")
#'   summary(ml)
#'   ml_partitioned <- split_data(ml)
#'   summary(ml_partitioned)
#'   ml_fs_obj <- ml_fs(ml_partitioned, fs_methods = "anova_fs",
#'                      selection_size = 50)
#'   summary(ml_fs_obj)
#' }
#' }
#'
#' @include fs_methods.R
#'
#' @export
ml_fs <- function(ml_data_obj, fs_methods, selection_size, use_all_data = TRUE,
                  cores = 1) {
  # Check if input object is from ml_data class.
  if (!("ml_data" %in% class(ml_data_obj))) {
    stop("ERROR: Input object is not from 'ml_data' class.")
  }
  # Check if data was partitioned and if all data should be used.
  if (ml_data_obj$data_partitioned & use_all_data) {
    # Merging training and test to select features.
    all_data <- rbind(ml_data_obj$training_data, ml_data_obj$test_data)
    all_labels <- factor(c(ml_data_obj$training_labels,
                           ml_data_obj$test_labels),
                         labels  = levels(ml_data_obj$training_labels))
  } else {
    all_data <- ml_data_obj$training_data
    all_labels <- ml_data_obj$training_labels
  }
  # Get the list of available feature selection methods.
  available_fs <- available_fs_methods()
  # Initialize a list of results, in which the name of the item is the method
  # and set it to FALSE, waiting for results.
  selection_results <- split(rep(FALSE, length(available_fs$all)),
                             f = available_fs$all)
  # Perform indicated feature selection methods.
  if (fs_methods[1] == "all") {
    fs_methods <- available_fs$all
  } else if (fs_methods[1] == "filters") {
    fs_methods <- available_fs$filters
  } else if (fs_methods[1] == "wrappers") {
    fs_methods <- available_fs$wrappers
  }
  # Var to store the methods applied.
  selection_methods <- fs_methods
  # Check whether indicated methods are available.
  availability <- vapply(fs_methods, function(x) x %in% available_fs$all,
                         logical(1))
  if (length(fs_methods) > sum(availability)) {
    stop(paste("Error: the following methods are not available:",
               names(availability)[availability == FALSE]))
  }
  # Set a variable to store computation time. It will be a named list.
  computation_time <- list()
  # Run first all the methods with built-in parallelization.
  parallelizable <- c("anova_fs", "mRMR_fs", "ga_fs")
  for (m in selection_methods) {
    if (m %in% parallelizable) {
      cat(paste0("Executing feature selection method: '", m, "'. Using ",
                 cores, " cores.\n"))
      method <- get(m)
      c_time <- system.time(
        sel_m <- method(all_data, all_labels, selection_size, cores)
      )
      # Fill results.
      selection_results[[m]] <- sel_m
      # Fill computation time.
      computation_time[[m]] <- as.numeric(c_time["elapsed"])
      cat(paste0("Execution of '", m, "' completed.\n"))
      # Remove this method from the list.
      fs_methods <- fs_methods[-which(fs_methods == m)]
    }
  }
  # Boruta is also parallelizable, but does not use "cores" argument.
  if ("boruta_fs" %in% fs_methods) {
    cat("Executing feature selection method: 'boruta_fs'. Using all available cores.\n")
    # This function utilize all available cores.
    c_time <- system.time(
      sel_boruta <- boruta_fs(all_data, all_labels, selection_size)
    )
    # Fill results.
    selection_results$boruta_fs <- sel_boruta
    # Fill computation time.
    computation_time[["boruta_fs"]] <- as.numeric(c_time["elapsed"])
    cat("Execution of 'boruta_fs' completed.\n")
    # Remove this method from the list.
    fs_methods <- fs_methods[-which(fs_methods == "boruta_fs")]
  }
  # random_fs also uses different arguments, so doing it sepparately.
  if ("random_fs" %in% fs_methods) {
    cat("Executing feature selection method: 'random_fs'. Using 1 core.\n")
    # This function utilize all available cores.
    c_time <- system.time(
      sel_random <- random_fs(colnames(all_data), selection_size)
    )
    # Fill results.
    selection_results$random_fs <- sel_random
    # Fill computation time.
    computation_time[["random_fs"]] <- as.numeric(c_time["elapsed"])
    cat("Execution of 'random_fs' completed.\n")
    # Remove this method from the list.
    fs_methods <- fs_methods[-which(fs_methods == "random_fs")]
  }
  if (length(fs_methods) == 1) {
    # Only one remaining method, no parallelization.
    cat(paste0("Executing feature selection method: '", fs_methods,
               "'. Single threaded execution.\n"))
    method <- get(fs_methods)
    c_time <- system.time(
      sel_x <- method(all_data, all_labels, selection_size)
    )
    # Fill_results.
    selection_results[[fs_methods]] <- sel_x
    # Fill computation time.
    computation_time[[fs_methods]] <- as.numeric(c_time["elapsed"])
    cat(paste0("Execution of '", fs_methods,"' completed.\n"))
  }
  if (length(fs_methods) > 1) {
    # Only single-threaded methods remain. Single core run.
    cat(paste("Executing the following feature selection methods:",
               paste(fs_methods, collapse = ", "), "using single core.\n"))
    filters_results <- lapply(fs_methods, function(x) {
      method <- get(x)
      c_time <- system.time(
        sel <- method(all_data, all_labels, selection_size)
      )
      return(list(selection = sel, computation = as.numeric(c_time["elapsed"])))
    })
    cat(paste("Execution of", paste(fs_methods, collapse = ", "),
              "completed.\n"))
    # Fill results
    i <- 1
    for (method in fs_methods) {
      selection_results[[method]] <- filters_results[[i]][["selection"]]
      computation_time[[method]] <- filters_results[[i]][["computation"]]
      i <- i + 1
    }
  }
  # Constructor.
  structure(list(classes = ml_data_obj$classes,
                 feature_types = ml_data_obj$feature_types,
                 sample_names = ml_data_obj$sample_names,
                 feature_names = ml_data_obj$feature_names,
                 training_data = ml_data_obj$training_data,
                 training_labels = ml_data_obj$training_labels,
                 training_sample_names = ml_data_obj$training_sample_names,
                 test_data = ml_data_obj$test_data,
                 test_labels = ml_data_obj$test_labels,
                 test_sample_names = ml_data_obj$test_sample_names,
                 data_partitioned = ml_data_obj$data_partitioned,
                 selection_methods = selection_methods,
                 selection_results = selection_results,
                 computation_time = computation_time),
            class = c("ml_fs", "ml_data"))
}

#===============================================================================
# Methods for 'ml_fs' class objects.
#===============================================================================
#' Prints a summary for a \code{ml_fs} object
#'
#' \code{summary} method for objects of the class \code{ml_fs}. It
#' prints all the information inherited from \code{ml_data} class, plus the
#' specific \code{ml_fs} data.
#'
#' @param object [ml_fs] an object of the class
#'   \code{ml_fs}.
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
#' data_sel <- ml_fs(data, fs_methods = "anova_fs", selection_size = 2)
#' summary(data_sel)
#'
#' @export
summary.ml_fs <- function(object, ...) {
  cat("Selection methods applied: \n")
  cat(paste("\t", paste(object$selection_methods, collapse = ", "), "\n"))
  cat("\n")
  cat("Number of features selected: \n")
  for (m in object$selection_methods){
    num <- length(object$selection_results[[m]])
    cat(paste0("\t", m,": ", num, " features selected\n"))
  }
  cat("\n")
  cat("Computation time for each feature selection method: \n")
  for (i in 1:length(object$computation_time)) {
    method_name <- names(object$computation_time)[i]
    method_time <- round(object$computation_time[[i]], 3)
    cat(paste0("\t", method_name, ": ", method_time, " seconds\n"))
  }
  cat("\n")
  NextMethod("summary", object)
}

#===============================================================================
#' Get Feature Selection methods used
#'
#' This function returns the feature selection methods already used to create
#' an object of the \code{ml_fs} class.
#'
#' @param fs_object [ml_fs] an object of \code{ml_fs}
#'   class.
#'
#' @return A \code{character} vector with the used feature selection methods.
#'
#' @examples
#' set.seed(1234)
#' df <- data.frame(a = runif(10, 0, 1),
#'                  b = runif(10, 0, 1),
#'                  c = as.factor(as.integer(runif(10, 0, 2))))
#' data <- ml_data(df, 3)
#' data_sel <- ml_fs(data, fs_methods = "filters", selection_size = 2)
#' m <- selection_methods(data_sel)
#' print(m)
#'
#' @export
selection_methods <- function(fs_object){
  UseMethod("selection_methods", fs_object)
}
#' @export
selection_methods.default <- function(fs_object){
  stop("Error: selection_methods is only defined for objects of ml_fs class.\n")
}
#' @export
selection_methods.ml_fs <- function(fs_object) {
  return(fs_object$selection_methods)
}

#===============================================================================
#' Get Feature Selection Results
#'
#' This function returns the results of previous feature selection, stored in an
#' object of the \code{ml_fs} class.
#'
#' @param fs_object [ml_fs] an object of \code{ml_fs}
#'   class.
#'
#' @return a \code{list} with the \code{character} vector of the selected
#'   features, for each of the used methods.
#'
#' @examples
#' set.seed(1234)
#' df <- data.frame(a = runif(10, 0, 1),
#'                  b = runif(10, 0, 1),
#'                  c = as.factor(as.integer(runif(10, 0, 2))))
#' data <- ml_data(df, 3)
#' data_sel <- ml_fs(data, fs_methods = "filters", selection_size = 2)
#' m <- selection_results(data_sel)
#' print(m)
#'
#' @export
selection_results <- function(fs_object){
  UseMethod("selection_results", fs_object)
}
#' @export
selection_results.default <- function(fs_object) {
  stop("Error: selection_results is only defined for objects of ml_fs class.\n")
}
#' @export
selection_results.ml_fs <- function(fs_object) {
  i <- 1  # To iterate over fs_object$selection_results.
  j <- 1  # To iterate over res.
  res <- list()
  while (i <= length(fs_object$selection_results)) {
    if (fs_object$selection_results[[i]][1] != FALSE) {
      res[[j]] <- fs_object$selection_results[[i]]
      names(res)[j] <- names(fs_object$selection_results)[i]
      j <- j + 1
    }
    i <- i + 1
  }
  return(res)
}

#===============================================================================
#' Calculate shared features
#'
#' This function calculates shared features between the indicated columns.
#' It was meant to be called from \code{plot_venn_fs}.
#'
#' @param df [data.frame] a \code{data.frame} where rows are features and
#'   columns are feature selection methods. Only \code{logical} values for the
#'   cell are allowed.
#' @param x [numeric] the selected column positions of \code{df}.
#'
#' @return the number of rows where columns in \code{x} are TRUE, i.e. the
#'   shared features (rows) between indicated columns (\code{x}, fs_methods).
#'
agreed <- function(df, x) {
  ag <- df
  for (i in 1:length(x)) {
    ag <- subset(ag, ag[, x[i]])
  }
  return(nrow(ag))
}

#' Plot Venn Diagram of Feature Selection Results
#'
#' This function uses \link[VennDiagram]{VennDiagram} package to plot Venn
#' diagrams of the Feature Selection results, for up to 5 methods.
#'
#' @details For 2 methods comparison \link[VennDiagram]{draw.pairwise.venn}
#'   is used. For 3 methods comparison \link[VennDiagram]{draw.triple.venn}
#'   is used. For 4 methods comparison \link[VennDiagram]{draw.quad.venn}
#'   is used. For 5 methods comparison \link[VennDiagram]{draw.quintuple.venn}
#'   is used. Additional arguments for each of the
#'   \link[VennDiagram]{VennDiagram} functions are available.
#'
#' @param fs_object [ml_fs] an object of \code{ml_fs} class.
#' @param ... additional arguments for \link[methylearning]{plot_venn_fs.ml_fs}.
#'
#' @return a Venn diagram plot.
#'
#' @examples
#' set.seed(1234)
#' df <- data.frame(a = runif(10, 0, 1),
#'                  b = runif(10, 0, 1),
#'                  c = as.factor(as.integer(runif(10, 0, 2))))
#' data <- ml_data(df, 3)
#' data_sel <- ml_fs(data, fs_methods = "filters", selection_size = 2)
#' plot_venn_fs(data_sel, fs_methods = c("anova_fs", "mRMR_fs"),
#'              fill = c("red", "blue"), category = c("anova_fs", "mRMR_fs"))
#'
#' @export
plot_venn_fs <- function(fs_object, ...){
  UseMethod("plot_venn_fs")
}
#' @export
plot_venn_fs.default <- function(fs_object, ...) {
  stop("Error: plot_venn_fs is only defined for objects of ml_fs class.\n")
}
#' Plot Venn Diagram of Feature Selection Results
#'
#' This function uses \link[VennDiagram]{VennDiagram} package to plot Venn
#' diagrams of the Feature Selection results, for up to 5 methods.
#'
#' @details For 2 methods comparison \link[VennDiagram]{draw.pairwise.venn}
#'   is used. For 3 methods comparison \link[VennDiagram]{draw.triple.venn}
#'   is used. For 4 methods comparison \link[VennDiagram]{draw.quad.venn}
#'   is used. For 5 methods comparison \link[VennDiagram]{draw.quintuple.venn}
#'   is used. Additional arguments for each of the
#'   \link[VennDiagram]{VennDiagram} functions are available.
#'
#' @param fs_object [ml_fs] an object of \code{ml_fs} class.
#' @param fs_methods [character] a \code{character} vector with the name of
#'   the feature selection methods to compare. From 1 to 5 methods are allowed.
#' @param ... additional optional arguments for
#'   \link[VennDiagram]{draw.pairwise.venn},
#'   \link[VennDiagram]{draw.triple.venn},
#'   \link[VennDiagram]{draw.quad.venn} or
#'   \link[VennDiagram]{draw.quintuple.venn}.
#'
#' @export
plot_venn_fs.ml_fs <- function(fs_object, fs_methods, ...) {
  if (length(fs_methods) > 5 | length(fs_methods) < 1) {
    stop("Venn diagrams can only be plotted from 1 up to 5 methods.")
  }
  grid::grid.newpage()
  df <- selection_df(fs_object)
  df <- df[,-dim(df)[2]]  # Remove agree (last) column.
  # Plot depending on length(fs_methods).
  if (length(fs_methods) == 1) {
    p <- VennDiagram::draw.single.venn(agreed(df, fs_methods), ...)
  }
  if (length(fs_methods) == 2) {
    p <- VennDiagram::draw.pairwise.venn(agreed(df, fs_methods[1]),
                                         agreed(df, fs_methods[2]),
                                         agreed(df, fs_methods[1:2]), ...)
  }
  if (length(fs_methods) == 3) {
    p <- VennDiagram::draw.triple.venn(agreed(df, fs_methods[1]),
                                       agreed(df, fs_methods[2]),
                                       agreed(df, fs_methods[3]),
                                       agreed(df, fs_methods[1:2]),
                                       agreed(df, fs_methods[2:3]),
                                       agreed(df, fs_methods[c(1, 3)]),
                                       agreed(df, fs_methods), ...)
  }
  if (length(fs_methods) == 4) {
    p <- VennDiagram::draw.quad.venn(agreed(df, fs_methods[1]),
                                     agreed(df, fs_methods[2]),
                                     agreed(df, fs_methods[3]),
                                     agreed(df, fs_methods[4]),
                                     agreed(df, fs_methods[1:2]),
                                     agreed(df, fs_methods[c(1, 3)]),
                                     agreed(df, fs_methods[c(1, 4)]),
                                     agreed(df, fs_methods[2:3]),
                                     agreed(df, fs_methods[c(2, 4)]),
                                     agreed(df, fs_methods[3:4]),
                                     agreed(df, fs_methods[1:3]),
                                     agreed(df, fs_methods[c(1, 2, 4)]),
                                     agreed(df, fs_methods[c(1, 3, 4)]),
                                     agreed(df, fs_methods[2:4]),
                                     agreed(df, fs_methods), ...)
  }
  if (length(fs_methods) == 5) {
    p <- VennDiagram::draw.quintuple.venn(agreed(df, fs_methods[1]),
                                          agreed(df, fs_methods[2]),
                                          agreed(df, fs_methods[3]),
                                          agreed(df, fs_methods[4]),
                                          agreed(df, fs_methods[5]),
                                          agreed(df, fs_methods[1:2]),
                                          agreed(df, fs_methods[c(1, 3)]),
                                          agreed(df, fs_methods[c(1, 4)]),
                                          agreed(df, fs_methods[c(1, 5)]),
                                          agreed(df, fs_methods[2:3]),
                                          agreed(df, fs_methods[c(2, 4)]),
                                          agreed(df, fs_methods[c(2, 5)]),
                                          agreed(df, fs_methods[3:4]),
                                          agreed(df, fs_methods[c(3, 5)]),
                                          agreed(df, fs_methods[4:5]),
                                          agreed(df, fs_methods[1:3]),
                                          agreed(df, fs_methods[c(1, 2, 4)]),
                                          agreed(df, fs_methods[c(1, 2, 5)]),
                                          agreed(df, fs_methods[c(1, 3, 4)]),
                                          agreed(df, fs_methods[c(1, 3, 5)]),
                                          agreed(df, fs_methods[c(1, 4, 5)]),
                                          agreed(df, fs_methods[2:4]),
                                          agreed(df, fs_methods[c(2, 3, 5)]),
                                          agreed(df, fs_methods[c(2, 4, 5)]),
                                          agreed(df, fs_methods[3:5]),
                                          agreed(df, fs_methods[1:4]),
                                          agreed(df, fs_methods[c(1, 2, 3, 5)]),
                                          agreed(df, fs_methods[c(1, 2, 4, 5)]),
                                          agreed(df, fs_methods[c(1, 3, 4, 5)]),
                                          agreed(df, fs_methods[2:5]),
                                          agreed(df, fs_methods), ...)
  }
  if (!exists("p")) {
    stop("Error: venn diagram not drawn.")
  }
  return(p)
}

#===============================================================================
#' Create a \code{data.frame} with feature selection results
#'
#' This function creates a \code{data.frame} with the feature selection results.
#'
#' @param fs_object [ml_fs] an object of \code{ml_fs} class.
#'
#' @return A \code{data.frame} where unique features selected are rows and
#'   feature selection methods used are columns. TRUE/FALSE indicated whether a
#'   mathod selected a feature. The last column is called \code{agree} and
#'   counts the number of methods that agreed in each feature selection.
#'   The returned \code{data.frame} is deacreasing ordered by \code{agree}
#'   column, so most selected features are first.
#'
#' @examples
#' set.seed(1234)
#' df <- data.frame(a = runif(10, 0, 1),
#'                  b = runif(10, 0, 1),
#'                  c = as.factor(as.integer(runif(10, 0, 2))))
#' data <- ml_data(df, 3)
#' data_sel <- ml_fs(data, fs_methods = "filters", selection_size = 2)
#' df <- selection_df(data_sel)
#' print(df)
#' # Select the most selected feature.
#' rownames(df)[1]
#'
#' @export
selection_df <- function(fs_object){
  UseMethod("selection_df", fs_object)
}
#' @export
selection_df.default <- function(fs_object) {
  stop("Error: selection_df is only defined for objects of ml_fs class.\n")
}
#' @export
selection_df.ml_fs <- function(fs_object) {
  sel <- selection_results(fs_object)
  # Generating a data.frame with whether a method selected a particular feature.
  methods <- names(sel)  # Columns.
  features <- unique(unlist(sel))  # Rows.
  df <- as.data.frame(lapply(methods, function(x) features %in% sel[[x]]))
  colnames(df) <- methods
  rownames(df) <- features
  # Add agree column.
  df$agree <- apply(df, 1, sum)
  # Decreasing ordering by agree column.
  df <- df[order(-df$agree), ]
  return(df)
}

#===============================================================================
#' Plot Computation Time for \code{ml_fs} object
#'
#' This function creates a plot of the computation time elapsed in each of the
#' feature selection calculations performed on a \code{ml_fs} object.
#'
#' @param fs_object [ml_fs] an object of \code{ml_fs} class.
#'
#' @return a barchart with the time elapsed in each feature selection
#'   computation.
#'
#' @examples
#' set.seed(1234)
#' df <- data.frame(a = runif(10, 0, 1),
#'                  b = runif(10, 0, 1),
#'                  c = as.factor(as.integer(runif(10, 0, 2))))
#' data <- ml_data(df, 3)
#' data_sel <- ml_fs(data, fs_methods = "filters", selection_size = 2)
#' plot_fs_time(data_sel)
#'
#' @export
plot_fs_time <- function(fs_object){
  UseMethod("plot_fs_time", fs_object)
}
#' @export
plot_fs_time.default <- function(fs_object) {
  stop("Error: plot_fs_time is only defined for objects of ml_fs class.\n")
}
#' @export
plot_fs_time.ml_fs <- function(fs_object) {
  fs_methods <- names(fs_object$computation_time)
  fs_times <- round(as.numeric(unlist(fs_object$computation_time)), 3)
  df2plot <- data.frame(fs_methods, fs_times)
  # Plot
  ggplot2::ggplot(df2plot,
                  ggplot2::aes(x=fs_methods, y=fs_times, fill=fs_methods)) +
    ggplot2::geom_bar(position=ggplot2::position_dodge(), stat="identity",
                      colour="black", # Use black outlines,
                      size=.3) +      # Thinner lines
    # Add values on top.
    ggplot2::geom_text(
      data = df2plot,
      ggplot2::aes(x = fs_methods, y = fs_times, label = fs_times),
      vjust = 0) +
    ggplot2::xlab("Feature selection method") +
    ggplot2::ylab("Elapsed time (seconds)") +
    ggplot2::scale_fill_hue(name="") +
    ggplot2::ggtitle("Computation time for each Feature Selection method") +
    ggplot2::theme_bw()
}
