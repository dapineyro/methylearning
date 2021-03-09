#===============================================================================
# Implementation of feature selection methods.
# Author: David Piñeyro
# License: GPL-3
# Date: 2018-04-30
#===============================================================================

## Filter methods ==============================================================

#' Random Feature Selection
#'
#' Implementation of random feature selection, for testing and performance
#' evaluation.
#'
#' @param features [character] a \code{character} vector with the data
#'   \code{colnames}, i.e. feature names (probes or CpGs).
#' @param n [numeric] the number of features to select.
#'
#' @return [character] A character vector with the features \code{colnames}
#'   selected.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' features <- colnames(input_data)
#' features_selected <- random_fs(features, 2)
#' print(features_selected)
#'
#' @export
random_fs <- function(features, n) {
  # Check whether the number of features to select is smaller than the total
  # number of features.
  if (length(features) <= n) {
    warning("The number of features to select is not smaller than the total
            number of features. No features were selected.")
    return(features)
  } else {
    selected <- sample(features, n)  # No replacement.
    return(selected)
  }
}

#' ANOVA feature selection
#'
#' Selects features based on ANOVA F-statistic. For each feature, a one-way
#' ANOVA model is fitted. Resulting p-values are ascending sorted and the
#' first \code{n} features are selected, based on this sorting.
#'
#' @param input_data [data.frame] a \code{data.frame} with samples as rows and
#'   features as columns. Each column is treated as a response (dependent
#'   variable) for an ANOVA model.
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels. It will be used as a predictor (independent
#'   variable) for each ANOVA model. NOTE: this variable should not be included
#'   as a column of \code{input_data}.
#' @param n [numeric] the number of features to select.
#' @param cores [numeric] number of cores to use for parallel ANOVA fitting.
#'   Default: 1.
#'
#' @return [character] A character vector with the features \code{colnames}
#'   selected.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- anova_fs(input_data, sample_labels, 2)
#' print(features_selected)
#'
#' @export
anova_fs <- function(input_data, sample_labels, n, cores = 1) {
  # First, check whether n is < total number of features.
  if (dim(input_data)[2] <= n) {
    warning("The number of features to select is not smaller than the total
            number of features. No features were selected.")
    return(colnames(input_data))
  } else {
    # Fits an anova to each probe as dependent (predicted) variable and
    # sample_labels as independent (predictor) variable. Save all results in a
    # list.
    pvals <- unlist(parallel::mclapply(input_data, function(x) {
      df <- data.frame(x, sample_labels)
      fit <- stats::aov(x ~ sample_labels, data = df)
      pval <- summary(fit)[[1]]$`Pr(>F)`[1]
      return(pval)
    }, mc.cores = cores))
    # Sort ascending.
    pvals_sorted <- sort(pvals)
    # Select first n features (most signifiant).
    return(names(pvals_sorted[1:n]))
  }
}

#' Select CpGs by differential methylation using limma
#'
#' This funtion uses limma to select top differentially methylated probes
#' (features) between groups (sample labels factor). The number of class labels
#' can be 2 or more, but no paired designs are implemented. It is equivalent to
#' anova_fs, but much faster. Only applicable to methylation array data.
#'
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (beta values or methylation ratios). The
#'   default behaviour is to assume a betas data.frame or matrix, but an
#'   already converted to m-values matrix is allowed, setting `are_m_vals` to
#'   TRUE.
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#' @param n [numeric] the number of features (top differentially methylated CpGs
#'   to select).
#' @param remove_dup [bool] Whether to return a unique set of features or not.
#'   This is important when there are more than 2 classes and any feature is
#'   repeatedly selected for more than one class. Default: TRUE.
#' @param are_m_vals [bool] Set to TRUE when an m-values matrix or data.frame
#'   is passed as `input_data`. Default: FALSE.
#'
#' @return [character] A character vector with the \code{n} top differentially
#'   methylated (sorted by p-value) feature names, for all contrasts tested.
#'   For more than two groups in \code{sample_labels}, each contrast is designed
#'   as one Vs the others.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- limma_fs(input_data, sample_labels, 2)
#' print(features_selected)
#'
#' @export
limma_fs <- function(input_data, sample_labels, n, remove_dupl = TRUE,
                     are_m_vals = FALSE) {
  # First, check whether n is < total number of features.
  if (dim(input_data)[2] <= n) {
    warning("The number of features to select is not smaller than the total
            number of features. No features were selected.")
    return(colnames(input_data))
  } else if (sum(input_data == 0) | sum(input_data == 1)) {
    # Check whether absolute 0 or 1 values exists. So it is usually a sign of
    # sequencing data.
    warning("0 and/or 1 are found in your data. Is this data generated by
            bisulfite sequencing?. Limma filtering is not applied.")
    return(colnames(input_data))
  } else {
    # Calculate M-values using M=log2(Beta/(1-Beta)). All statistics will be
    # performed on M-values. Transpose to have rows as features and columns as
    # samples.
    sample_names <- rownames(input_data)
    input_data <- t(input_data)
    if (!are_m_vals) {
      m_vals <- log2(input_data / (1-input_data))
    } else {
      m_vals <- input_data
    }
    # Create a targets dataframe.
    pheno_data <- data.frame(sample_names, sample_labels)
    rownames(pheno_data) <- sample_names
    targets <- stats::model.frame(sample_names ~ sample_labels, pheno_data)
    # Design matrix (only unpaired test supported).
    design <- stats::model.matrix(~0+sample_labels, data=targets)
    colnames(design) <- levels(sample_labels)
    # Contrast matrix (one vs the others).
    if (nlevels(sample_labels) == 2) {
      contr <- paste0(levels(sample_labels)[1], "-", levels(sample_labels)[2])
      contMatrix <- limma::makeContrasts(contrasts = contr, levels = design)
    } else {
      # More than 2 groups. Using One Vs the Others contrasts.
      i <- 1
      contr <- character()
      while (i <= nlevels(sample_labels)) {
        one <- levels(sample_labels)[i]
        the_others <- levels(sample_labels)[-i]
        # We will make contrasts such as "A-(B+C)/2".
        the_others_mean <- paste0("(", paste(the_others, collapse = "+"),
                                  ")/", length(the_others))
        contr <- c(contr, paste0(one, "-", the_others_mean))
        i <- i + 1
      }
      contMatrix <- limma::makeContrasts(contrasts = contr, levels = design)
    }
    # fit the linear model
    fit <- limma::lmFit(m_vals, design)
    # fit the contrasts
    fit2 <- limma::contrasts.fit(fit, contMatrix)
    fit2 <- limma::eBayes(fit2)
    # Toptable of all contrasts, selected by p-value.
    if (nlevels(sample_labels) < 3) {
      list_of_results <- lapply(1:(nlevels(sample_labels) - 1), function(x) {
        limma::topTable(fit2, coef = x, number=Inf, sort.by = "P")
      })
    } else {
      list_of_results <- lapply(1:(nlevels(sample_labels)), function(x) {
        limma::topTable(fit2, coef = x, number=Inf, sort.by = "P")
      })
    }
    # Selecting the "best" CpGs of each contrast. If
    # n/nlevels(sample_labels) has no integer result, round
    # approximation is taken.
    if (nlevels(sample_labels) < 3) {
      each_contrast_n <- round(n/(nlevels(sample_labels) - 1), 0)
    } else {
      each_contrast_n <- round(n/(nlevels(sample_labels)), 0)
    }
    feature_selection <- character()
    for (r in list_of_results) {
      feature_selection <- c(feature_selection, rownames(r)[1:each_contrast_n])
    }
    if (remove_dupl){
      return(unique(feature_selection))
    } else {
      return(feature_selection)
    }
  }
}

#' Correlation-based Feature Selection
#'
#' This funtion uses FSelector::cfs implementation to select features based
#' on correlation with \code{sample_labels} and inter-feature correlation.
#'
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (with beta values or methylation ratios).
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#' @param n [numeric] the max number of features to output.
#'
#' @return [character] A character vector with at most \code{n} selected
#'   features.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- correlation_based_fs(input_data, sample_labels, 2)
#' print(features_selected)
#'
#'
#' @export
correlation_based_fs <- function(input_data, sample_labels, n) {
  # First, check whether n is < total number of features.
  if (dim(input_data)[2] <= n) {
    warning("The number of features to select is not smaller than the total
            number of features. No features were selected.")
    return(colnames(input_data))
  } else {
    # Use FSelector::cfs algorithm.
    selection <- FSelector::cfs(sample_labels ~ .,
                                data = cbind(input_data, sample_labels))
    if (length(selection) > n) {
      # Select the best n features.
      return(selection[1:n])
    } else {
      # Return all the selection.
      return(selection)
    }
  }
}

#' Entropy-based information gain feature selection
#'
#' This funtion uses FSelector::information.gain implementation to select
#' features based on entropy information gain between feature and
#' \code{sample_labels}.
#'
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (with beta values or methylation ratios).
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#' @param n [numeric] the max number of features to output.
#'
#' @return [character] A character vector with at most \code{n} selected
#'   features.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- information_gain_fs(input_data, sample_labels, 2)
#' print(features_selected)
#'
#' @export
information_gain_fs <- function(input_data, sample_labels, n) {
  # First, check whether n is < total number of features.
  if (dim(input_data)[2] <= n) {
    warning("The number of features to select is not smaller than the total
            number of features. No features were selected.")
    return(colnames(input_data))
  } else {
    # Using FSelector::information.gain() algorithm.
    selection <- FSelector::information.gain(sample_labels ~ .,
                                  data = cbind(input_data, sample_labels))
    # Attach rownames (feature names) as a column.
    feature <- rownames(selection)
    selection <- cbind(feature, selection)
    # Select features with ig > 0.
    selection <- selection[selection$attr_importance > 0, ]
    # Order by decreasing imporance.
    selection <- selection[order(selection$attr_importance,
                                 decreasing = TRUE), ]
    if (dim(selection)[1] > n) {
      # Select the best n.
      return(as.character(selection$feature[1:n]))
    } else {
      # Return all as there are less than n.
      return(as.character(selection$feature))
    }
  }
}

#' RELIEF algorithm for feature selection
#'
#' This function implements FSelector::relief algorithm for feature selection.
#' The algorithm finds weights of continous and discrete attributes basing on a
#' distance between instances.
#'
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (with beta values or methylation ratios).
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#' @param n [numeric] the max number of features to output.
#' @param neighbours_count [numeric] number of neighbours to find for every
#'   sampled instance. Default: 5.
#' @param sample_size [numeric] number of instances to sample. Default: 10.
#'
#' @return [character] A character vector with at most \code{n} selected
#'   features.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- relief_fs(input_data, sample_labels, 2)
#' print(features_selected)
#'
#' @export
relief_fs <- function(input_data, sample_labels, n, neighbours_count = 5,
                      sample_size = 10) {
  # First, check whether n is < total number of features.
  if (dim(input_data)[2] <= n) {
    warning("The number of features to select is not smaller than the total
            number of features. No features were selected.")
    return(colnames(input_data))
  } else {
    # Using FSelector::relief algorithm.
    selection <- FSelector::relief(sample_labels ~ .,
                                   data = cbind(input_data, sample_labels),
                                   neighbours.count = neighbours_count,
                                   sample.size = sample_size)
    # Select the best n features.
    selection_best <- FSelector::cutoff.k(selection, n)
    return(selection_best)
  }
}

#' Minimum redundancy, maximum relevance feature filtering
#'
#' This function implements mRMRe::mRMR.classic algorithm for feature selection.
#' The algorithm finds weights of continous and discrete attributes basing on a
#' distance between instances.
#'
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (with beta values or methylation ratios).
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#' @param n [numeric] the max number of features to output.
#' @param cores [numeric] number of cores to use.
#'
#' @return [character] A character vector with at most \code{n} selected
#'   features.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- relief_fs(input_data, sample_labels, 2)
#' print(features_selected)
#'
#' @export
mRMR_fs <- function(input_data, sample_labels, n, cores) {
  # First, check whether n is < total number of features.
  if (dim(input_data)[2] <= n) {
    warning("The number of features to select is not smaller than the total
            number of features. No features were selected.")
    return(colnames(input_data))
  } else {
    # Use classic mRMR algorithm.
    # Set the number of cores to use.
    mRMRe::set.thread.count(cores)
    # Prepare mRMRe.Data object.
    mRMR_data <- mRMRe::mRMR.data(data = cbind(as.numeric(sample_labels),
                                               input_data))
    # Create an mRMRe.Filter object, with the results.
    mRMR_res <- mRMRe::mRMR.classic(data = mRMR_data, target_indices = 1,
                                    feature_count = n)
    # Convert to feature names.
    selection <- colnames(input_data)[mRMRe::solutions(mRMR_res)[[1]][,1]]
    return(selection)
  }
}

## Wrappers ====================================================================

#' Calculate feature subset fitness
#'
#' This function calculate feature subset fitness, based on classification
#' accuracy. Classification is performed using the three implemented algorithms:
#'  "knn", "C5.0" and "lda", averaging their results. Each of the algorithms is
#' performed with its default parameters.
#'
#' @param subset [character] a character vector of the feature names of the
#'   particular subset to test accuracy.
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (with beta values or methylation ratios).
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#'
#' @return the mean accuracy achieved by the three algorithms.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' subset <- sample(colnames(input_data), 2)
#' accuracy_subset <- cls_fitness(subset, input_data, sample_labels)
#' print(accuracy_subset)
#'
#' @export
cls_fitness <- function(subset, input_data, sample_labels) {
  ## KNN
  # Use k = sqrt(n), where n is the number of samples.
  k <- as.integer(sqrt(length(sample_labels)))
  knn_predictions <- class::knn(train = input_data[subset],
                                test = input_data[subset],
                                cl = sample_labels, k = k)
  knn_accuracy <- sum(knn_predictions == sample_labels) / length(sample_labels)

  ## C5.0
  # Training decision tree with by default parameters.
  C50_model <- C50::C5.0(x = input_data[subset], y = sample_labels)
  # Predict.
  C50_predictions <- stats::predict(C50_model, input_data[subset])
  C50_accuracy <- sum(C50_predictions == sample_labels) / length(sample_labels)

  ## Linear Discriminant Analysis
  # Create the model.
  lda_model <- MASS::lda(x = input_data[subset], grouping = sample_labels)
  # Make predictions.
  lda_p <- stats::predict(lda_model, input_data[subset])
  lda_accuracy <- sum(lda_p$class == sample_labels) / length(sample_labels)

  return(mean(c(knn_accuracy, C50_accuracy, lda_accuracy)))
}

#' Generate random mutations for Genetic Algorithm
#'
#' This function is intented to be called from ga_fs. Produce mutations, i.e.
#' feature swapping between remaining features and subsets, of each generation.
#'
#' @param subsets [list] a \code{list} of \code{character} vectors with the
#'   subsets selected.
#' @param remaining [character] a \code{character} vector with the features not
#'   selected.
#' @param m_rate [numeric] a number between 0 and 0.2. It will be taken as the
#'   mutation rate. Default: 0.01.
#'
#' @return a \code{list} of two elements: the mutated \code{subsets list} and
#'   the mutated \code{remaining character} vector.
#'
#' @export
ga_mutation <- function(subsets, remaining, m_rate = 0.01) {
  # Check whether m_rate is 0 or too high.
  if (m_rate == 0) {
    return(subsets, remaining)
  } else if (m_rate > 0.2) {
    # A hard limit of 20% mutations is imposed.
    warning("Mutation rate was set too high, use a value between 0 and 0.2.
            No mutation was produced.")
    return(subsets, remaining)
  } else if (m_rate < 0) {
    warning("Mutation rate incorrect, use a value between 0 and 0.2.
            No mutation was produced.")
    return(subsets, remaining)
  } else {
    # Number of features in subsets.
    n_cases_subsets <- length(subsets) * length(subsets[[1]])
    # Number of features to select from both remaining and subsets. This number
    # of features will be the ones that will be swapped between remaining and
    # subsets.
    features_to_swap <- ceiling(n_cases_subsets * m_rate)
    # Randomly, select features to swap from remaining.
    remaining_swap <- sample(1:length(remaining), features_to_swap)
    for (remaining_idx in remaining_swap) {
      # Randomly select the subset.
      subset_idx <- sample(1:length(subsets), 1)
      # Randomly select feature from the selected subset.
      feature_idx <- sample(1:length(subsets[[1]]), 1)
      # Make the swap (the actual mutation).
      remaining_feature <- remaining[remaining_idx]
      subset_feature <- subsets[[subset_idx]][feature_idx]
      subsets[[subset_idx]][feature_idx] <- remaining_feature
      remaining[remaining_idx] <- subset_feature
    }
    return(list(subsets, remaining))
  }
}

#' Wrapped Genetic Algorithm - knn/C5.0/lda for feature selection
#'
#' This function implements a custom genetic algorithm wrapped to several
#' classification methods, to perform feature selection.
#'
#' @details This is a very computational expensive feature selection method that
#' uses a genetic algorithm (GA) to select features. As fitness function to
#' optimize, this GA uses the mean accuracy obtained by the 3 implemented
#' classification algorithms (knn, C5.0 and lda). A recombination and mutation
#' step is performed in each generation.
#'
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (with beta values or methylation ratios).
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#' @param n [numeric] the max number of features to output. It is also the
#'   number of features in each GA "chromosome".
#' @param cores [numeric] number of cores to use in the parallel parts of the
#'   algorithm.
#' @param generations [numeric] number of generations to perform. Default: 500.
#' @param m_rate [numeric] a number between 0 and 0.2. It will be taken as the
#'   mutation rate. Default: 0.01.
#' @param verbose [logical] a logical indicating whether to status reports.
#'   Default: FALSE.
#' @param return_history [logical] if \code{TRUE}, a list of two elements
#'   (\code{$selection} and \code{$max_ac}) is returned. If \code{FALSE}, only
#'   a vector with the selected features is returned. Default: FALSE.
#'
#' @return [character] A character vector with \code{n} selected features. If
#'   \code{return_history = TRUE}, a list is returned instead, with 3 elements:
#'   \code{$selection} [character] (the selected n features); \code{$max_ac}
#'   [numerical] a numerical vector with the max accuracy achieaved in each
#'   generation; \code{$mean_ac} [numerical] a numerical vector with the mean
#'   accuracy of th best 10% subsets selected in each generation.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- ga_fs(input_data, sample_labels, 50, 2)
#' print(features_selected)
#'
#' @export
ga_fs <- function(input_data, sample_labels, n, cores, generations = 500,
                  m_rate = 0.01, verbose = FALSE, return_history = FALSE) {
  # First, check whether n is < total number of features.
  if (dim(input_data)[2] <= n) {
    warning("The number of features to select is not smaller than the total
            number of features. No features were selected.")
    return(colnames(input_data))
  } else {
    # Start time
    start_time <- Sys.time()
    # Initial message.
    if (verbose) {
      cat("Genetic Algorithm coupled to knn/C5.0/lda classifiers\n")
    }
    # Random generation of the initial population.
    # Initialization of the first population.
    all_features <- colnames(input_data)
    # Choromosmes (subsets), which are the individuals of the population, are
    # composed of genes (features). The first initialization is random, and
    # create substets from 3/4 of the total features.
    n_subsets <- ceiling((length(all_features) * 0.75) / n)  # always an int.
    if ((n_subsets * n) >= length(all_features)) {
      warning("The ratio features to select / total features is too high to
               use this Genetic Algorithm. No features were selected.")
      return(colnames(input_data))
    }
    # Random generation of subsets.
    random_idx <- sample(unlist(lapply(0:n_subsets, function(x) {
      if(x == 0) {
        # 0 marks the features that are not going to be used in this iteration.
        remaining <- length(all_features) - n_subsets * n
        rep(x, remaining)
      } else {
        rep(x, n)
      }
    })))
    subsets <- split(all_features, f = as.factor(random_idx))
    # separate remainder features.
    remaining_features <- subsets[[1]]
    subsets <- subsets[2:length(subsets)]
    # Var to store the max accuracy of each generation.
    max_ac <- numeric()
    # Var to store the mean accuracy of the best 10% subsets.
    mean_ac_best_10_all <- numeric()
    # Start evolution.
    for (i in 1:generations) {
      # Calculate fitness of each subset.
      ac <- parallel::mclapply(subsets[1:length(subsets)],
                               cls_fitness, input_data, sample_labels,
                               mc.cores = cores)
      # Storing the max_ac of this generation.
      max_ac <- c(max_ac, max(unlist(ac)))
      # To accelerate evolution, each generation the best 10% are always
      # selected.
      best_10_idx <- sort(unlist(ac),
                          decreasing = TRUE,
                          index.return = TRUE)$ix[1:(ceiling(n_subsets/10))]
      # Calculate the mean accuracy of these best 10 performant subsets.
      mean_ac_best_10 <- mean(sort(unlist(ac),
                              decreasing = TRUE)[1:(ceiling(n_subsets/10))])
      mean_ac_best_10_all <- c(mean_ac_best_10_all, mean_ac_best_10)
      # The last generation, select the best subset.
      if (i == generations) {
        # Order subsets by accuracy.
        best_idx <- best_10_idx[1]
        # Print final status.
        if (verbose) {
          cat("Final generation\n")
          cat("\tSubset of features selected:\n")
          cat(paste(subsets[[best_idx]], "\n"))
          cat(paste("Mean accuracy of the last best 10% feature subsets:",
                    mean_ac_best_10, "\n"))
          cat(paste("Accuracy achieved for the selected subset:",
                    ac[[best_idx]], "\n"))
          # Run time reporting.
          end_time <- Sys.time()
          cat("Elapsed time:\n")
          print(end_time - start_time)
        }
        # Returning depending on return_history value.
        if (return_history) {
          return(list(selection = subsets[[best_idx]], max_ac = max_ac,
                      mean_ac = mean_ac_best_10_all))
        } else {
          return(subsets[[best_idx]])
        }
      }
      # Select for reproduction the best 10% subsets.
      selected_subsets <- list()
      unselected <- character()
      ac_counter <- 1
      selected_counter <- 1
      for (s in subsets) {
        # Making the best 10% performant subsets selected.
        if (ac_counter %in% best_10_idx) {
          fitness <- TRUE
        } else {
          fitness <- FALSE
        }
        if (fitness) {
          # Selected.
          selected_subsets[[selected_counter]] <- s
          selected_counter <- selected_counter + 1
        } else {
          unselected <- c(unselected, s)
        }
        ac_counter <- ac_counter + 1
      }
      # Control to avoid to select all of the subsets (which would stop
      # evolution).
      if (length(selected_subsets) == n_subsets) {
        to_drop <- sample(1:n_subsets, 1)
        unselected <- c(unselected, selected_subsets[[to_drop]])
        selected_subsets <- selected_subsets[-to_drop]
      }
      # Merge unselected and remaining_features, for the next generation
      # newcommers.
      remaining_features <- c(remaining_features, unselected)
      # Generate offspring.
      # Generate couples.
      # Check whether any subset was selected.
      if (length(selected_subsets) == 0){
        couples <- numeric()
      } else {
        couples <- sample(1:length(selected_subsets))  # Each pair is a couple.
      }
      # If there is an unpaired subset, the first one do not recombine.
      if (length(couples) %% 2) {
        couples <- couples[-1]
      }
      j <- 1
      while (j < length(couples)) {
        a_idx <- couples[j]
        b_idx <- couples[j+1]
        # Recombination.
        recomb_idx <- sort(sample(1:n, 2))
        subset_a <- selected_subsets[[a_idx]]
        subset_b <- selected_subsets[[b_idx]]
        selected_subsets[[a_idx]][recomb_idx[1]:recomb_idx[2]] <- subset_b[
          recomb_idx[1]:recomb_idx[2]]
        selected_subsets[[b_idx]][recomb_idx[1]:recomb_idx[2]] <- subset_a[
          recomb_idx[1]:recomb_idx[2]]
        j <- j + 2
      }
      # Generating new subsets from remaining features.
      n_to_generate <- n_subsets - length(selected_subsets)
      # Random generation of subsets.
      random_idx <- sample(unlist(lapply(0:n_to_generate, function(x) {
        if (x == 0) {
          # 0 marks the features that are not going to be used in this iteration
          remaining <- length(all_features) - n_subsets * n
          rep(x, remaining)
        } else {
          rep(x, n)
        }
      })))
      subsets_new <- split(remaining_features, f = as.factor(random_idx))
      # separate remainder features.
      remaining_features <- subsets_new[[1]]
      # Append the recombined offspring to the new ones.
      subsets <- c(selected_subsets, subsets_new[2:length(subsets_new)])
      names(subsets) <- 1:length(subsets)
      ## Perform mutation.
      mutation_results <- ga_mutation(subsets, remaining_features, m_rate)
      subsets <- mutation_results[[1]]
      remaining_features <- mutation_results[[2]]
      # Print a little summary.
      if (verbose) {
        cat(paste("Generation:", i, "\n"))
        cat(paste("\tPopulation size:", n_subsets, "\n"))
        cat(paste("\tSelected:", length(selected_subsets), "\n"))
        cat(paste("\tMean accuracy of the best selected subsets:",
                  mean_ac_best_10, "\n"))
        cat(paste("\tMaximum accuracy achieved:", max(unlist(ac)), "\n"))
      }
    }
  }
}

#' Modified Sequential Feature Selection
#'
#' A modified version of Sequential Feature Selection algorithm.
#'
#' @details This funtion implements a custom version of the Sequential Feature
#'   Selection algorithm. This wrapper is coupled to three classification
#'   algorithms (knn, C5.0 and lda) using \link[methylearning]{cls_fitness}
#'   function. An initial step is performing, selecting only those features
#'   obtaining an accuracy over random expect (0.5 in 2 class problems, 0.33 in
#'   3 class, etc.). On this initial selection, each feature is coupled with
#'   each other and the best performant couples are selected. The rest of the
#'   iterations adds one of the other initial features to each couple, with the
#'   only constrain that no repetition is allowed in each group. Best performant
#'   groups are selected in each iteration. Iteration stops when no improvement
#'   is achieaved in any of the additions or when the groups are of size
#'   \code{n}, which is first meet. If more than one group has the same accuracy
#'   by the \link[methylearning]{cls_fitness} function, only one is chosen at
#'   random.
#'
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (with beta values or methylation ratios).
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#' @param n [numeric] the max number of features to output and, at the same time
#'   the max number of iterations of the algorithm.
#' @param verbose [logical] a logical indicating whether to status reports.
#'   Default: FALSE.
#'
#' @return [character] A character vector with at most \code{n} selected
#'   features.
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- mseq_fs(input_data, sample_labels, 50, 2)
#' print(features_selected)
#'
#' @export
mseq_fs <- function(input_data, sample_labels, n, verbose = FALSE) {
  features <- colnames(input_data)
  # First feature selection based on one feature models. Parallelization
  # disabled due to doTryCatch warnings.
  results <- lapply(features, cls_fitness, input_data, sample_labels)
  # Select only those features performing best than random.
  random_accuracy <- 1/nlevels(sample_labels)
  initial_selection <- vapply(results, function(x) {
    ifelse(x > random_accuracy, TRUE, FALSE)
  }, logical(1))
  features <- features[initial_selection]
  # Generate a list with all the candidates.
  f_list <- split(features, f = as.factor(1:length(features)))
  ### Loop for sequential feature selection.
  # In each iteration, an stop criteria will be no improvement over past
  # iteration accuracy.
  prev_max_ac <- max(unlist(results))
  # If stop condition is meet, the best group (at random if there are more than
  # just one) will be returned by the algorithm.
  max_groups <-f_list[which(unlist(results) == prev_max_ac)]
  if (length(max_groups) > 1) {
    # Randomly select one best, but just in case the next iteration will be the
    # last and a group to return is needed.
    r_num <- sample(1:length(max_groups), 1)
    prev_max_group <- max_groups[[r_num]]
  } else {
    prev_max_group <- max_groups[[1]]  # the only element.
  }
  # Print verbose message.
  if (verbose) {
    cat("Initial selection:\n")
    cat(paste("\tFeatures selected:", length(features), "\n"))
    cat(paste("\tMax accuracy:", prev_max_ac, "\n"))
  }
  # This loop has a max of n iterations, as in each iteration a feature is
  # incorporated.
  for (iter in 1:n) {
    # Initializing variables.
    new_f_list <- list()
    new_ac_list <- list()
    i <- 1
    for (f in f_list) {
      # For each feature group:
      # Remove the group members from initial feature selection.
      features_no_f <- features[which(features != f)]
      # Calculate accuracy adding one feature. Parallelization disabled
      # due doTryCatch warnings.
      res <- lapply(features_no_f, function(x) {
        new_group <- c(f, x)
        ac <- cls_fitness(new_group, input_data, sample_labels)
        return(list(new_group = new_group, ac = ac))
      })
      # New feature groups and their accuracy.
      f_groups <- lapply(res, function(x) x$new_group)
      ac_groups <- lapply(res, function(x) x$ac)
      # Add only the most performant.
      groups_to_add <- which(unlist(ac_groups) == max(unlist(ac_groups)))
      for (j in groups_to_add) {
        new_f_list[[i]] <- f_groups[[j]]
        new_ac_list[[i]] <- ac_groups[[j]]
        i <- i + 1
      }
    }
    # Select only the best feature groups.
    max_ac <- max(unlist(new_ac_list))
    max_groups <- new_f_list[which(unlist(new_ac_list) == max_ac)]
    if (max_ac <= prev_max_ac) {
      # Stop criteria meet, as no improvement found.
      if (verbose) {
        cat(paste("No improvement found at iteration:", iter, "\n"))
        cat(paste("\tSelected features:",
                  paste(prev_max_group, collapse = ", "), "\n"))
        cat(paste("\tMax accuracy:", prev_max_ac, "\n"))
      }
      # Return previous saved copy of the best group (at random if there is more
      # than one).
      return(prev_max_group)
    }
    # Setting values for the next iteration.
    prev_max_ac <- max_ac
    if (length(max_groups) > 1) {
      # Randomly select one best, just in case the next iteration will be the
      # last and a group to return is needed.
      r_num <- sample(1:length(max_groups), 1)
      prev_max_group <- max_groups[[r_num]]
    } else {
      prev_max_group <- max_groups[[1]]  # the only element.
    }
    if (verbose) {
      cat(paste("Iteration number:", iter, "\n"))
      cat(paste("\tFeature groups selected:", length(max_groups) , "\n"))
      cat(paste("\tLength of feature groups:", length(max_groups[[1]]), "\n"))
      cat(paste("\tMax accuracy:", prev_max_ac, "\n"))
    }
  }
  if (verbose) {
    cat("No stop condition was found.\n")
    cat(paste("\tSelected features:",
              paste(prev_max_group, collapse = ", "), "\n"))
    cat(paste("\tMax accuracy:", prev_max_ac, "\n"))
  }
  # If after n iterations stop condition was not meet, the previous best group
  # will be returned.
  return(prev_max_group)
}

#' Boruta Algorithm for Feature Selection
#'
#' This function implements a Boruta algorithm for feature selection. See
#' \code{Boruta::Boruta} documentation for further details.
#'
#' @param input_data [data.frame] a \code{data.frame} rows represent samples and
#'   columns represent features (with beta values or methylation ratios).
#' @param sample_labels [factor] a \code{factor} of length == nrow(input_data),
#'   containing the sample labels.
#' @param n [numeric] the max number of features to output. It is also the
#'   number of features in each GA "chromosome".
#' @param doTrace [numeric] controls \code{Boruta} verbosity. 0 means no
#'   tracing, 1 means reporting decision about each attribute as soon as it is
#'   justified, 2 means same as 1, plus reporting each importance source run.
#'   Default = 0.
#'
#' @return [character] A character vector with \code{n} selected features.
#'
#' @references Miron B. Kursa, Witold R. Rudnicki (2010). Feature Selection with
#'  the Boruta Package. Journal of Statistical Software, 36(11), p. 1-13.
#'  URL: \url{http://www.jstatsoft.org/v36/i11/}
#'
#' @examples
#' set.seed(1234)
#' input_data <- data.frame(a = runif(10, min = 0, max = 1),
#'                          b = runif(10, min = 0, max = 1),
#'                          c = runif(10, min = 0, max = 1),
#'                          d = runif(10, min = 0, max = 1))
#' sample_labels <- as.factor(c(rep("classA", 5), rep("classB", 5)))
#' features_selected <- boruta_fs(input_data, sample_labels, 2)
#' print(features_selected)
#'
#' @export
boruta_fs <- function(input_data, sample_labels, n, doTrace = 0) {
  # First, check whether n is < total number of features.
  if (dim(input_data)[2] <= n) {
    warning("The number of features to select is not smaller than the total
            number of features. No features were selected.")
    return(colnames(input_data))
  } else {
    # Run Boruta algorithm.
    boruta_train <- Boruta::Boruta(x = input_data, y = sample_labels,
                                   doTrace = doTrace)
    # Try to fix tentatives, if any.
    boruta_results <- Boruta::attStats(boruta_train)
    if ("Tentative" %in% boruta_results$decision) {
      boruta_train <- Boruta::TentativeRoughFix(boruta_train)
      boruta_results <- Boruta::attStats(boruta_train)
    }
    # Sort features by meanImp.
    boruta_results <- boruta_results[order(boruta_results$meanImp,
                                           decreasing = TRUE), ]
    # Select all confirmed important attributes, ordered by meanImp.
    selection <- rownames(boruta_results)[boruta_results$decision=="Confirmed"]
    # Return a max of n features.
    if (length(selection) > n) {
      return(selection[1:n])
    } else {
      return(selection)
    }
  }
}

#===============================================================================

#' Prints available feature selection methods
#'
#' This function returns a list with all the available feature selection
#' methods. If \code{verbose = TRUE}, a summary list of the available methods is
#' also printed.
#'
#' @param verbose [logical] if TRUE, a summary list of the available methods is
#'   returned. Default: FALSE.
#'
#' @return A list of three \code{character} vectors is returned: \code{all},
#'   with all the available feature selection methods; \code{filters}, with all
#'   the filter methods; \code{wrappers}, with all the wrapper methods.#'
#'
#' @examples
#' available_fs_methods()
#'
#' @export
available_fs_methods <- function(verbose = FALSE) {
  if (verbose) {
    cat("Available Feature Selection Methods.
Filter methods:
\t.random_fs: Random Feature Selection. For testing purposes.
\t.anova_fs: ANOVA feature selection.
\t.limma_fs: Select CpGs by differential methylation using limma. Similar to
\t           anova_fs but optimized for methylation arrays. Do not use with
\t           bisulfite seq.
\t.correlation_based_fs: Correlation-based Feature Selection. From FSelector
\t                       package.
\t.information_gain_fs: Entropy-based information gain feature selection. From
\t           FSelector package.
\t.relief_fs: RELIEF algorithm for feature selection. From FSelector.
\t.mRMR_fs: Minimum redundancy, maximum relevance feature filtering. Classic
\t          version from mRMRe package. Due to implementation details, the
\t          number of features should be <= 46340.
Wrapper methods:
\t.ga_fs: Wrapped Genetic Algorithm - knn/C5.0/lda for feature selection.
\t.mseq_fs: Modified Sequential Feature Selection.
\t.boruta_fs: Boruta Algorithm for Feature Selection. From Boruta package.\n")
  } else {
    filters <- c("random_fs", "anova_fs", "limma_fs", "correlation_based_fs",
                 "information_gain_fs", "relief_fs", "mRMR_fs")
    wrappers <- c("ga_fs", "mseq_fs", "boruta_fs")
    all <- c(filters, wrappers)
    return(list(all = all, filters = filters, wrappers = wrappers))
  }
}
