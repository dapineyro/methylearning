#===============================================================================
# Functions to create ml_data objects from GSE objects (see GEOquery::getGEO).
# Author: David Piñeyro
# License: GPL-3
# Date: 2018-04-22
#===============================================================================
# Golbal variables
globalVariables(c("crossreactive27k", "crossreactive450k", "crossreactiveEPIC",
                  "polymorphic27k", "polymorphic450k", "polymorphicEPIC",
                  "probesxy27k", "probesxy450k", "probesxyEPIC"))
#' Detect Illumina Array platform
#'
#' This function returns either "27k", "450k", "EPIC" or NULL, depending on
#' the correct platform detection.
#'
#' @param gse [GSE] a GSE object created with getGEO, containing an Illumina
#'   array dataset.
#'
#' @return [character] A character vector containing either "27k", "450k" or
#'   "EPIC", depending on the correct platform detection or NULL if platform is
#'   not one of those.
#'
detect_platform <- function(gse) {
  #   Allowed GPLs are:
  # GPL8490: Illumina HumanMethylation 27K
  # GPL13534: Illumina HumanMethylation450_15017482
  # GPL16304: Illumina HumanMethylation450 BeadChip_UBC enhanced annotation v1.0
  # GPL18809: Illumina HumanMethylation450 BeadChip (v1.2, extended annotation)
  # GPL21145: Infinium MethylationEPIC
  #   Pending of testing for annotation compatibility:
  # GPL23976: Illumina Infinium HumanMethylation850 BeadChip
  platform_id <- unique(gse$platform_id)
  # Controlling only one platform per series.
  if (length(platform_id) > 1) {
    stop(paste("More than one platform detected:",
               paste(platform_id, collapse = "; ")))
  }
  if (platform_id == "GPL8490") {
    return("27k")
  } else if (platform_id %in% c("GPL13534", "GPL16304", "GPL18809")) {
    return("450k")
  } else if (platform_id == "GPL21145") {
    return("EPIC")
  } else {
    return(NULL)
  }
}

#===============================================================================

#' Get Sex Chromosome Probes
#'
#' From a valid platform, this function returns a character vector with all the
#' probes located at sex chromosomes, based on platform annotation.
#'
#' @param platform [character] a character vector of lenght 1, with either
#'   "27k", "450k" or "EPIC".
#'
#' @return [character] A character vector all those probes located at sex
#'   chromosomes.
#'
get_chrXY <- function(platform) {
  if (platform == "27k") {
    return(probesxy27k)
  } else if (platform == "450k") {
    return(probesxy450k)
  } else if (platform == "EPIC") {
    return(probesxyEPIC)
  } else {
    warning("Illumina Methylation platform was not detected. Sex chromosome probes not removed.")
    return(NULL)
  }
}

#===============================================================================

#' Select CpGs by differential methylation
#'
#' This funtion uses limma to select top differentially methylated probes
#' between groups (of class labels factor). The number of classes can be 2 or
#' more, but no paired designs are yet implemented.
#'
#' @param betas [matrix] a matrix with normalized beta-values, where rows are
#'   probes and columns are samples.
#' @param class_labels [factor] a factor with sample class labels.
#' @param limma_select_best [numeric] the number of top differentially
#'   methylated CpGs to select.
#'
#' @return [character] A character vector with the \code{limma_select_best}
#'   top differentially methylated (sorted by p-value) probe names. If there
#'   are more than 2 groups, the same number of top CpGs are returned, but
#'   divided between each one vs the others contrast.
#'
apply_limma_filter <- function(betas, class_labels, limma_select_best) {
  # Calculate M-values using M=log2(Beta/(1-Beta)). All statistics will be
  # performed on M-values.
  m_vals <- log2(betas / (1-betas))
  # Create a targets dataframe.
  sample_names <- colnames(betas)
  pheno_data <- data.frame(sample_names, class_labels)
  rownames(pheno_data) <- sample_names
  targets <- stats::model.frame(sample_names ~ class_labels, pheno_data)
  # Design matrix (only unpaired test supported).
  design <- stats::model.matrix(~0+class_labels, data=targets)
  colnames(design) <- levels(class_labels)
  # Contrast matrix (one vs the others).
  if (nlevels(class_labels) == 2) {
    contr <- paste0(levels(class_labels)[1], "-", levels(class_labels)[2])
    contMatrix <- limma::makeContrasts(contrasts = contr, levels = design)
  } else {
    # More than 2 groups. Using One Vs the Others contrasts.
    i <- 1
    contr <- character()
    while (i <= nlevels(class_labels)) {
      one <- levels(class_labels)[i]
      the_others <- levels(class_labels)[-i]
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
  list_of_results <- lapply(1:(nlevels(class_labels) - 1), function(x) {
    limma::topTable(fit2, coef = x, number=Inf, sort.by = "P")
  })
  # Selecting the "best" CpGs of each contrast. If
  # limma_select_best/nlevels(class_lables) has no integer result, round
  # approximation is taken.
  each_contrast_n <- round(limma_select_best/(nlevels(class_labels) - 1), 0)
  gene_selection <- character()
  for (r in list_of_results) {
    gene_selection <- c(gene_selection, rownames(r)[1:each_contrast_n])
  }
  return(gene_selection)
}

#===============================================================================

#' Prepare Methylation Array Data From GEO
#'
#' This function uses a \code{GSE} object (a \code{list} object) created with
#' \link[GEOquery]{getGEO}, containing only methylation array data,  and process
#' it to create an \code{ml_data} object.
#'
#' @details By default data is preprocessed to remove sex chromosome,
#'   cross-reactive (multi-mapping) and polymorphic (with known SNP) probes.
#'   Moreover, one of th following filters can be applied to reduce the number
#'   of features (probes): limma_filter, based on limma statistical
#'   significance; random_filter, selecting random probes.
#'
#' @param gse [GSE] an object of \code{GSE} class, created using
#'   \code{GEOquery::getGEO}.
#' @param target_name [character] a character vector of the target (class
#'   labels) column name in \code{gse phenoData data.frame}.
#' @param exclude_sex_probes [logical] should chrX and chrY probes be removed?
#'   Default: TRUE.
#' @param exclude_cross_reactive [logical] should cross-reactive probes be
#'   removed? Default: TRUE.
#' @param exclude_polymorphic [logical] should polymorphic (i.e. probes with
#'   known SNP at CpG position) be removed? Default: TRUE.
#' @param limma_filter [logical/numeric] should limma filtering be performed?
#'   Default: FALSE. If an integer is specified, this will be the number of
#'   features selected.
#' @param random_filter [logical/numeric] should random filtering be performed?
#'   Default: FALSE. If an integer is specified, this will be the number of
#'   features selected
#'
#' @return [ml_data] an object of ml_data class.
#'
#' @examples
#' \dontrun{
#' if (require(GEOquery) &
#'     require(Biobase)) {
#'   gse <- getGEO("GSE69229")
#'   ml <- get_GEO_methylarray(gse = gse, target_name = "outcome:ch1")
#'   summary(ml)
#' }
#' }
#'
#' @include ml_data.R
#' @export
get_GEO_methylarray <- function(gse, target_name,
                                exclude_sex_probes = TRUE,
                                exclude_cross_reactive = TRUE,
                                exclude_polymorphic = TRUE,
                                limma_filter = FALSE,
                                random_filter = FALSE) {
  # First detect number of series present in the experiment. As series
  # usually imply different technologies, more than one series present
  # (Superseries) is not allowed.
  if (length(gse) > 1) {
    stop("Error: your experiment has more than one series. Please, select only one series from this superseries (i.e. using slicing operator in your superseries list) and re-run this command.")
  }
  # Check only one filter is indicated.
  if (limma_filter & random_filter) {
    stop("Error: limma and random filters cannot be selected at the same time.")
  }
  gse <- gse[[1]]  # Select the only series.
  # Detect Illumina platform.
  platform <- detect_platform(gse)
  if (exclude_sex_probes) {
    # Generate a character vector with probe names of sex-chromosome probes.
    sex_chr_probes <- get_chrXY(platform)
    keep <- !(minfi::featureNames(gse) %in% sex_chr_probes)
    # Remove sex-chromosome probes.
    gse <- gse[keep, ]
  }
  # Remove cross-reactive probes.
  if (exclude_cross_reactive) {
    if (platform == "27k") {
      cr_probes <- crossreactive27k
    } else if (platform == "450k") {
      cr_probes <- crossreactive450k
    } else if (platform == "EPIC") {
      cr_probes <- crossreactiveEPIC
    } else {
      warning("Illumina Methylation platform was not detected. Cross-reactive probes not removed.")
    }
    keep <- !(minfi::featureNames(gse) %in% cr_probes)
    gse <- gse[keep, ]
  }
  # Remove polymorphic probes (with known SNPs in CpG site).
  if (exclude_polymorphic) {
    if (platform == "27k") {
      pm_probes <- polymorphic27k
    } else if (platform == "450k") {
      pm_probes <- polymorphic450k
    } else if (platform == "EPIC") {
      pm_probes <- polymorphicEPIC
    } else {
      warning("Illumina Methylation platform was not detected. Polymorphic probes not removed.")
    }
    keep <- !(minfi::featureNames(gse) %in% pm_probes)
    gse <- gse[keep, ]
  }
  # Get normalized beta-values.
  betas <- Biobase::exprs(gse)
  # Get class labels.
  class_labels <- as.factor(
    gse@phenoData@data[, which(colnames(gse@phenoData@data) == target_name)])
  # Make syntactically valid names for labels.
  levels(class_labels) <- make.names(levels(class_labels))
  # Check limma_filter or random_filter are in range.
  if (limma_filter) {
    if (limma_filter > dim(betas)[1]) {
      stop("Error: features indicated to keep by limma filter exceed total number of features")
    }
  }
  if (random_filter) {
    if (random_filter > dim(betas)[1]) {
      stop("Error: features indicated to keep by rando filter exceed total number of features")
    }
  }
  # Perform limma filtering.
  if (limma_filter) {
    cpg_selection <- apply_limma_filter(betas, class_labels, limma_filter)
    gse <- gse[cpg_selection, ]
  }
  # Perform random filtering.
  if (random_filter) {
    cpg_selection  <- sample(dim(betas)[1], random_filter)
    gse <- gse[cpg_selection, ]
  }
  # Generate final data.frame.
  betas_final <- Biobase::exprs(gse)
  input_data <- as.data.frame(t(betas_final))
  input_data <- cbind(class_labels, input_data)
  # Generate ml_data object.
  ml_data_object <- ml_data(input_data, labels_column = 1)

  return(ml_data_object)
}

#' Documentation for datasets used by get_GEO_methylarray
#'
#'
#' @name crossreactive27k
#' @docType data
#' @keywords data
NULL

#' Documentation for datasets used by get_GEO_methylarray
#'
#'
#' @name crossreactive450k
#' @docType data
#' @keywords data
NULL

#' Documentation for datasets used by get_GEO_methylarray
#'
#'
#' @name crossreactiveEPIC
#' @docType data
#' @keywords data
NULL

#' Documentation for datasets used by get_GEO_methylarray
#'
#'
#' @name polymorphic27k
#' @docType data
#' @keywords data
NULL

#' Documentation for datasets used by get_GEO_methylarray
#'
#'
#' @name polymorphic450k
#' @docType data
#' @keywords data
NULL

#' Documentation for datasets used by get_GEO_methylarray
#'
#'
#' @name polymorphicEPIC
#' @docType data
#' @keywords data
NULL

#' Documentation for datasets used by get_GEO_methylarray
#'
#' 27k chrX and chrY probes based on:
#' IlluminaHumanMethylation27k.db::IlluminaHumanMethylation27kCHR
#'
#' @name probesxy27k
#' @docType data
#' @keywords data
NULL

#' Documentation for datasets used by get_GEO_methylarray
#'
#' 450k chrX and chrY probes based on:
#' IlluminaHumanMethylation450kanno.ilmn12.hg19
#'
#' @name probesxy450k
#' @docType data
#' @keywords data
NULL

#' Documentation for datasets used by get_GEO_methylarray
#'
#' EPIC chrX and chrY probes based on:
#' IlluminaHumanMethylationEPICanno.ilm10b4.hg19
#'
#' @name probesxyEPIC
#' @docType data
#' @keywords data
NULL
