% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml_fs.R
\name{selection_df}
\alias{selection_df}
\title{Create a \code{data.frame} with feature selection results}
\usage{
selection_df(fs_object)
}
\arguments{
\item{fs_object}{[ml_fs] an object of \code{ml_fs} class.}
}
\value{
A \code{data.frame} where unique features selected are rows and
  feature selection methods used are columns. TRUE/FALSE indicated whether a
  mathod selected a feature. The last column is called \code{agree} and
  counts the number of methods that agreed in each feature selection.
  The returned \code{data.frame} is deacreasing ordered by \code{agree}
  column, so most selected features are first.
}
\description{
This function creates a \code{data.frame} with the feature selection results.
}
\examples{
set.seed(1234)
df <- data.frame(a = runif(10, 0, 1),
                 b = runif(10, 0, 1),
                 c = as.factor(as.integer(runif(10, 0, 2))))
data <- ml_data(df, 3)
data_sel <- ml_fs(data, fs_methods = "filters", selection_size = 2)
df <- selection_df(data_sel)
print(df)
# Select the most selected feature.
rownames(df)[1]

}
