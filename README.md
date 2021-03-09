
# Methylearning package
 
* Author: David Piñeyro
* License: GPL-3
* Date: 2018-06-01

The **methylearning** package is designed to provide a comprehensive yet easy-to-use framework to test several feature selection and classification methods. This package was primarily designed for Illumina DNA methylation arrays and Whole Genome Bisulfite Sequencing data, but any dataset composed of continuos values ranged between 0 and 1 should work. **VERY IMPORTANT:** no missing values (aka NAs) are allowed, so preprocess your input data to remove samples with missing values, or impute them.

# Installation

This package is created for R 3.4.4 or above. It make use of many [Bioconductor](http://bioconductor.org/) and [CRAN](https://cran.r-project.org/) packages. Make sure that you have upgraded Bioconductor to newest version (3.6). You can upgrade Bioconductor from `R` using:

```{r}
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()
```

Note that for a correct Bioconductor installation, additional system packages may be required. For instance, in Ubuntu, at least the following system packages should be installed before Bioconductor installation (root access required).

```
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
```

After you have R and Bioconductor installed properly, you can download `methylearning` from this repositoy. The installation can be achieved using `devtools::install` function in your R session:

```{r}
# Set one level up methylearning folder.
setwd("path/to/folder/containing/methylearning/folder")
devtools::install("methylearning")
```

Some errors may occur during the installation. Most of them may be due to errors loading recursively depending R packages and/or required system packages. If you got an error, please read carefully the error log which usually already suggests some solution. In Ubuntu 17.10, the package requires that Java was correctly installed and configured to work with R and also requires the `mysqlclient library` to install `RMySQL` package. To do so, in the command shell, type (it requires root access):

```
sudo apt-get update
sudo apt-get upgrade
# install and configure Java.
sudo apt-get install default-jre default-jdk
sudo R CMD javareconf
# install MySQL server.
sudo apt-get install libmariadbclient-dev
```

After installation, load the package in your R session using:

```{r}
library(methylearning)
```

`methylearning` package make use of many other packages. It is possible that you reach your maximum number of allowed DLLs in your `R` session (i.e. more than 100). If you experience any error realated to a failed loaded package please, start `methylearning` from a clean `R` session. If you still experience problems, type the code shown below in your command shell (Linux or similar OS) to increase your DLL limit.

```
# Create an .Renviron in your home directory, which contains the 
# new value for the R_MAX_NUM_DLLS variable (it will be increased to 150).
echo "R_MAX_NUM_DLLS=150" > ~/.Renviron
```