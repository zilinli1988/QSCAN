# QSCAN 
This is an R package for performing Q-SCAN procedure in whole genome sequencing studies.
## Description
QSCAN is an R package for performing a computationally efficient quadratic scan (Q-SCAN) statistic based method to 
detect the existence and the locations of signal regions by  scanning the genome continuously. 
The proposed method accounts for the correlation (linkage disequilibrium) among genetic variants, 
and allows for signal regions to have both causal and neutral variants, 
and the effects of signal variants to be in different directions.
## Dependencies
QSCAN links to R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a> and 
<a href="https://cran.r-project.org/web/packages/RcppArmadillo/index.html">RcppArmadillo</a>, 
and also imports R packages <a href="https://cran.r-project.org/web/packages/Rcpp/index.html">Rcpp</a>, 
<a href="https://cran.r-project.org/web/packages/RcppArmadillo/index.html">RcppArmadillo</a> and
<a href="https://cran.r-project.org/web/packages/Matrix/index.html">Matrix</a>, 
These dependencies should be installed before installing QSCAN.
## Installation
```
library(devtools)
devtools::install_github("zilinli1988/QSCAN")
```
## Usage
Please see the <a href="docs/QSCAN.pdf">**QSCAN** user manual</a> for detailed usage of QSCAN package. 
Please see the <a href="docs/QSCAN_Example.pdf">**QSCAN** Example</a> for an example of analyzing sequencing data using Q-SCAN procedure.

## License
This software is licensed under GPLv3.
