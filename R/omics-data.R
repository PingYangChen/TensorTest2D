
#' Lung-cancer cell lines data in cancer cell line encyclopedia (CCLE) dataset
#'
#' The omics data is a subset of the dataset provided by cancer cell line 
#' encyclopedia (CCLE) project (Barretina et al., 2012; 
#' \url{https://portals.broadinstitute.org/ccle/about}{https://portals.broadinstitute.org/ccle/about}).
#' This data consists of one response variable and ten genes evaluated under three different platforms. 
#' The response variable measures the log-transformed activity area of taking Vandertanib, 
#' a drug targeting on EGFR gene for lung cancer. 
#' The three platforms are DNA copy number variation (CNV), methylation and mRNA expression. 
#' Among the 10 genes, 7 of them (EGFR, EREG, HRAS, KRAS, PTPN11, STAT3, and TGFA) are involved 
#' in the protein-protein interaction network of EGFR (\href{https://string-db.org}{https://string-db.org}) 
#' and the rest (ACTB, GAPDH, and PPIA) are arbitrarily chosen housekeeping genes and 
#' play the role of negative control. 
#' Detailed pre-processing procedure is available in Chang et al. (2021).  
#'
#' @docType data
#'
#' @usage data(omics)
#'
#' @format A list of two objects, \code{omics$omics} and \code{omics$Y}.
#' \code{omics$omics} is a 3-dimensional array of size (3, 10, 68).
#' \code{omics$Y} is a vector of length 68 representing the response variable.
#'
#' @keywords datasets
#'
#' @references Barretina, J., Caponigro, G., Stransky, N. et al. 
#' The Cancer Cell Line Encyclopedia enables predictive modelling of anticancer drug sensitivity. Nature 483, 603–607 (2012).
#' (\href{https://doi.org/10.1038/nature11003}{URL})
#' 
#' @references Sheng-Mao Chang, Meng Yang, Wenbin Lu, Yu-Jyun Huang, Yueyang Huang, Hung Hung, 
#' Jeffrey C Miecznikowski, Tzu-Pin Lu, Jung-Ying Tzeng, 
#' Gene-set integrative analysis of multi-omics data using tensor-based association test, 
#' Bioinformatics, 2021;, btab125, (\href{https://doi.org/10.1093/bioinformatics/btab125}{URL})
#'
#' @source \href{https://string-db.org}{https://string-db.org}
#'
#' @examples
#' data(omics)
#' names(omics)
#' dim(omics$omics)
#' # 3 10 68
"omics"


