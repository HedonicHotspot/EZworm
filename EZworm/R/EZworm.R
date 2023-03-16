#' @name EZworm
#' @title EZworm: RNA-seq Alignment for the Planarian Flatworm
#' @docType package
#'
#' @description
#' EZworm is an R package for analyzing/aligning RNA-seq data for the Planarian genome. It simplifies the RNA-seq data analysis process for researchers by providing a user-friendly interface and efficient performance. The package includes various functionalities, such as downloading necessary data such as genomic sequences and gene predictions, as well as providing options for manipulating and organizing the data.
#' @author
#'   Maintainer: Yacoub Innabi | yinnabi@ucmerced.edu
#'
#'   Author:
#'   - Yacoub Innabi | yinnabi@ucmerced.edu
#'
#'   Contributor:
#'   - Ju-Won Lee | juwon7lee@gmail.com
#'
#'
#' @seealso
#' Useful links:
#' - https://github.com/HedonicHotspot/EZworm
#' - Report bugs at https://github.com/HedonicHotspot/EZworm/issues
#'
#' @examples
#' #Download necessary files
#' EZgenome()
#' EZgenePredict()
#' EZsampleData()  #For demonstration
#' EZannotations()
#'
#' #Build Index
#' EZindex()
#'
#' #Align and Summarize
#' EZalign(NameCSV = "EZfastaNames.csv")
#' EZbamSummary()
#'
#' #Generate transcript counts, convert to gene counts, and annotate
#' rawcounts <- read_RDS("rawcounts.RDS") |> EZsumTrans() |> EZcol2namerow()
#' rawcounts <- EZgene2sym(rawcounts, IDinRow = T)
NULL





