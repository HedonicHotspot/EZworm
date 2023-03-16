#' Bam File Summary
#'
#' Summarize the alignment statistics for a set BAM files
#'
#' @details The function EZBamSummary() reads in summary files with a ".summary"
#' extension in the working directory and converts them into a neat data frame.
#' The resulting data frame is then saved as a CSV file with the file name
#' "align_summary.csv".
#'
#' @return A csv file named "align_summary.csv" is created in the working directory containing the alignment summary information.
#'
#' @export
EZbamSummary <- function() {
  bam.files.summary <- list.files(path = getwd(),
                                  pattern = ".summary$",
                                  full.names = TRUE)

  align_summary <- lapply(bam.files.summary, read.table,
                          sep="\t", header=FALSE)
  align_summary <- as.data.frame(align_summary)
  rownames(align_summary) <- align_summary[,1]
  align_summary <- align_summary[ ,seq(2,length(align_summary),2)]
  colnames(align_summary) <- BAMfilenames
  align_sum_file <- paste0(getwd(), "/align_summary.csv", BAMfilenames)
  write_csv(align_summary, align_sum_file)
}

#' Order Heatmap Rows
#'
#' Organize data using z-score.
#'
#' @param df A data frame
#' @param cols A vector of column indices to use for ordering. The default is the first column (1).
#' @param decr A logical indicating whether the ordering should be in decreasing order. The default is FALSE.
#'
#' @return A data frame with the same columns but sorted rows based on the mean and standard deviation of the specified columns.
#'
#' @details The 'EZhm_order' function takes a data frame (df) and orders the rows based on
#' the calculated z-score. The selected
#' columns can be specified with the 'cols' parameter, which is a vector of
#' column indices. The default value is the first column (1). The order of
#' the rows can be either increasing (default) or decreasing, which is specified
#' with the 'decr' parameter. If multiple columns are selected, the function
#' calculates the z-score of each row relative to the entire data frame, averages
#' the z-scores of the selected columns, and sorts the data frame based on these values.
#' If only one column is selected, the function orders the data frame based on the z-scores
#' of that single column.
#'
#' @export
EZhm_order <- function(df, cols = c(1), decr = FALSE){
  if(length(cols) > 1){
    hmz <- apply(df, 2, FUN = function(x)
    {return( (rowMeans(df) - x)/rowSds(df) )} )
    hmrow <- rowSums(hmz[,cols])
    df[order(hmrow, decreasing = decr),]
  }
  else {
    hmz <- apply(df, 2, FUN = function(x)
    {return( (rowMeans(df) - x)/rowSds(df) )} )
    df[order(hmz[,cols], decreasing = decr),]
  }
}

#' Convert Column to Row Names in a Data Frame
#'
#' Convert a specified column of a data frame into row names, and remove the specified column.
#' @param df A data frame object.
#' @param col An integer specifying which column in the data frame to convert to row names. Default is 1.
#'
#' @return The modified data frame.
#'
#' @examples
#' data(df)
#' head(df)
#' EZcol2namerow(df, 1)
#'
#' @export
EZcol2namerow <- function(df, col=1) {
  rownames(df) <- df[,col]
  df <- df[,-col]
  return(df)
}

#' Sum Transcript Counts to Gene Counts
#'
#' Summarize Transcript Counts by Gene ID
#'
#' @param counts A data frame of count data
#'
#' @return A data frame with summed counts for each sample, grouped by 'GeneID', and with row names as column names.
#'
#' @details The EZSumTrans function takes a 'counts' data frame and merges it with a saved data frame 'PlAnnotation.RDS' based on a common column. The merged data is then grouped by the 'GeneID' column and the sum of each sample column is calculated. The resulting data is converted to a data frame and the column names are converted to row names
#'
#' @export
EZSumTrans <- function(counts) {
  merge(counts, read_rds("PlAnnotation.RDS")[,1:2], by.x = 0, by.y = 2) |>
    group_by(GeneID) |>
    summarise(across(colnames(counts),sum))|> as.data.frame() |>
    col2namerow()
}

#' Rename Gene IDs to Gene Symbols
#'
#' Map gene IDs to gene symbols in a data frame.
#'
#' @param df A data frame with gene IDs.
#' @param IDinRow A logical indicating whether the gene IDs are in the row names. The default is TRUE.
#'
#' @return A data frame with added gene symbol column, either next to the gene ID column or in the row names.
#'
#' @details The function takes in a data frame, df, and an optional argument IDinRow. If IDinRow is set to TRUE, the gene IDs are assumed to be in the row names, and the symbol column is added to the data frame. If IDinRow is set to FALSE, the gene IDs are assumed to be in a column, and the symbol column is added next to it.
#'
#' @examples
#' df_mapped <- EZgene2sym(df, IDinRow = FALSE)
#'
#' @export
EZgene2sym <- function (df,IDinRow = T){
  if (IDinRow) {
    df <- merge(df, read_rds("PlAnnotation.RDS")[,c("GeneID", "Symbol")],
                by.x = 0, by.y = "GeneID") |>
      col2namerow("Symbol") |> dplyr::select(-"Row.names")
  }
  else {
    df <- merge(df, read_rds("PlAnnotation.RDS")[,c("GeneID", "Symbol")],
                by = "GeneID") |> dplyr::select(-"GeneID")

  }
  return(df)
}

#' Convert List of Data Frames to a Single CSV File
#'
#' Take a list of data frames and write each data frame into a single comma separated values (CSV) file.
#'
#' @param DFlist List of data frames to be combined and written to a single CSV file.
#' @param File Name of the CSV file to be created. Default is 'DFlist.csv'.
#' @return The function returns a CSV file with each data frame in the input list written as a separate sheet within the file.
#'
#' @details Note: Useful for reporting GO results of several samples and genes grouped by functions.
#'
#' @examples
#' DFlist <- list(mtcars, iris)
#' EZlist2csv(DFlist)
#'
#' @export
EZlist2csv <- function(DFlist, File = 'DFlist.csv') {
  for (ind in 1:length(DFlist)) {
    write.table(paste0(names(DFlist)[ind]), file = File,
                sep = ",", col.names = FALSE, append = T, row.names = FALSE)
    write.table( DFlist[[ind]], file = File,
                 append= T, sep=',', col.names = NA)
  }
}
