% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/organization.R
\name{EZlist2csv}
\alias{EZlist2csv}
\title{Convert List of Data Frames to a Single CSV File}
\usage{
EZlist2csv(DFlist, File = "DFlist.csv")
}
\arguments{
\item{DFlist}{List of data frames to be combined and written to a single CSV file.}

\item{File}{Name of the CSV file to be created. Default is 'DFlist.csv'.}
}
\value{
The function returns a CSV file with each data frame in the input list written as a separate sheet within the file.
}
\description{
Take a list of data frames and write each data frame into a single comma separated values (CSV) file.
}
\details{
Note: Useful for reporting GO results of several samples and genes grouped by functions.
}
\examples{
DFlist <- list(mtcars, iris)
EZlist2csv(DFlist)

}
