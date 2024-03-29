% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align.R
\name{EZalign}
\alias{EZalign}
\title{Align RNA-seq data}
\usage{
EZalign(NameCSV = FALSE)
}
\arguments{
\item{NameCSV}{A character string indicating the file name of a CSV file
containing the original and renamed names of the RNA-seq data files.}
}
\value{
Returns the gene-level read counts and saves the result as a RDS file.
}
\description{
Perform alignment and generate counts (Rsubread) of RNA-seq data to the reference genome
and annotation of genes to obtain gene-level read counts.
}
\details{
The function checks for the presence of a reference index, reference genome, and gene annotations and creates them if they do
not exist. It is recommended to use the EZworm functions beforehand for troubleshooting purposes. The function also checks for the existence of RNA-seq data files. The function performs paired-end or single-end
alignment depending on the availability of the data. To prepare your fastq files for analysis,
download all files to the RNAseqData folder within the RNAseqAlignment folder. Ensure that you
have sufficient storage capacity as fastq files are large. If your data is paired-end, the file
names should contain "_R1.fa" or "_R2.fa" for the respective files. Avoid adding "." elsewhere.
Additionally, make sure the file names indicate the animal experimental condition and timepoint.
For example, the format could be Cond1_T1_R1.fa. If your data is single-end, ensure the file names
contain the .fa extension, animal experimental condition and timepoint too.
For example, the format could be Cond1_T1.fa.
}
\examples{
EZalign()
EZalign(NameCSV = "EZfastaNames.csv")

}
