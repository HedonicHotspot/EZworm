#' Align RNA-seq data
#'
#' Perform alignment and generate counts (Rsubread) of RNA-seq data to the reference genome
#' and annotation of genes to obtain gene-level read counts.
#'
#' @param NameCSV A character string indicating the file name of a CSV file
#' containing the original and renamed names of the RNA-seq data files.
#'
#' @return Returns the gene-level read counts and saves the result as a RDS file.
#'
#' @examples
#' EZalign()
#' EZalign(NameCSV = "EZfastaNames.csv")
#'
#' @details The function checks for the presence of a reference index, reference genome, and gene annotations and creates them if they do
#' not exist. It is recommended to use the EZworm functions beforehand for troubleshooting purposes. The function also checks for the existence of RNA-seq data files. The function performs paired-end or single-end
#' alignment depending on the availability of the data. To prepare your fastq files for analysis,
#' download all files to the RNAseqData folder within the RNAseqAlignment folder. Ensure that you
#' have sufficient storage capacity as fastq files are large. If your data is paired-end, the file
#' names should contain "_R1.fa" or "_R2.fa" for the respective files. Avoid adding "." elsewhere.
#' Additionally, make sure the file names indicate the animal experimental condition and timepoint.
#' For example, the format could be Cond1_T1_R1.fa. If your data is single-end, ensure the file names
#' contain the .fa extension, animal experimental condition and timepoint too.
#' For example, the format could be Cond1_T1.fa..
#'
#' @export

EZalign <- function(NameCSV = FALSE) {
  refIndexFiles <- list.files(pattern="reference_index")
  if (length(refIndexFiles)!=5){
    EZindex()
  }


  if (class(NameCSV) == "character"){
    if(prod(read.csv("EZfastaNames.csv", header = F)[,2] %in%
            list.files("RNAseqData")) == 0 ) {
              file.rename(
                paste0("RNAseqData/",read.csv(NameCSV, header = F)[,1]),
                paste0("RNAseqData/",read.csv(NameCSV, header = F)[,2])
              )
    
    }
  }
  

  gene_predict <- paste0(getwd(),'/smes_v2_repeatfil_YAI.saf')
  RNAseqData <- paste0(getwd(), "/RNAseqData/")
  fastqfiles <- list.files(path = "RNAseqData", full.names = TRUE)
  fastqfiles <- grep(".fa", fastqfiles, value = TRUE)


  if (sum(grepl("R1.fa|R2.fa", ignore.case = T, fastqfiles)) > 2) {
    if ("smes_v2_repeatfil_YAI.saf" %in% list.files() == F) {EZgene_predict()}

    fastqfiles1 <- grep("R1.fa", fastqfiles, value = TRUE)
    fastqfiles1 <- sort(fastqfiles1)
    fastqfiles2 <- grep("R2.fa", fastqfiles, value = TRUE)
    fastqfiles2 <- sort(fastqfiles2)

    BAMfilenames <- gsub("RNAseqData/|_R1.fasta.gz|_R1.fa", "",
                         fastqfiles1, ignore.case = T)

    if ( prod(BAMfilenames %in% list.files()) == 0) {
      align(index="reference_index",
            readfile1=fastqfiles1,
            readfile2 = fastqfiles2,
            type="rna",
            output_file = BAMfilenames,
            isGTF = TRUE)
    }


    rawcounts <- featureCounts(BAMfilenames,
                               annot.ext = read.table(gene_predict, header = T),
                               isPairedEnd = TRUE,
                               allowMultiOverlap = TRUE,
                               fraction = TRUE)

    write_rds(rawcounts, "rawcounts.RDS")
  }


  else {

    BAMfilenames <- gsub("RNAseqData/|.fasta.gz|.fasta|.fa|.fa.gz", "",
                         fastqfiles, ignore.case = T)
    if ( prod(BAMfilenames %in% list.files()) == 0) {

      align(index="reference_index",
            readfile1=fastqfiles,
            type="rna",
            output_file = BAMfilenames,
            isGTF = TRUE)
    }

    rawcounts <- featureCounts(BAMfilenames,
                               annot.ext = read.table(gene_predict, header = T),
                               isPairedEnd = FALSE,
                               allowMultiOverlap = TRUE,
                               fraction = T)

    write_rds(rawcounts, "rawcounts.RDS")
  }
}
