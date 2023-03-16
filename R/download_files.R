#' Download Genome
#'
#' Download the planmine Schmidtea mediterranea flatworm genome (SMESG.1) pre-formatted with new lines.
#'
#' @return Returns a .fa file containing genome sequence.
#' @export
EZgenome <- function() {
  gitlink <- "https://github.com/HedonicHotspot/EZworm/raw/master/genomeYAI.fa.zip"
  getOption('timeout')
  options(timeout=1500)
  download.file(gitlink, "genomeYAI.fa.zip")
  unzip("genomeYAI.fa.zip")
  file.remove("genomeYAI.fa.zip")
}

#' Acquire Gene Predictions
#'
#' Match pre-formatted saf file RNASEQ data (SMESG-high confidence release 1) that complements the genome from EZgenome( ).
#'
#' @return This function does not return any values.
#' @export
EZgenePredict <- function() {
  gitlink <- "https://github.com/HedonicHotspot/EZworm/raw/master/smes_v2_repeatfil_YAI.saf"
  download.file(gitlink, "smes_v2_repeatfil_YAI.saf")
}

#' Download Sample Data
#'
#' Download sample pair-read data (fastq) that is shortened for convenience.
#'
#' @export
EZsampleData <- function() {
  gitlinks <- c("https://github.com/HedonicHotspot/EZworm/raw/master/RNAseqData/C71_S76_L006_R1_001.fastq.gz",
                "https://github.com/HedonicHotspot/EZworm/raw/master/RNAseqData/C71_S76_L006_R2_001.fastq.gz",
                "https://github.com/HedonicHotspot/EZworm/raw/master/RNAseqData/P42_S71_L006_R1_001.fastq.gz",
                "https://github.com/HedonicHotspot/EZworm/raw/master/RNAseqData/P42_S71_L006_R2_001.fastq.gz")
  FileNames <- c("C71_S76_L006_R1_001.fastq.gz", "C71_S76_L006_R2_001.fastq.gz",
                 "P42_S71_L006_R1_001.fastq.gz", "P42_S71_L006_R2_001.fastq.gz")
  for (ind in 1:length(gitlinks)) {
    download.file(gitlinks[ind], FileNames[ind])
  }
  dir.create("RNAseqData")
  file.copy(FileNames, "RNAseqData")
  file.remove(FileNames)
  csvlink <- "https://raw.githubusercontent.com/HedonicHotspot/EZworm/master/EZfastaNames.csv"
  download.file(csvlink, "EZfastaNames.csv")
}

#' Download gene annotations
#'
#' The function downloads the Schmidtea mediterranea gene annotations. The gene symbols consist of Augustus and Blast predictions, respectively.
#'
#' @return Downloads the file "PlAnnotation.RDS" to the working directory.
#'
#' @export
EZannotations <- function() {
  gitlink <- "https://github.com/HedonicHotspot/EZworm/raw/master/PlAnnotation.RDS"
  download.file(gitlink, "PlAnnotation.RDS")
}

#' Generate an Index
#'
#' Create an index for the reference genome.
#'
#' @details If the reference genome is not present in the current working directory, it runs
#' the EZgenome() function to download it. The index is then built using the
#' reference genome. This index can be used for further analysis of the
#' reference genome by aiding in retrieval of desired sequences.
#'
#' (Note: This step will take 15-20 minutes.)
#'
#' @export
EZindex <- function() {
  if ("genomeYAI.fa" %in% list.files() == F){EZgenome()}
  genome <- paste0(getwd(), '/genomeYAI.fa')
  buildindex(basename="reference_index", reference= genome)
}
