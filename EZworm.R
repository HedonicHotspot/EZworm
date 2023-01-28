packages = c("BiocManager","tidyverse", "ggplot2", "dplyr", 
             "ggthemes", "formattable", "gplots", "readr",  "RColorBrewer", 
             "limma", "edgeR", "Rsubread", "topGO", "circlize", "genefilter",
              "scriptName")

## Load R and BioConductoR packages
lapply(packages, library, character.only = T)

EZ <- list()  # Stores functions


# Alignment Functions

# Download all your fastq files to the RNAseqData folder located within the 
# RNAseqAlignment folder. Fastq files are large, so make sure you have enough 
# storage on your device. 
# IF your data is paired-end, the file names need to contain "_R1.fa" or 
#"_R2.fa" for the respective files. Do add a "." anywhere else.
# Make sure your file names also indicate the animal experimental condition, 
# timepoint, etc. For example, follow this format:    
#   Cond1_T1_R1.fa or else the code will not work.
# IF your data is single-end, you just need to make sure your files contain the .fa   
# extension and the animal experimental condition, timepoint, etc. For example,      
# follow this format: Cond1_T1.fa or else the code will not work.

# Download Genome
# Planmine Schmidtea mediterranea S2F2 genome (FASTA format)
# Smed genome formatted with newlines (i.e. "\n")
EZ$DLgenome <- function() {
  gitlink <- "https://github.com/HedonicHotspot/EZworm/raw/master/genomeYAI.fa.zip"
  getOption('timeout')
  options(timeout=1500)
  download.file(gitlink, "genomeYAI.fa.zip")
  unzip("genomeYAI.fa.zip")
  }



# Fetch Gene Predictions
EZ$DLgene_predict <- function() {
  gitlink <- "https://github.com/HedonicHotspot/EZworm/raw/master/smes_v2_repeatfil_YAI.saf"
  write.table(as.data.frame(read_table(gitlink)), "smes_v2_repeatfil_YAI.saf")
}

# Download Schmidtea Mediterranea gene annotations
# Gene symbols consist of Augustus and Blast predictions
EZ$DLannotations <- function() {
  gitlink <- "https://github.com/HedonicHotspot/EZworm/raw/master/PlAnnotation.RDS"
  download.file(gitlink, "PlAnnotation.RDS")
}

EZ$sampleData <- function() {
  gitlinks <- c("https://github.com/HedonicHotspot/EZworm/raw/master/RNAseqData/C71_S76_L006_R1_001.fastq.zip",
                "https://github.com/HedonicHotspot/EZworm/raw/master/RNAseqData/C71_S76_L006_R2_001.fastq.zip",
                "https://github.com/HedonicHotspot/EZworm/raw/master/RNAseqData/P42_S71_L006_R1_001.fastq.zip",
                "https://github.com/HedonicHotspot/EZworm/raw/master/RNAseqData/P42_S71_L006_R2_001.fastq.zip")
  FileNames <- c("RNAseqData/C71_S76_L006_R1_001.fastq.zip", "RNAseqData/C71_S76_L006_R2_001.fastq.zip",
                 "RNAseqData/P42_S71_L006_R1_001.fastq.zip", "RNAseqData/P42_S71_L006_R2_001.fastq.zip")
  for (ind in 1:length(gitlinks)) {
    download.file(gitlinks[ind], FileNames[ind])
    unzip(FileNames[ind])
  }
  csvlink <- "https://raw.githubusercontent.com/HedonicHotspot/EZworm/master/EZfastaNames.csv"
  download.file(csvlink, "EZfastaNames.csv")
}

# For building reference index
EZ$index <- function() {
  if ("genomeYAI.fa" %in% list.files() == F){EZ$DLgenome()}
  genome <- paste0(getwd(), '/genomeYAI.fa')
  buildindex(basename="reference_index", reference= genome)
}

# For aligning RNAseq data
# May input a csv with first column as old RNAseq filenames
# and second column new filenames EZformatted (i.e includes R1.fa/R2.fa)
# Ex. EZ$align(NameCSV = "EZfastaNames.csv")
EZ$align <- function(NameCSV = FALSE) {
  refIndexFiles <- list.files(pattern="reference_index")
  if (length(refIndexFiles)!=5){ # Download index if not present
    EZ$index()
  }
  
  # Rename Fastq Files
  if ( (class(NameCSV) == "character") &
       (prod(read.csv("EZfastaNames.csv", header = F)[,2] %in% 
             list.files("RNAseqData")) == 0 )) {
    file.rename(
      paste0("RNAseqData/",read.csv(NameCSV, header = F)[,1]), 
      paste0("RNAseqData/",read.csv(NameCSV, header = F)[,2])
    )
  }
  
  gene_predict <- paste0(getwd(),'/smes_v2_repeatfil_YAI.saf') 
  RNAseqData <- paste0(getwd(), "/RNAseqData/")
  fastqfiles <- list.files(path = "RNAseqData", full.names = TRUE)
  fastqfiles <- grep(".fa", fastqfiles, value = TRUE)
  
  # For paired-end
  if (sum(grepl("R1.fa|R2.fa", ignore.case = T, fastqfiles)) > 2) {
    if ("smes_v2_repeatfil_YAI.saf" %in% list.files() == F) {EZ$DLgene_predict()}
    #Each pair read file in seperate vectors
    fastqfiles1 <- grep("R1.fa", fastqfiles, value = TRUE)
    fastqfiles1 <- sort(fastqfiles1)
    fastqfiles2 <- grep("R2.fa", fastqfiles, value = TRUE)
    fastqfiles2 <- sort(fastqfiles2)
    # Make BAM file names and remove path/extention from names
    BAMfilenames <- gsub("RNAseqData/|_R1.fasta.gz|_R1.fa", "", 
                         fastqfiles1, ignore.case = T)
    # Align reads to genome
    if ( prod(BAMfilenames %in% list.files()) == 0) {
      align(index="reference_index",
            readfile1=fastqfiles1,
            readfile2 = fastqfiles2,
            type="rna", 
            output_file = BAMfilenames,
            annot.ext = read.table(gene_predict, header = T),
            isGTF = TRUE)      
    }

    # Get annotated raw counts
    rawcounts <- featureCounts(BAMfilenames, 
                               annot.ext = read.table(gene_predict, header = T),
                               isPairedEnd = TRUE)
    
    write_rds(rawcounts, "rawcounts.RDS")
  }
  
  # For single-end
  else {
    #Make BAM file names and remove path/extention from names
    BAMfilenames <- gsub("RNAseqData/|.fasta.gz|.fasta|.fa|.fa.gz", "", 
                         fastqfiles, ignore.case = T)
    if ( prod(BAMfilenames %in% list.files()) == 0) {
    #Align reads to genome
    align(index="reference_index",
          readfile1=fastqfiles,
          type="rna", 
          output_file = BAMfilenames,
          annot.ext = read.table(gene_predict, header = T),
          isGTF = TRUE)
    }  

    rawcounts <- featureCounts(BAMfilenames, 
                               annot.ext = read.table(gene_predict, header = T),
                               isPairedEnd = FALSE)
    
    write_rds(rawcounts, "rawcounts.RDS")
  }
}

  

EZ$BamSummary <- function() {
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


# Post Alignment Functions

# Order dataframe/matrix by selected columns
#Automatically reorders heatmaps by whichever columns selected. For example, 
#to order dataframe (df) in increasing order by the first 3 columns all you need 
#to do is type hm_order(df, c(1,2,3))

#input dataframe/matrix and selected columns
EZ$hm_order <- function(df, cols = c(1), decr = FALSE){
  if(length(cols) > 1){
    hmz <- apply(df, 2, FUN = function(x) 
    {return( (rowMeans(df) - x)/rowSds(df) )} )
    hmrow <- rowSums(hmz[,cols]) #Columns 1 to 3 are control
    df[order(hmrow, decreasing = decr),]
  }
  else {
    hmz <- apply(df, 2, FUN = function(x) 
    {return( (rowMeans(df) - x)/rowSds(df) )} )
    df[order(hmz[,cols], decreasing = decr),]
  }
}


# Column to Rownames
# Takes the first column and makes it a row.
# rename dataframe rows from selected column
EZ$col2namerow <- function(df, col=1) {
  rownames(df) <- df[,col]
  df <- df[,-col]
  return(df)
}

# Sum transcript counts to gene counts
EZ$SumTrans <- function(counts) {
  merge(counts, read_rds("PlAnnotation.RDS")[,1:2], by.x = 0, by.y = 2) %>%
    group_by(GeneID) %>%
    summarise(across(metadata$samples,sum)) %>% as.data.frame() %>%
    col2namerow()  
}



# Rename rows from Gene IDs to gene symbols
EZ$gene2sym <- function (df,IDinRow = T){
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


# Print List of Dataframes to single csv page
# Useful for reporting GO results of several samples
# Useful for reporting genes grouped by functions
EZ$list2csv <- function(DFlist, File = 'DFlist.csv') {
  pdf(File)
  for (ind in 1:length(DFlist)) {
    write.table(paste0(names(DFlist)[ind]), file = File, 
                sep = ",", col.names = FALSE, append = T, row.names = FALSE)
    write.table( data.frame(allRes[[DFlist]]), file = File,  
                 append= T, sep=',', col.names = NA)
  }
  dev.off()
}


