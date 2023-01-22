packages = c("BiocManager","tidyverse", "ggplot2", "dplyr", "RColorBrewer", 
             "ggthemes", "formattable", "gplots", "readr",  "RColorBrewer", 
             "limma", "edgeR", "Rsubread", "topGO", "circlize", "genefilter",
             "RSelenium")



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
EZ$DLgenome <- function() 
  {
  gitlink <- "http://planmine.mpibpc.mpg.de/planmine/model/bulkdata/dd_Smes_g4.fasta.zip"
  getOption('timeout')
  options(timeout=1500)
  download.file(gitlink, "genome.fa")
  }

source("")

# Fetch Gene Predictions
EZ$DLgene_predict <- function() 
{
  gitlink <- "https://raw.githubusercontent.com/HedonicHotspot/EZworm/main/smes_v2_repeatfil_YAI.saf"
  as.data.frame(read_table(gitlink))
}

# For aligning RNAseq data
EZ$index <- function( ) 
  {
  if ("genome.fa" %in% list.files() == F){DLgenome()}
  genome <- paste0(getwd(), '/genome.fa')
  if ("smes_v2_repeatfil_YAI.saf" %in% list.files() == F) {GETgene_predict()}
  gene_predict <- paste0(getwd(),'/smes_v2_repeatfil_YAI.saf') 
  RNAseqData <- paste0(getwd(), "/RNAseqData/")
  fastqfiles <- list.files(path = RNAseqData, full.names = TRUE)
  fastqfiles <- grep(".fa", fastqfiles, value = TRUE)
  buildindex(basename="reference_index", reference= genome)
  }

EZ$align <- function() 
{
  # For paired-end
  if (sum(grepl("_R1.fa|_R2.fa", fastqfiles)) > 2) 
  {
    #Each pair read file in seperate vectors
    fastqfiles1 <- grep("_R1.fa", fastqfiles, value = TRUE)
    fastqfiles1 <- sort(fastqfiles1)
    fastqfiles2 <- grep("_R2.fa", fastqfiles, value = TRUE)
    fastqfiles2 <- sort(fastqfiles2)
    # Make BAM file names and remove path/extention from names
    BAMfilenames <- substr(fastqfiles1, nchar(RNAseqData), nchar(fastqfiles)-6)
    # Align reads to genome
    align(index="reference_index",
          readfile1=fastqfiles1,
          readfile2 = fastqfiles2,
          type="rna", 
          output_file = BAMfilenames,
          annot.ext = gene_predict,
          isGTF = TRUE)
    # Get annotated raw counts
    rawcounts <- featureCounts(bam.files, 
                               annot.ext = gene_predict,
                               isPairedEnd = TRUE)
    
    write_delim(as.data.frame(rawcounts$counts), paste(getwd(), "/rawcounts.txt", sep = ""))
  }
  # For single end
  else 
  {
    #Make BAM file names and remove path/extention from names
    BAMfilenames <- substr(fastqfiles, nchar(RNAseqData), nchar(fastqfiles)-3)
    #Align reads to genome
    align(index="reference_index",
          readfile1=fastqfiles,
          type="rna", 
          output_file = BAMfilenames,
          annot.ext = gene_predict,
          isGTF = TRUE)
    bam.files <- paste0(getwd(), "/RNAseqAlignment/", BAMfilenames)
    bam.path <- paste0(getwd(), "/RNAseqAlignment/")
    rawcounts <- featureCounts(bam.files, 
                               annot.ext = gene_predict,
                               isPairedEnd = FALSE)
    
    write_delim(as.data.frame(rawcounts$counts), 
                paste0(getwd(), "/rawcounts.txt"))
  }
}



  

EZ$BamSummary <- function() 
  {
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
#Takes the first column and makes it a row.

# rename dataframe rows from selected column
EZ$col2namerow <- function(df, col=1) {
  rownames(df) <- df[,col]
  df <- df[,-col]
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

