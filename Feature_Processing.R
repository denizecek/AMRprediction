

#Nucleotide k-mer count from contig.fasta files 

library(ape)

library(kmer)


kmer_input <- "#directory of your fasta file"

all_kmer <- list.files(kmer_input, full.names = TRUE)

for (sample_file in all_kmer) {
  
  
  cat("Working on:", which(all_kmer == sample_file), "out of", length(all_kmer), "       \r")
  contig <- read.dna(sample_file,format = "fasta", as.matrix = NULL)
  x=kcount(contig, k=10)
  write.table(x,sub("fasta","RDS", sample_file), append = FALSE, sep = " ", dec = ".",row.names = TRUE, col.names = TRUE)
}

# Combine all nucleatide k-mer files in one file

kmer_output <- "#directory of your kmer files"

all_samples <- list.files(kmer_output, full.names = TRUE)

all_tots_per_kmer <- data.frame()


for (sample_file in all_samples) {
  sample_name <- sub("_contigs\\.txt", "", sub("all_kmer", "", sample_file))
  
  cat("Working on:", sample_name, which(all_samples == sample_file), "out of", length(all_samples), "       \r")
  
  tmp <- data.table::fread(sample_file)
  tmp <- as.data.frame(tmp)
  tmp$V1 <- NULL
  
  total_per_kmer <- apply(tmp, 2, sum)
  
  all_tots_per_kmer <- rbind(all_tots_per_kmer,
                             data.frame(Sample = sample_name,
                                        t(total_per_kmer)))
  
  if (!all(names(total_per_kmer) == colnames(all_tots_per_kmer)[-1])) {
    stop("kmers do not match!")
  }
  
}

# SNP information from vcf files.

library(dplyr)

library(reshape2)

library(data.table)

library(tidyverse)

all_vcf <- "#directory of your vcf files"

files <- list.files(all_vcf ,full.names = TRUE)

for (sample_file in files) {
  cat("Working on:", which(files == sample_file), "out of", length(files), "       \r")
  new <- read.table(sample_file,header = TRUE)
  x   <- data.frame(new[,2])
  write.table(x,sub("RDS","txt",sample_file), append = FALSE, sep = " ", dec = ".",row.names = TRUE, col.names = TRUE)
}

all_snp <- "#directory of your snp files"

files <- list.files(all_snp ,full.names = TRUE)

data <- lapply(files, read.table, header=TRUE, sep=" ")

value <- "1"

for (i in 1:length(data)){data[[i]]<-cbind(data[[i]],files[i],value[i])}
data_rbind <- do.call("rbind", data)

data_snp <- acast(data_rbind, pos ~ files, value.var = "value" )

all_snp <- list.files("all_kmer", full.names = TRUE)

# Combine all SNP locations in one file

for (sample_file in all_snp) {
  
  
  cat("Working on:", which(all_snp == sample_file), "out of", length(all_snp), "       \r")
  new <- read.table(sample_file,quote="\"",header = FALSE)
  x=subset(data2, V3)
  write.table(x,sub("txt","RDS",sample_file), append = FALSE, sep = " ", dec = ".",row.names = TRUE, col.names = TRUE)
}

#AA k-mer count from protein.fasta files


library(seqinr)

kmer_aa_input <- "#directory of your protein.fasta files"

all_aa_kmer <- list.files(kmer_aa_input, full.names = TRUE)

for (sample_file in all_aa_kmer) {
  
  
  cat("Working on:", which(all_aa_kmer == sample_file), "out of", length(all_aa_kmer), "       \r")
  contig <- read.FASTA(sample_file, type = "AA" )
  x=kcount(contig, k = 5, residues = NULL, gap = "-", named = TRUE,
           compress = TRUE, encode = FALSE)
  new_x <- data.frame(colSums(x))
  new_x <- t(new_x)
  
  write.table(new_x,sub("fasta","txt", sample_file), append = FALSE, sep = " ", dec = ".",row.names = TRUE, col.names = TRUE)
}


# Combine all aa k-mer files in one file

kmer_aa_output <- "#directory of your aa k-mer files"

all_samples <- list.files(kmer_aa_output, full.names = TRUE)

all_tots_per_kmer <- data.frame()


for (sample_file in all_samples) {
  sample_name <- sub("\\.txt", "", sub(kmer_aa_output, "", sample_file))
  
  cat("Working on:", sample_name, which(all_samples == sample_file), "out of", length(all_samples), "       \r")
  
  tmp <- data.table::fread(sample_file)
  tmp <- as.data.frame(tmp)
  tmp$V1 <- NULL
  
  total_per_kmer <- apply(tmp, 2, sum)
  
  all_tots_per_kmer <- rbind(all_tots_per_kmer,
                             data.frame(Sample = sample_name,
                                        t(total_per_kmer)))
  
  if (!all(names(total_per_kmer) == colnames(all_tots_per_kmer)[-1])) {
    stop("kmers do not match!")
  }
  
}



