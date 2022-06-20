
#Nucleotide k-mer feature selection

library(dplyr)

library(reshape2)

library(data.table)


# all dataset is your combine nucleotide kmer counts for your all samples

all_dataset <- "#directory of your combine nucleotide kmer counts file"

# antibiotic file (sample name and MIC results (R,I,S))

antibiotic <- "#directory of your antibiotic file" 

#turn kmer counts to binary results

binary_kmer_mat <- all_dataset[, -1] > 0 
rownames(binary_kmer_mat) <- sub("/","",all_dataset[,1])

# Feature elimination ( keep samples absolute mean difference 0.2 threshold between 3 classes)

ant_r <- antibiotic$sample[antibiotic$antibiotic == "R"]
ant_s <- antibiotic$sample[antibiotic$antibiotic == "S"]
ant_i <- antibiotic$sample[antibiotic$antibiotic == "I"]

ant_r_idx <- which(rownames(binary_kmer_mat) %in% ant_r)
ant_s_idx <- which(rownames(binary_kmer_mat) %in% ant_s)
ant_i_idx <- which(rownames(binary_kmer_mat) %in% ant_i)


filter_column <- function(x, thr = 0.2){
  
  r_prop <- mean(x[ant_r_idx])
  s_prop <- mean(x[ant_s_idx])
  i_prop <- mean(x[ant_i_idx])
  
  
  cond1 <- abs(r_prop - s_prop) >= thr
  cond2 <- abs(r_prop - i_prop) >= thr
  cond3 <- abs(s_prop - i_prop) >= thr
  
  return(cond1|cond2|cond3) 
  
}

# Feature elimination ( keep samples absolute mean difference 0.4 threshold between 2 classes)

filter_column_rs <- function(x, thr = 0.4){
  
  r_prop <- mean(x[ant_r_idx])
  s_prop <- mean(x[ant_s_idx])
  
  
  cond1 <- abs(r_prop - s_prop) >= thr
  
  
  return(cond1) 
  
}

keep <- apply(binary_kmer_mat, 2 , filter_column_rs)


filtered_kmer_total_bin <- all_dataset[, c(TRUE, keep)]
filtered_kmer_total_bin[,1] <- sub("/","",filtered_kmer_total_bin[,1])
filtered_kmer_total_bin[filtered_kmer_total_bin > 0] <- 1
filtered_kmer_total_bin$Sample <- antibiotic[,1]


#antibiotic and nucleotide kmers combined by their name


colnames(antibiotic) <- c("Sample","antibiotic")
combine_all_ant <- inner_join(filtered_kmer_total_bin,antibiotic,by = "Sample")


#exclude Intermediate sample

combine_all_ant_rs <- combine_all_ant[combine_all_ant$antibiotic %in% c("R", "S"), ]



#Boruta feature Selection

library(Boruta)

set.seed(111)

combine_all_ant_rs$antibiotic <- as.factor(combine_all_ant_rs$antibiotic)

boruta_all_ant_rs <- Boruta(antibiotic~., data = combine_all_ant_rs, doTrace = 1)


boruta_kmer_stat <- attStats(boruta_all_ant_rs)



boruta_confirmed <- subset(boruta_kmer_stat, subset = boruta_kmer_stat$decision == "Confirmed")
boruta_tentative <- subset(boruta_kmer_stat, subset = boruta_kmer_stat$decision == "Tentative")


plot(boruta_all_ant_rs, las = 2, cex.axis = 0.7)


#green color-important attributes, 
#yellow boxes-tentative attributes 
#red boxes-unimportant.

plotImpHistory(boruta_all_ant_rs)


combine_final_kmer <- combine_all_ant_rs[,which(colnames(combine_all_ant_rs) %in% boruta_confirmed)]


# SNP feature selection

library(dplyr)

library(reshape2)

library(data.table)

#turn kmers to binary results


# data_snp_all contains all samples SNP locations (I-R-S)

data_snp_all <- "#directory of your combine SNP locations file"

#transpose

data_snp_all_t <- as.data.frame(t(data_snp_all))


data_snp_all_t$files <- gsub(".vcf","",as.character(data_rbind$files))

# antibiotic file (sample name and MIC results (R,I,S))

antibiotic <- "#directory of your antibiotic file" 

#turn kmer counts to binary results

binary_snp_mat <- data_snp_all_t[, -1] > 0 
rownames(binary_snp_mat) <- sub("/","",data_snp_all_t[,1])

# Feature elimination ( keep samples absolute mean difference 0.05 threshold between 3 classes)

ant_r <- antibiotic$sample[antibiotic$antibiotic == "R"]
ant_s <- antibiotic$sample[antibiotic$antibiotic == "S"]
ant_i <- antibiotic$sample[antibiotic$antibiotic == "I"]

ant_r_idx <- which(rownames(binary_snp_mat) %in% ant_r)
ant_s_idx <- which(rownames(binary_snp_mat) %in% ant_s)
ant_i_idx <- which(rownames(binary_snp_mat) %in% ant_i)


filter_column <- function(x, thr = 0.05){
  
  r_prop <- mean(x[ant_r_idx])
  s_prop <- mean(x[ant_s_idx])
  i_prop <- mean(x[ant_i_idx])
  
  
  cond1 <- abs(r_prop - s_prop) >= thr
  cond2 <- abs(r_prop - i_prop) >= thr
  cond3 <- abs(s_prop - i_prop) >= thr
  
  return(cond1|cond2|cond3) 
  
}

# Feature elimination ( keep samples absolute mean difference 0.05 threshold between 2 classes)

filter_column_rs <- function(x, thr = 0.05){
  
  r_prop <- mean(x[ant_r_idx])
  s_prop <- mean(x[ant_s_idx])
  
  
  cond1 <- abs(r_prop - s_prop) >= thr
  
  
  return(cond1) 
  
}

keep <- apply(binary_snp_mat, 2 , filter_column_rs)

filtered_snp_total <- data_snp_all_t[, c(TRUE, keep)]

filtered_snp_total$sample <- antibiotic[,1]


#antibiotic and SNPs combined by their name

combine_all_ant <- inner_join(filtered_snp_total,antibiotic,by = "sample")


#exclude Intermediate sample

combine_all_ant_rs <- combine_all_ant[combine_all_ant$antibiotic %in% c("R", "S"), ]


#Boruta feature Selection

library(Boruta)


combine_all_ant_rs$antibiotic <- as.factor(combine_all_ant_rs$antibiotic)

boruta_all_ant_rs <- Boruta(antibiotic~., data = combine_all_ant_rs, doTrace = 1)


boruta_kmer_stat <- attStats(boruta_all_ant_rs)



boruta_confirmed <- subset(boruta_kmer_stat, subset = boruta_kmer_stat$decision == "Confirmed")
boruta_tentative <- subset(boruta_kmer_stat, subset = boruta_kmer_stat$decision == "Tentative")


plot(boruta_all_ant_rs, las = 2, cex.axis = 0.7)


#green color-important attributes, 
#yellow boxes-tentative attributes 
#red boxes-unimportant.

plotImpHistory(boruta_all_ant_rs)

combine_final_snp <- combine_all_ant_rs[,which(colnames(combine_all_ant_rs) %in% boruta_confirmed)]


#AA k-mer feature selection

# all dataset is your combine aa kmer counts for your all samples

all_dataset <- "#directory of your combine aa kmer counts file"

# antibiotic file (sample name and MIC results (R,I,S))

antibiotic <- "#directory of your antibiotic file" 


#antibiotic and aa kmers combined by their name

combine_all_ant <- inner_join(all_dataset,antibiotic,by = "Sample")


#exclude Intermediate sample

combine_all_ant_rs <- combine_all_ant[combine_all_ant$antibiotic %in% c("R", "S"), ]


#Boruta feature Selection

library(Boruta)

set.seed(111)

combine_all_ant_rs$antibiotic <- as.factor(combine_all_ant_rs$antibiotic)

boruta_all_ant_rs <- Boruta(antibiotic~., data = combine_all_ant_rs, doTrace = 1)


boruta_kmer_stat <- attStats(boruta_all_ant_rs)


boruta_confirmed <- subset(boruta_kmer_stat, subset = boruta_kmer_stat$decision == "Confirmed")
boruta_tentative <- subset(boruta_kmer_stat, subset = boruta_kmer_stat$decision == "Tentative")


plot(boruta_all_ant_rs, las = 2, cex.axis = 0.7)


#green color-important attributes, 
#yellow boxes-tentative attributes 
#red boxes-unimportant.

plotImpHistory(boruta_all_ant_rs)


combine_final_aa_kmer <- combine_all_ant_rs[,which(colnames(combine_all_ant_rs) %in% boruta_confirmed)]



