# charge functions 
source("./extract_outliers_functions.R")

###################################
## 1. FILTERING SIGNIFICANT SNPs ##
###################################

# open GRoSS results 
setwd("")

h1 <- read.table("human.tsv", header=T, sep="\t")

# min p-value for each genetic position (SNP)
# create an empty vector 
mins <- numeric() 
for(i in 1:nrow(h1)){
  # save the min p-value for a given position (SNP)
  mins[i] <- min(h1[i,38:71]) 
  print(i)
}

# add min p-value to the previous tables only with columns one ("CHR"=nº chromosome) and two ("START"=SNP position)
h1f <- as.data.frame(cbind(h1[,1], h1[,2], mins))
colnames(h1f) <- c("CHR", "START", "min_Pvalue")

# plot histogram of p-values for Scree plot
x11()
hist(h1f$min_Pvalue, breaks=5000, xlim=c(0,0.01))

# number of SNPs with p-value < 10E-4
sum(h1f$min_Pvalue<0.00001)

# number of SNPs with p-value < 10E-5
sum(h1f$min_Pvalue<0.000001)

# P-value=10E-5 was selected as cut-off for giving the lowest number of variables (significant SNPs)

# save only significant SNPs (with P-value < 10E-5
h1s <- h1[which(mins < 0.00001), ]

# filter allele frequency files by significant SNPs ()
# enter the directory where allele frequency files (.afreq) are
# (each file corresponds to a population --> 18 in total)
setwd("./Allele_freq_by_pop")

# create output empty lists. Each element of the lists will correspond to the filtering  
# of allele frequency tables of each population (fr) by significant SNPs (those in h1s and h2s)
frh1_list <- list()

for(i in 1:length(list.files(pattern="afreq"))){
  
  fr <- read.table(list.files(pattern="afreq")[[i]])
  
  # save results as elements of output lists (=one population)
  frh1_list[[i]] <- filtering(fr, h1s)
  
  print(i)
}

# convert previous results from list to table
# first column = genetic position (SNP) (the same in all the elements of the lists)
frh1_final <- frh1_list[[1]][,1]
for(i in 1:length(frh1_list)){
  
  # paste the second column (allelic frequencies) of each element of the list (each population) 
  # to the first column (genetic position)
  frh1_final <- as.data.frame(cbind(frh1_final, frh1_list[[i]][, 2]))
  colnames(frh1_final)[i + 1] <- colnames(frh1_list[[i]])[2]
}

#write.table(rownames(frh1_final), "SNPs_under_sel.txt", quote = F, col.names = F, row.names = F)
