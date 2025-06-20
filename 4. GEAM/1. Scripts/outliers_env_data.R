##################################
## ADDING ENVIRONMENTAL DATA ##
##################################

# to join them with environmental data, in which population names appear in first column, 
# transpose genetic data tables and 
# establish population names in the first column 
rownames(frh1_final) <- frh1_final[,1]; frh1_final <- frh1_final[, -1] # genetic positions (SNPs) as rownames to make the transposition
gen1 <- as.data.frame(t(frh1_final)) # transpose 
gen1 <- as.data.frame(cbind("pop" = rownames(gen1), gen1)) # as first column population names 
rownames(gen1) <- NULL # remove rownames 

# open environmental data 
setwd("")
# gen.env = environment of populations in 1000genomes database 
load("3.environment_1000genomes_populations.RData")
# sim.env = environment of populations in Simons database
load("3.environment_simons_populations.RData") 

# change 3 letter codes by complete population names in genomic 
# tables to be able to join to the environmental ones by merge function 
setwd("")
corr <- read.table("correspondence_names_populations.txt", header=T, sep="\t")

gen1 <- merge(gen1, corr, by.x="pop", by.y="Code")
gen1 <- as.data.frame(cbind(gen1$pop.y, gen1[, 2 : 259])) # column with complete pop names the first and join with genetic position data 
colnames(gen1)[1] <- "pop" # rename the column with complete names 

# check if all populations are in both databases 
all(gen1$pop %in% gen.env$pop)
all(gen1$pop %in% sim.env$pop)

# not all populations are in Simons database
# use 1000genomes database
gen1 <- merge(gen1, gen.env, by="pop")

# save complete databases with significant SNPs and environmental data 
write.csv(gen1, "./gen1.csv", row.names=FALSE)

# there are too much variables to include in a PCA 
# make 1000 random subsamples of 5% of the total of SNPs
(5 * (ncol(gen1)-6)) / 100 # 5% of the gen1 SNPs are aprox. 13 SNPs by table

# create output empty lists
gs1 <- list()
for(i in 1 : 1000){
  
  # remove from tables all that do not correspond to genetic positions 
  # and make the subsample 
  gs1[[i]] <- gen1[, -c(1, 260 : 264)][, sample(1 : ncol(gen1[, -c(1, 260 : 264)]), 13)] 
  
  # add again environmental variables 
  gs1[[i]] <- as.data.frame(cbind(gen1[, 262:264], gs1[[i]]))
  
  # add X to the name of genetic position (easier work with them as a character
  # and not as an integer)
  colnames(gs1[[i]])[-c(1 : 3)] <- paste("X", colnames(gs1[[i]])[-c(1 : 3)], sep="")
  
  # save each subsample as an element of the output lists
  gs1[[i]] <- as.data.frame(apply(gs1[[i]], 2, as.numeric))
  
}

# save random tables (input of GEAM)
setwd("")
save(gs1, file = "1000_random_SNPs_gen1.RData")
