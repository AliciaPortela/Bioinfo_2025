#### CHARGE FUNCTIONS AND LIBRARIES ####
library(sm) # for density plots 
library(RColorBrewer) # color package
library(MASS) # statistic library with lm function for linear regression 
library(raster) # for spacial (raster) data

# my functions 
source("C:/Users/Alicia/Desktop/AntropoGeo/GEAM_inference_functions.R")

# set working directory 
setwd("C:/Users/Alicia/Desktop/AntropoGeo")

# charge input tables with 1000 random samples 
load("1000_random_SNPs_gen1.RData")

# 1. Find the best model
# implement a multiple linear regression for each random sample (1000)
# save AICs and Adj_Rsq of each of the 1000 models 

res <- as.data.frame(matrix(ncol = 3, nrow = 0))
for(i in 1 : length(gs1)){
  
  outp <- GEAM(gs = gs1[[i]], gen = "G1EAM", env = "rad")
  res <- as.data.frame(rbind(res, cbind("id" = i, outp)))
  
  outp <- GEAM(gs = gs1[[i]], gen = "G1EAM", env = "th")
  res <- as.data.frame(rbind(res, cbind("id" = i, outp)))
  
  outp <- GEAM(gs = gs1[[i]], gen = "G1EAM", env = "ele")
  res <- as.data.frame(rbind(res, cbind("id" = i, outp)))
  
  print(i)
}

# save table with AIC and Adj_Rsq for each model 
setwd("Master-Thesis-MBC/7. GEAM and distribution range inference/Results/1. Searching the best model")
write.table(res, "AIC_adjRsq.txt", row.names=FALSE, quote=F, sep="\t")

# save plots of AIC and Adj_Rsq distributions
plots_save(Label="G1EAM_rad", var="adj_Rsq")
plots_save(Label="G1EAM_th", var="adj_Rsq")
plots_save(Label="G1EAM_ele", var="adj_Rsq")
plots_save(Label="G1EAM_rad", var="AIC")
plots_save(Label="G1EAM_th", var="AIC")
plots_save(Label="G1EAM_ele", var="AIC")



# 2. Obtaining the best model
# select the regression model with the lowest value of AIC for each environmental variable and genomic database (6 models)

setwd("Master-Thesis-MBC/7. GEAM and distribution range inference/Results/2. Best model/")
bg_results <- best_GEAM()

# save results of the obtained models in table format 
write.table(bg_results$table, "GEAM_results.txt", row.names = F, quote = F, sep = "\t")
# save coefficient tables in a variable to work with them in the range inference 
coeffs <- bg_results$Coeffs


# 3. Inference of distribution range applying best regression models  
# Define the environmental limits of the distribution range
# and fill the cells of the raster that are between the minimum and maximum values

# create a table to store the values of the extremes of the distribution range
extremes <- as.data.frame(matrix(ncol=3, nrow=3))
colnames(extremes) <- c("Label", "Zeros", "Ones")
extremes$Label <- bg_results$table$Label # fill the first column with the gs-env.var id (e.g. "G1EAM_rad")

# loop to obtain extremes in each model (6)
for(i in 1:nrow(extremes)){
  
  # identify coefficients table for that model
  coeffi <- coeffs[[which(names(coeffs) == extremes$Label[i])]]
  
  # ZEROS = First extreme = allele frequency of 0 in all SNPs: 
  # Environmental value = Intercept + 0*SNP1 + 0*SNP2 + 0*SNP3 ... + 0*SNPn. = Intercept
  # Example for G1EAM ~ Radiation
  # 5.0909185 + 0*(-10.9351638) +  0*8.8152773  + 0*(-2.4859137) +
  # 0 *(-4.7480158) + 0*3.0492881 + 0*6.7565186 + 0*(-2.7323701) +
  # 0*(-10.0658544) + 0*0.534166 = 5.0909185 = Intercept 
  extremes[extremes$Label == extremes$Label[i], "Zeros"] <- coeffi[1, 1] # pass the intercept value to the table 
  
  # ONES = Second extrem = allele frequency of 1 in all SNPs: 
  # Environmental value = Intercept + 1*SNP1 + 1*SNP2 + 1*SNP3 ... + 1*SNPn = Intercept + SNP1 + SNP2 + SNP3 ... + SNPn
  # Example for G1EAM ~ Radiation:
  # 5.0909185 + 1*(-10.9351638) +  1*8.8152773  + 1*(-2.4859137) +
  # 1*(-4.7480158) + 1*3.0492881 + 1*6.7565186 + 1*(-2.7323701) +
  # 1*(-10.0658544) + 1*0.534166 = Sum of the Estimate column.
  extremes[extremes$Label == extremes$Label[i], "Ones"] <- sum(coeffi[, 1]) # pass the sum of Estimate column to the table 
  
}

# open rasters of environmental variables to plot the distribution range 
setwd("Master-Thesis-MBC/1. Environmental data retrieval/Environmental rasters")
rrad <- raster("rad.tif")
rth <- raster("th.tif")
rele <- raster("ele.tif")

# working directory to save range distribution inference plot results 
setwd("Master-Thesis-MBC/7. GEAM and distribution range inference/Results/2. Best model")

# function to plot distribution range
# needs to be debugged
GEAM_extremes_plot <- function(ras, label, env.label){
  
  # comparing extreme values of rasters with 
  # the predicted by GEAM and then decide to plot or not:
  
  # 1. if they remain inside the raster margins (at least one)...
  # plot
  ext_pred <- extremes[extremes$Label == label, 2:3]
  
  min_ras <- min(values(ras)[-which(is.na(values(ras)))])
  max_ras <- max(values(ras)[-which(is.na(values(ras)))])
  
  if(min(ext_pred) > min_ras | max(ext_pred) < max_ras){
    
    # change to '1000' NAs (NA=sea) to operate 
    # with them
    values(ras)[which(is.na(values(ras)))] <- 1000
    
    values(ras) <- round(values(ras), digits=6)
    presences <- which(values(ras) >= min(ext_pred) & values(ras) <= max(ext_pred))
    
    # if an environmental value is predicted in the sea, remove it
    if(any(presences %in% which(values(ras) == 1000))){
      presences <- presences[-which(presences %in% which(values(ras)==0))]
    }
    values(ras)[presences] <- 950
    values(ras)[! values(ras) %in% c(1000, 950)] <- 500
    
    cuts = c(0,600,960,1100) # set breaks for colors
    pal <- colorRampPalette(c("gray","darkred", "white"))
    
    # make and save maps
    x11()
    tiff(paste(label, ".tiff", sep=""), units="in", width=15, height=9, res=300)
    plot(ras, breaks=cuts, col = pal(3), legend=F, main=label, 
         xlab=paste("from", round(ext_pred[1,1], digits=2),
                    "to", round(ext_pred[1,2], digits=2)))
    dev.off()
    
  }else{ 
    
    # 2. if they remain outside the raster margins, (=prediction occupies all the environment),
    # the map does not make sense (=all). 
    # Do not plot it
    # Message:
    print("Prediction occupies the whole environment")
  }
  
}

GEAM_extremes_plot(ras=rth, label="G1EAM_th")
GEAM_extremes_plot(ras=rrad, label="G1EAM_rad")
GEAM_extremes_plot(ras=rele, label="G1EAM_ele")