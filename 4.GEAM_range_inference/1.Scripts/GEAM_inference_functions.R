# 1. "GEAM" function implements multiple linear regression (lm)
# for each gs element (each 1000 random samples in list); genomic database and environmental variable
# giving as output only the preliminar values of AIC (Akaike Information Criterion) 
# and Adj_Rsq of the best models (1000)
GEAM <- function(gs, gen, env){
  
  # remove environmental variables that are not being used 
  if(env == "rad"){
    gs <- gs[, -c(2, 3)]
  }
  if(env == "th"){
    gs <- gs[, -c(1, 3)]
  }
  if(env == "ele"){
    gs <- gs[, -c(1, 2)]
  }
  
  colnames(gs)[1] <- "var"
  
  # create multiple linear regression model 
  # with environmental variable as response/dependent variable 
  # and all SNPs as predictive/independent variables 
  full <- lm(formula = var ~ ., data = gs)
  
  # if any of the coefficients = NA, remove it 
  while(any(is.na(full$coefficients))){
    remove <- names(which(is.na(full$coefficients)))
    remove <- which(colnames(gs) %in% remove)
    # repeat the model without it 
    full <- lm(formula = var ~ ., data = gs[, -c(remove)])
  }
  # find the best model with lowest AIC: stepwise regression
  # starting point: previous model (full)
  # searching in both directions (forward and backward): adding and removing variables 
  step <- stepAIC(full, direction = "both", trace = FALSE)
  
  # save results (AIC and adjusted Rsq) in table format 
  outs <- as.data.frame(matrix(ncol = 3, nrow = 1))
  colnames(outs) <- c("Label", "AIC", "adj_Rsq")
  outs[1, "Label"] <- paste(gen, "_", env, sep = "")
  outs[1, "AIC"] <- extractAIC(step)[2]
  outs[1, "adj_Rsq"] <- summary(step)[9]
  
  return(outs)
}

# 2. "plots_save" function to save distributions of AIC and Adj_Rsq of each of the 1000 models 
# in the exploratory searching of the best among the 1000 random samples 
plots_save <- function(Label, var){
  
  titu <- paste(var, gsub("_", "~", Label))
  
  if(var == "adj_Rsq"){
    jpeg(paste(var, Label, ".jpg", sep="_"), units="cm", width=15, height=15, res=300)
    plot(density(res[res$Label==Label, "adj_Rsq"]), main=titu, lwd=3, lty=c(1,1,1,1))
    dev.off()
  }
  if(var == "AIC"){
    jpeg(paste(var, Label, ".jpg", sep="_"), units="cm", width=15, height=15, res=300)
    plot(sort(res[res$Label==Label, "AIC"]), pch=20, ylab = titu)
    dev.off()
  }
}

# 3. "lmp" function to obtain the p-value in multiple linear regression models
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

# 4. "best_GEAM_in" function finds the random sampling of the SNPs 
# with which the best multiple linear regression model was obtained (lowest AIC) for a given 
# genomic database (gs) and environmental variable (rad, th, ele)
# and saves information of that model (AIC, p-value, adj_Rsq, non-sig SNPs and coefficients)
best_GEAM_in <- function(res, gs, label){
  
  # 1. Get input table with the best regression (lowest AIC)
  # for the correspondent gs and environmental variable (label)
  best <- min(res[res$Label == label, "AIC"]) # min AIC 
  table <- gs[[res[which(res$AIC == best), "id"]]] # obtaining id of the best model (with lowest AIC) from res table which corresponds to the index of gs list 
  # (where the 1000 subsamples (input tables) are)

  # 2. GEAM
  ## 2.1. Running again initial GEAM (full) for the obtained input table 
  ## in order to identify non-significant SNPs for the three environmental variables:
  
  ### 2.1.1. Radiation
  if(label %in% c("G1EAM_rad")){
    
    table <- table[,-which(colnames(table) %in% c("th", "ele"))] # remove th and ele environmental variables
    full0 <- lm(formula = rad ~ ., data = table) # create the multiple linear regression model for that input table 
    step0 <- stepAIC(full0, direction = "both", trace = FALSE) # stepwise regression for the previous model 
    names_nosig <- names(which(summary(step0)$coefficients[,4] >= 0.05)) # get the names of non-significant SNPs
    
    # if there is any non-significant SNP, remove it and run GEAM again without it/them 
    if(length(names_nosig > 0)){
      table <- table[,-which(colnames(table) %in% names_nosig)] # remove non-significant SNPs from input table 
      full <- lm(formula = rad ~ ., data = table) # # run lm again without non-significant SNPs
    }else{
      full <- full0 # change the name of the model in case there is not any non-significant SNP
    }
  }
  
  ### 2.1.2. Temperature & Humidity
  #if(label %in% c("G1EAM_th")){
    
    #table <- table[,-which(colnames(table) %in% c("rad", "ele"))] 
    #full0 <- lm(formula = th ~ ., data = table) 
    #step0 <- stepAIC(full0, direction = "both", trace = FALSE) 
    #names_nosig <- names(which(summary(step0)$coefficients[,4] >= 0.05)) 
    
    #if(length(names_nosig > 0)){ 
     # table <- table[,-which(colnames(table) %in% names_nosig)]
    #  full <- lm(formula = th ~ ., data = table)
    #}else{
     # full <- full0
    #}
  #}
  
  ### 2.1.3. Elevation
  #if(label %in% c("G1EAM_ele")){
    
   # table <- table[,-which(colnames(table) %in% c("th", "rad"))]
    #full0 <- lm(formula = ele ~ ., data = table)
    #step0 <- stepAIC(full0, direction = "both", trace = FALSE)
    #names_nosig <- names(which(summary(step0)$coefficients[,4] >= 0.05))
    #if(length(names_nosig > 0)){
     # table <- table[,-which(colnames(table) %in% names_nosig)]
      #full <- lm(formula = ele ~ ., data = table)
    #}else{
     # full <- full0
    #}
  #}
  
  ## 2.2. Run Stepwise Regression for the previous model (full) 
  step <- stepAIC(full, direction = "both", trace = FALSE) 
  
  # 3. Save results as output of the function 
  outs <- list("best_AIC" = 0, "p_value" = 0, "adj_R_sq" = 0, "x" = 0, "Coeffs." = 0) # create a list with 5 elements to save information of interest 
  outs$best_AIC <- best # save AIC
  outs$p_value <- lmp(step) # save p-value 
  outs$adj_R_sq <- unlist(summary(step)[9]) # convert list in vector: adj_Rsq is in list inside a list/save adj_Rsq
  names(outs[[2]]) <- NULL # make the name of the element NULL (lmp is giving already the name "p-value", avoid double the name)
  outs$x <- names_nosig # save non-sig SNPs
  names(outs)[4] <- "SNPs>0.05" # rename element (outs$SNPs>0.05 gives problems)
  outs$Coeffs. <- summary(step)[[4]] # save coefficients 
  
  ## save coefficients as table
  write.table(outs$Coeffs., paste(label, ".txt", sep=""), quote = F)
  
  return(outs)
}


# 5. "best_GEAM" function applies "best_GEAM_in" function (to find the best model) 
# for each genomic database (gs) and environmental variable (rad, th, ele): 6 times 
# and saves global results in table
best_GEAM <- function(){
  
  # create empty table to save the results 
  results <- as.data.frame(matrix(ncol=5, nrow=1)) 
  colnames(results) <- c("Label", "AIC", "p_value", "adj_R_sq", "SNPs>0.05") # save information except from coefficients 
  # (saved in tables with the "best_GEAM_in" function)
  results$Label <- c("G1EAM_rad")
  
  # create empty list to save the coefficients tables 
  coeffs <- list()
  
  # run GEAM for gs-env.var combination (6 times):
  ## G1EAM ~ Radiation
  geam <- best_GEAM_in(res, gs1, "G1EAM_rad")
  # save information in global results table (except from coefficients)
  results[results$Label=="G1EAM_rad", "AIC"] <- geam$best_AIC 
  results[results$Label=="G1EAM_rad", "p_value"] <- geam$p_value
  results[results$Label=="G1EAM_rad", "adj_R_sq"] <- geam$adj_R_sq
  results[results$Label=="G1EAM_rad", "SNPs>0.05"] <- length(geam[[4]])
  # save coefficients in list 
  coeffs[[1]] <- geam$Coeffs.
  names(coeffs)[1] <- "G1EAM_rad" 
  
  ## G1EAM ~ Temperature & Humidity
 # geam <- best_GEAM_in(res, gs1, "G1EAM_th")
  #results[results$Label=="G1EAM_th", "AIC"] <- geam$best_AIC
  #results[results$Label=="G1EAM_th", "p_value"] <- geam$p_value
  #results[results$Label=="G1EAM_th", "adj_R_sq"] <- geam$adj_R_sq
  #results[results$Label=="G1EAM_th", "SNPs>0.05"] <- length(geam[[4]])
  #coeffs[[2]] <- geam$Coeffs.
  #names(coeffs)[2] <- "G1EAM_th"
  
  ## G1EAM ~ Elevation
  #geam <- best_GEAM_in(res, gs1, "G1EAM_ele")
  #results[results$Label=="G1EAM_ele", "AIC"] <- geam$best_AIC
  #results[results$Label=="G1EAM_ele", "p_value"] <- geam$p_value
  #results[results$Label=="G1EAM_ele", "adj_R_sq"] <- geam$adj_R_sq
  #results[results$Label=="G1EAM_ele", "SNPs>0.05"] <- length(geam[[4]])
  #coeffs[[3]] <- geam$Coeffs.
  #names(coeffs)[3] <- "G1EAM_ele"
  
  out <- list("Coeffs" = coeffs, "table" = results)
  
  return(out)
}
