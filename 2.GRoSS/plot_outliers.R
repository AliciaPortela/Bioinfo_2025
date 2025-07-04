library("ggplot2")
library("reshape")

#popnames <- "1KG_MSL_ESN_CDX_JPT_CEU_TSI_CHB"
setwd("")

table <- read.table("human.tsv",header=TRUE)
table <- as.data.frame(table)


# WINDOWS

#average <- c()
#for(chr in seq(1,22)){
#print(chr)
#onechr <- table[which(table[,1]==chr),]
#pvalpos <- which(mapply(grep,"Pval",colnames(onechr)) == 1)
#chipos <-  which(colnames(onechr) %in% colnames(onechr)[-c(pvalpos)][-c(1,2,3)])
#onechr <- as.data.frame(t(apply(onechr,1,function(row){row[pvalpos] <- -log10(row[pvalpos]); return(row)})))
#winsize <- 10
#winstep <- 1
#startvec <- seq(1,dim(onechr)[1]-winsize+1,winstep)
#endvec <- seq(winsize,dim(onechr)[1],winstep)
#winlims <- cbind(startvec,endvec)


# AVERAGE CHI-SQUARED (10 SNPs)
#averagetemp <- as.data.frame(t(apply(winlims,1,function(window){     
#block <- onechr[seq(window[1],window[2]),]
#chisqstats <- block[,chipos]
#avgstats <- apply(chisqstats,2,mean)
#avgpval <- 1 - pchisq(avgstats,1)
#return(c(chr,onechr[window[1],2],onechr[window[2],3],avgstats,avgpval))
#})))
#colnames(averagetemp) <- colnames(onechr)
#average <- rbind(average,averagetemp)
#}


# PLOT

finaltab <- table
#finaltab <- average

CHR <- finaltab[,1]
START <- finaltab[,2]
END <- finaltab[,3]
MIDPOINT <- (finaltab[,2]+finaltab[,3])/2
newtab <- finaltab[,seq(4,dim(finaltab)[2])]
newtab <- newtab[,seq(dim(newtab)[2]/2+1,dim(newtab)[2])]
POS <- seq(1,dim(newtab)[1])
melttab <- melt(newtab)
melttab[,2] <- -log10(melttab[,2])
melttab <- cbind(melttab, POS,CHR,START,END)
melttab[,1] <- factor(melttab[,1])
names(melttab) <- c("Branches","Pvals","SNPID","CHR","START","END")
totalbranches <- length(unique(melttab[,1]))

plot1 <- ggplot(melttab, aes(x=SNPID, y=Pvals, fill=Branches, shape=Branches,colour=Branches))+
geom_point()+
scale_shape_manual(values=rep(c(4,15,17,18,19),100)[seq(1,totalbranches)])
print(plot1)
ggsave("outliers.png", plot = plot1, width = 12, height = 10, dpi = 300)
