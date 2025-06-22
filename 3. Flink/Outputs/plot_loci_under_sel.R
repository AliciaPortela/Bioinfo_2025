setwd("C:/Users/Alicia/Desktop/AntropoGeo/")


args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
rm(args)

pdf(paste("plotSel", name, ".pdf", sep=""), width = 14, height = 3.6 * 8)  # más alto
layout(matrix(rep(1:8, each=3), ncol=1, byrow=TRUE))  # 8 paneles: 1 jerarquía + 7 grupos
par(mar=c(1,6.5,1.0,2.0), mgp=c(3.5,1.5,0), oma=c(6,0,0,0))

library(data.table)

# Leer posiciones de loci
freq_pos <- read.table("input.flink.txt", header=FALSE, skip=2)
numloci <- nrow(freq_pos)

# Funciones de selección
calcDivSel <- function(data, numloci) {
  rows <- nrow(data)
  (colSums(data[,2:(numloci+1)] > 0)) / (rows + 1)
}

calcBalSel <- function(data, numloci) {
  rows <- nrow(data)
  (colSums(data[,2:(numloci+1)] < 0)) / (rows + 1)
}

# Función para plotear
plotdata <- function(freq, divPostProb, balPostProb, namepop) {
  plot(-log10(1.0 - divPostProb), main="", ylab="", xlab="", type='l', ylim=c(0,4),
       las=1, cex.lab=2.5, cex.axis=2, col="orange2", xaxs="i", yaxs="i", xaxt="n", yaxt="n")
  lines(-log10(1.0 - balPostProb), col="dodgerblue")
  legend("topright", legend=namepop, bty="n", pch=NA, cex=1.6, inset=c(-0.0001,-0.15))
  mtext(expression(paste(-log[10]~q)), side=2, line=4, cex=1.3)
  axis(2, seq(0, 4, 2), las=2, cex.axis=2)
  abline(h=2.0, lty="dashed")
}

# 1. Jerarquía superior
data_flink_g <- fread("Posterior_A.txt", header=TRUE)
divPostProb <- calcDivSel(data_flink_g, numloci)
balPostProb <- calcBalSel(data_flink_g, numloci)
plotdata(freq_pos$V2, divPostProb, balPostProb, "Higher hierarchy")

# 2. Grupos (0 a 6)
for (i in 0:6) {
  file <- paste("Posterior_alphas_group_", i, ".txt", sep="")
  data_flink_g <- fread(file, header=TRUE)
  divPostProb <- calcDivSel(data_flink_g, numloci)
  balPostProb <- calcBalSel(data_flink_g, numloci)
  plotdata(freq_pos$V2, divPostProb, balPostProb, paste("Group", i))
}

# Eje X en el último panel
axis(1, las=TRUE, cex.axis=2.0)
mtext("Locus", side=1, line=4, cex=1.6, outer=TRUE)

dev.off()