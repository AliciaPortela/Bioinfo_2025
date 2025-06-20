# funtion to filter table f by genetic position of table h
# (in this case: f=fr, h=hs1/hs2)

filtering <- function(f, h){
  # V2 is the SNP position. 1: is the chromosome number, it is always 1 because 
  # we are working with only one chromosome
  # remove 1:
  f$V2 <- as.numeric(gsub("1:", "", f$V2))
  
  # join fr with h1s only by the genetic position with merge 
  # keeping only those that are in both tables = significant SNPs 
  fh <- merge(h, f, by.x="START", by.y="V2")
  
  # keep only in the table genetic position (col.1) and 
  # allelic frequencies (col.75)
  fh <- fh[,c(1, 75)]
  
  # rename allele frequency col. with the population name 
  # which is in the name of the .afreq file 
  # remove everything that is not the name of the population and replace
  # by nothing 
  colnames(fh)[2] <- gsub(".afreq", "", gsub("allele-freq_chr12_all_SNPs.","",list.files(pattern="afreq")[[i]]))
  
  return(fh)
}
