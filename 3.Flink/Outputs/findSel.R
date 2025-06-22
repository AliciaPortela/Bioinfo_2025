args <- commandArgs(trailingOnly = TRUE);
name <- args[1];  #chromosome or group of scaffolds
group <- args[2];  #give a name for the group (dolphin example: Coastal or Pelagic)
index <- args[3];  #index of the group in the Flink output. It starts from 0
sel <-args[4];  #divergent or balancing
prob <-as.numeric(args[5]); #suggested value 0.99 that corresponds to FDR=0.01
rm(args);

if((sel != "divergent") && (sel != "balancing")) {
  stop("type a proper kind of selection: divergent or balancing!")
}

library(data.table)
freq_pos<-read.table("input.flink.txt",sep="",header=F,skip = 2)
numloci<-nrow(freq_pos)
if(group == "Hierarchy"){
  data_flink<-read.table("Posterior_A.txt",sep="",header=T)
}else{
  data_flink<-read.table(paste0("Posterior_alphas_group_",index,".txt"),sep="",header=T)
}

flink_sel<-vector(mode="numeric",length=numloci)
rows<-nrow(data_flink)

if(sel == "divergent"){
  flink_sel=colSums(data_flink[,2:(numloci+1)]>0)/rows
}else{
  flink_sel=colSums(data_flink[,2:(numloci+1)]<0)/rows
}

lim=0
locus=1
numreg=0
while(locus <= (numloci-1)){
  j=0
  if(flink_sel[locus] >= prob){
    selection=sel
    j=j+1
    numreg=numreg+1
    while((flink_sel[locus+j] >= 0.99) & (locus+j<=numloci)){
      j=j+1
    }
    start=as.numeric(freq_pos$V2[locus])
    end=as.numeric(freq_pos$V2[locus+j-1])
    diff=end-start
    if(diff>=lim){
      i=1
      while(i<=j){
        write(c(name, numreg, as.character(freq_pos$V1[locus]), freq_pos$V2[locus], flink_sel[locus], selection, group), file=paste("coefSel",sel,group,".txt",sep=""), ncolumns=7, append=TRUE, sep = "\t")
        locus=locus+1
        i=i+1
      }
    }
  }
  locus=locus+1
}