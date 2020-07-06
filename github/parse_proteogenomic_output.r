#!/usr/bin/env Rscript
#script to parse PSMs from pipeline
args = commandArgs(trailingOnly=TRUE)
# read input 
#files in the directory should be formatted as: cancer_psm.tsv
input_dir=args[1]
cancer=args[2]
list_of_files=list.files(input_dir)
for (f in 1:length(list_of_files))
{
  file_name=list_of_files[f]
  
  psm=read.csv(file_name, sep='\t', header=T, stringsAsFactors = F)	
  #remove decoy sequences 
  mapped_prots=subset(psm, !grepl("XXX_", Protein))
  #get mapped proteins 
  s=strsplit(as.character(psm$Mapped.Proteins), ",")  
  mapped_prots=mapped_prots[,c(17,25)]
  psm_prob=data.frame(PeptideProphet.Probability = rep(psm$PeptideProphet.Probability, sapply(s, length)), Protein = unlist(s))
  psm=unique(rbind(mapped_prots,psm_prob))
  psm$Protein <- sub("_[^_]+$", "", psm$Protein)
  colnames(psm)[2]='transcript'
  #write parsed output 
  write.table(psm,paste0('cancer','.tsv'),sep='\t', col.names=T,row.names=F, quote=F)
  
}
