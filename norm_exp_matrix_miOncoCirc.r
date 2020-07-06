#!/usr/bin/Rscript

#Read miOncoCirc data, sample info and sequencing depth files
miOncoCirc_output <- read.table("v0.1.release.txt",
                                header=T,sep="\t",stringsAsFactors=F,quote="");
sample_info <- read.table("Data/miOncoCirc_sample_info.txt",
                          header=T,sep="\t",stringsAsFactors=F,quote="");
seq_depth <- read.table("Data/miOncoCirc_seq_depth.txt",
                        header=T,sep="\t",stringsAsFactors=F,quote="");

#Normalize read numbers by sequencing depth of each individual sample
colnames(seq_depth)[1]<-"sample";
miOncoCirc_output<-merge(miOncoCirc_output,seq_depth,by="sample",all.x=T,all.y=F);
miOncoCirc_output$readNumberNormalized<-miOncoCirc_output$reads/miOncoCirc_output$Uniquely.Mapped.Reads*1000000;

miOncoCirc_output <- merge(miOncoCirc_output,sample_info,by.x="sample",by.y="ID",all.x=T);
miOncoCirc_output <- subset(miOncoCirc_output, !is.na(miOncoCirc_output$Sample.name));

miOncoCirc_output$release <- NULL;
miOncoCirc_output$Input.Reads <- NULL;

#Read exon annotation
exons = read.table("Data/exons_cleaned.txt",
                   header=T,sep="\t",stringsAsFactors=F,quote="");

#For each cancer cohort, get the read number and normalized read number of each sample
for (i in 1:length(unique(miOncoCirc_output$Analysis.Cohort)))
{
  cohort <- unique(miOncoCirc_output$Analysis.Cohort)[i];
  print(cohort);
  
  circRNA_exon_pos <- subset(miOncoCirc_output, miOncoCirc_output$Analysis.Cohort==cohort);
  colnames(circRNA_exon_pos)[3:4] <- c("five_prime_start","three_prime_end");
  
  circRNA_exon_pos_three_prime <- merge(circRNA_exon_pos,exons,
                                        by.x=c("chr","three_prime_end","symbol"),
                                        by.y=c("chrom","end","gene_name"),
                                        all.x=T, all.y=F)
  
  circRNA_exon_pos_five_prime <- merge(circRNA_exon_pos,exons,
                                       by.x=c("chr","five_prime_start","symbol"),
                                       by.y=c("chrom","start","gene_name"),
                                       all.x=T, all.y=F)
  
  merge_pos <- merge(circRNA_exon_pos_three_prime,circRNA_exon_pos_five_prime,by=c("chr","strand","symbol","transcript_id","sample","gene_id","transcript_name","Analysis.Cohort","five_prime_start","three_prime_end"));
  merge_pos <- unique(subset(merge_pos, !is.na(merge_pos$start)&!is.na(merge_pos$end)));
  
  #Generate unique identifiers for each circular RNA
  cols <- c("chr","five_prime_start","three_prime_end","strand","symbol","transcript_id");

  merge_pos$uniq_identifier <- apply( merge_pos[ , cols ] , 1 , paste , collapse = "|" );
  merge_pos$uniq_identifier <- gsub(" ","",merge_pos$uniq_identifier,fixed=TRUE);

  merge_pos <- merge_pos[,c("uniq_identifier","sample",
                            "reads.x","readNumberNormalized.x","Uniquely.Mapped.Reads.x")]
  colnames(merge_pos) <- c("uniq_identifier","sample",
                           "reads","readNumberNormalized","Uniquely.Mapped.Reads")
  
  #Write expression matrix for each cohort
  write.table(merge_pos,paste0("expr_matrix/",
                                                      cohort,".norm_expr_matrix.txt"),sep="\t",col.names = T,row.names = F,quote = F)
}


