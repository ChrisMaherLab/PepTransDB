#!/usr/bin/Rscript

#script to extract coordinates of circRNAs

#Read miOncoCirc data and cohort information
#miOncoCirc_v0.1.release.txt can be downloaded from https://mioncocirc.github.io/download/
#miOncoCirc_sample_info.txt and  exons_cleaned.txt are available from https://github.com/ChrisMaherLab/pepTransDB
miOncoCirc_output <- read.table("miOncoCirc_v0.1.release.txt",
                                header=T,sep="\t",stringsAsFactors=F,quote="");
sample_info <- read.table("miOncoCirc_sample_info.txt",
                          header=T,sep="\t",stringsAsFactors=F,quote="");
miOncoCirc_output <- merge(miOncoCirc_output,sample_info,by.x="sample",by.y="ID",all.x=T);
miOncoCirc_output <- subset(miOncoCirc_output, !is.na(miOncoCirc_output$Sample.name));

miOncoCirc_output$release <- NULL;

#Read exon annotation
exons = read.table("exons_cleaned.txt",
                   header=T,sep="\t",stringsAsFactors=F,quote="");

#For each cancer cohort, merge exon and circRNA data to obtain the individual start and end coordinates of
#the first and last exons of each circular transcript
for (i in 1:length(unique(miOncoCirc_output$Analysis.Cohort)))
{
  cohort <- unique(miOncoCirc_output$Analysis.Cohort)[i];
  print(cohort);
  
  miOncoCirc_output_cohort <- subset(miOncoCirc_output, miOncoCirc_output$Analysis.Cohort==cohort);
  circRNA_exon_pos <- miOncoCirc_output_cohort[,c("chr","start","end","symbol","sample","Analysis.Cohort")];
  colnames(circRNA_exon_pos)[2:3] <- c("five_prime_start","three_prime_end");
  
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
  
  #Generate unique identifier for each circular RNA
  cols <- c("chr","five_prime_start","three_prime_end","strand","symbol","transcript_id");
  merge_pos$uniq_identifier <- apply( merge_pos[ , cols ] , 1 , paste , collapse = "|" );
  merge_pos$uniq_identifier <- gsub(" ","",merge_pos$uniq_identifier,fixed=TRUE);
  
  #Generate bed files with exon coordinates of the 3' and 5' exons separately for each circular RNA
  circRNA_exon_pos_three_prime_bed <- merge_pos[,c(1,11,10,13,2:4)];
  circRNA_exon_pos_three_prime_bed$score <- ".";
  circRNA_exon_pos_three_prime_bed <- unique(circRNA_exon_pos_three_prime_bed[,c(1:4,8,5:7)]);
  circRNA_exon_pos_five_prime_bed <- merge_pos[,c(1,9,12,13,2:4)];
  circRNA_exon_pos_five_prime_bed$score <- ".";
  circRNA_exon_pos_five_prime_bed <- unique(circRNA_exon_pos_five_prime_bed[,c(1:4,8,5:7)]);

  write.table(circRNA_exon_pos_three_prime_bed,paste0("junction_bed/",
                                                      cohort,".circRNA_exon_pos_three_prime_uniq.bed"),sep="\t",col.names = F,row.names = F,quote = F)
  write.table(circRNA_exon_pos_five_prime_bed,paste0("junction_bed/",
                                                     cohort,".circRNA_exon_pos_five_prime_uniq.bed"),sep="\t",col.names = F,row.names = F,quote = F)
  
}



