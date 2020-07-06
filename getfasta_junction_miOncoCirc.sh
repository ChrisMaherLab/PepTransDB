#!/bin/bash

#Loop through each cancer cohort
for s in 'PROSTATE' 'BRCA' 'COLO' 'BONE_MARROW' 'COLON' 'LIVER' 'SPINAL_CORD' 'THYMUS' 'FETAL_LIVER' 'STOMACH' 'SMALL_INTESTINE' 'HEART' 'PLACENTA' 'SPLEEN' 'LUNG' 'BRAIN' 'KIDNEY' 'TESTIS' 'UTERUS' 'BONE' 'LYMPH_NODE' 'MM' 'PRAD' 'SARC' 'BLCA' 'CHOL' 'ACC' 'STAD' 'PAAD' 'SKCM' 'SECR' 'HNSC' 'OV' 'THCA' 'THYM' 'KDNY' 'ESCA' 'HCC' 'TGCT' 'MPN' 'NHL' 'LEUK' 'GBM' 'UCEC' 'PANC' 'MESO' 'ALL' 'RHABDO' 'LYMP' 'NRBL' 'NPBL' 'AML' 'HPBL' 'TLYM' 'LCH' 'MBL' 'JMML' 'ATRT' 'PNET' 'LGG' 'PLATELET' 'SKIN' 'OTHER'

#Use getfasta to obtain fasta sequences of both the 3' and 5' exons of each circular RNA from BED files generated by <extractSeq_miOncoCirc.r>

do
echo $s
cmd="bedtools getfasta -s -name -fi /gscmnt/gc2608/maherlab/sidizhao/circ_rna/cabanski_cdna_capture/circexplorer2_output/all-chrs.fa -bed /gscmnt/gc2601/maherlab/sidizhao/circ_rna/proteomics_miOncoCirc/junction_bed/"$s".circRNA_exon_pos_three_prime_uniq.bed > /gscmnt/gc2601/maherlab/sidizhao/circ_rna/proteomics_miOncoCirc/junction_bed/"$s".circRNA_exon_pos_three_prime_uniq.fasta;
bedtools getfasta -s -name -fi /gscmnt/gc2608/maherlab/sidizhao/circ_rna/cabanski_cdna_capture/circexplorer2_output/all-chrs.fa -bed /gscmnt/gc2601/maherlab/sidizhao/circ_rna/proteomics_miOncoCirc/junction_bed/"$s".circRNA_exon_pos_five_prime_uniq.bed > /gscmnt/gc2601/maherlab/sidizhao/circ_rna/proteomics_miOncoCirc/junction_bed/"$s".circRNA_exon_pos_five_prime_uniq.fasta"

LSF_DOCKER_PRESERVE_ENVIRONMENT=true bsub -J $s -oo logs/$s.getfasta.out -eo logs/$s.getfasta.err -R"select[mem>32000 && gtmp>50] rusage[mem=32000]" -M 32000000 -q research-hpc -a "docker(registry.gsc.wustl.edu/genome/genome_perl_environment)" /bin/bash -c "$cmd"


done
