#!/bin/bash 
#Bash script to concatenate exon fasta files 
#Loop through each cancer cohort
for s in 'PROSTATE' 'BRCA' 'COLO' 'BONE_MARROW' 'COLON' 'LIVER' 'SPINAL_CORD' 'THYMUS' 'FETAL_LIVER' 'STOMACH' 'SMALL_INTESTINE' 'HEART' 'PLACENTA' 'SPLEEN' 'LUNG' 'BRAIN' 'KIDNEY' 'TESTIS' 'UTERUS' 'BONE' 'LYMPH_NODE' 'MM' 'PRAD' 'SARC' 'BLCA' 'CHOL' 'ACC' 'STAD' 'PAAD' 'SKCM' 'SECR' 'HNSC' 'OV' 'THCA' 'THYM' 'KDNY' 'ESCA' 'HCC' 'TGCT' 'MPN' 'NHL' 'LEUK' 'GBM' 'UCEC' 'PANC' 'MESO' 'ALL' 'RHABDO' 'LYMP' 'NRBL' 'NPBL' 'AML' 'HPBL' 'TLYM' 'LCH' 'MBL' 'JMML' 'ATRT' 'PNET' 'LGG' 'PLATELET' 'SKIN' 'OTHER'

do

echo $s

#Separate the FASTA sequences based on strand.

grep -A1 -- "(+)" "junction_bed/"$s".circRNA_exon_pos_three_prime_uniq.fasta" | grep -v -- '--' > "junction_bed/"$s".circRNA_exon_pos_three_prime_uniq.positive.fasta"

grep -A1 -- "(-)" "junction_bed/"$s".circRNA_exon_pos_three_prime_uniq.fasta" | grep -v -- '--' > "junction_bed/"$s".circRNA_exon_pos_three_prime_uniq.negative.fasta"

grep -A1 -- "(+)" "junction_bed/"$s".circRNA_exon_pos_five_prime_uniq.fasta" | grep -v -- '--' > "junction_bed/"$s".circRNA_exon_pos_five_prime_uniq.positive.fasta"

grep -A1 -- "(-)" "junction_bed/"$s".circRNA_exon_pos_five_prime_uniq.fasta" | grep -v -- '--' > "junction_bed/"$s".circRNA_exon_pos_five_prime_uniq.negative.fasta"

#3' exon on positive strand will be concatenated to 5' exon on positive strand; 5' exon on negative strand will be concatenated to 3' exon on negative strand.

paste --delimiters='' "junction_bed/"$s".circRNA_exon_pos_three_prime_uniq.positive.fasta" "junction_bed/"$s".circRNA_exon_pos_five_prime_uniq.positive.fasta" | sed 's/>[^>]*//2g' > "junction_fasta/"$s".backsplice.junction.positive.fasta"

paste --delimiters='' "junction_bed/"$s".circRNA_exon_pos_five_prime_uniq.negative.fasta" "junction_bed/"$s".circRNA_exon_pos_three_prime_uniq.negative.fasta" | sed 's/>[^>]*//2g' > "junction_fasta/"$s".backsplice.junction.negative.fasta"

#Combine both strands to generate backspace junction FASTA file.
cat "junction_fasta/"$s".backsplice.junction.positive.fasta" "junction_fasta/"$s".backsplice.junction.negative.fasta" > "junction_fasta/final_fasta/"$s".backsplice.junction.fasta"


done


