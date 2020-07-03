# PepTransDB
PepTransDB (Peptides in Transcripts Database, https://www.maherlab.com/peptransdb-lncrna) is a comprehensive integrative analyses of mass spectrometry-based proteomics and transcriptomic sequencing data from greater than 1,000 patients across ten cancer types. The resource has 21,027 novel peptides derived from 9,318 long noncoding RNAs (lncRNAs). We also exploited open reading frames overlapping the backspliced region of circular RNAs (circRNAs) to identify 3,515 peptides that are uniquely derived from 3,002 circRNAs and not their corresponding linear RNAs. We hope this pan-cancer proteogenomic analysis will serve as a resource for evaluating the coding potential of lncRNAs and circRNAs that could aide future mechanistic studies exploring their function in cancer. Here we provide scripts and steps used to generate the data.
## proteogenomic search and downstream analysis
pipeline.sh is a bash script that utilizes msFragger, peptideProphet and proteinProphet to generate .psm and .protein files for downstream analysis.  

Usage: sh pipeline.sh <directory_to_mzML_files> <database_file> <parameter_file> 
<directory_to_mzML_files> = directory for mzML files for each cancer.
<database_file> = database of target sequences with reverse decoy sequences with .fasta extension.
<parameter_file> = parameter file formatted in the format specified by msFragger (2019) and with .params extension. 

parse_proteogenomic_output.r is an R script used to parse the psm output for transcript-specific peptides for each cancer. 

Usage: Rscript parse_proteogenomic_output.r <directory_to_psm_files>
Args: <directory_to_psm_files> = location of .psm files directory by pipeline.sh

circ_orf_to_genome.py calculates coordinates of circular ORFs from EMBOSS getORF FASTA based on transcript locations (hg38).

sorf_to_genome.py calculates coordinates of short ORFs from EMBOSS getORF FASTA based on transcript locations (hg38).

Usage: python circ_orf_to_genome.py <ORF file> <exon file>

Usage: python sorf_to_genome.py <ORF file> <exon file>

Args: <ORF file> = reformatted FASTA file of ORFs from getORF e.g. headers should look like this: >ENST00000456328.2 7 - 63 not >ENST00000456328.2_1 [7 - 63]

<exon file> = contains exon sequence coordinates of transcripts ran through getORF columns (0 based numbering) 0:chromosome, 1:gene name, 2:transcript ID, 3:backsplice identifier, 4:strand, 5:exons (format start;stop,start;stop etc. )
## Supplementary steps

To obtain open reading frames through EMBOSS getORF
getorf  -sequence input_mRNA.fa -minsize 100  -outseq ORF.fa -table 1 -find 1 -noreverse -nomethionine
To obtain PhyloCSF scores
bigWigAverageOverBed https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg38/latest/PhyloCSF+3.bw  -minMax ORF.bed phylocsf_frame3
bigWigAverageOverBed https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg38/latest/PhyloCSF+2.bw  -minMax ORF.bed phylocsf_frame2
bigWigAverageOverBed https://data.broadinstitute.org/compbio1/PhyloCSFtracks/hg38/latest/PhyloCSF+1.bw  -minMax ORF.bed phylocsf_frame1  

## Generate junction sequences for CircRNA analysis
To extract coordinates of circRNAs from miOncoCirc
Rscript extractSeq_miOncoCirc.r
To get sequence of junction coordinates
sh getfasta_junction_miOncoCirc.sh
To concatinate junction exon sequences
sh getfasta_junction_miOncoCirc.sh 
 
## How to cite pepTransDB?
Coming soon...

## Contact
gothoum  (at) wustl (dot) edu
