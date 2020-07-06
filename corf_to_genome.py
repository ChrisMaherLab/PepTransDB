#!usr/bin/env python3

"""
Calculates coordinates of ORFs from getORF FASTA based on transcript locations (hg38). 

Usage: python corf_to_genome.py <ORF file> <exon file> 

Args: <ORF file> = reformatted FASTA file of ORFs from getORF 
ex headers should look like this: >ENST00000456328.2 7 - 63 not >ENST00000456328.2_1 [7 - 63]

<exon file> = contains exon sequence coordinates of transcripts ran through getORF 
columns (0 based numbering) 0:chromosome, 1:gene name, 2:transcript ID, 3:backsplice identifier, 4:strand, 5:exons (format start;stop,start;stop etc. )

"""
#import sys, matplotlib, numpy, os
import sys
import os
#used to save dictionaries to make script run faster
import pickle

#Check if an argument was passed to the python script
#If not exit exit and print the documentation
if (len(sys.argv)!=3):
        sys.exit(__doc__)

#to save dictionaries
def save_obj(obj, name ):
    with open('obj/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

#to load saved dictionaries
def load_obj(name ):
    with open('obj/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)

#exon_list format [[start,stop],[start,stop]] (inclusive of both)
#orf_location format [start,stop] (inclusive of both)
#takes in exon locations and open reading frame locations (based on one based 
#numbering from beginning of feature) and gives list of base locations of open reading frames
def find_regions(exon_list, orf_location, strand_sign):
	positions = []
	#add exon positions to positions list
	for i in exon_list:
		positions += range(i[0],i[1]+1)
	#get orfs based on position values
	if strand_sign == '-':
		#takes care of negatives being backwards since all transcripts written as 5' to 3' in origional FASTA
		positions = positions[::-1]
	orf = [p for p in positions[orf_location[0]-1:orf_location[1]]]
	return(orf)

#use for file naming later
basename = os.path.splitext(os.path.basename(sys.argv[1]))[0]

#list of  IDs
name_list = []
#dictionary of  IDs (key) and chromosomes (value)
chrom_dict = {}
#dictionary of  IDs (key) and their exons in a list of lists (value)
exon_dict = {}
#dictionary of  IDs (key) and strands (value)
strand_dict = {}
#dictionary of IDs (key) and gene name (value)
geneid_dict = {} #comment out if load saved instead
#get information from <exon file>
with open(sys.argv[2], 'r') as bed:
	for line in bed:
		linesplit = line.split()
		name_list.append(linesplit[3])
		chrom_dict[linesplit[3]] = linesplit[0]
		strand_dict[linesplit[3]] = linesplit[4]
		geneid_dict[linesplit[3]] = linesplit[1]
		#split by exons
		exon = linesplit[5].split(',')
		exon_list = []
		for i in exon:
			#format for use in find_regions later
			exon1 = i.split(';')
			exon2 = list(map(int,exon1))
			exon_list.append(exon2)	
		exon_dict[linesplit[3]] = exon_list

#if not saved comment out
#geneid_dict = load_obj('circ_geneid_dict_for_orfs')

#if have saved comment out this section ##################################################
#otherwise save dictionary  ID (value) geneid (key)
save_obj(geneid_dict, 'circ_geneid_dict_for_orfs')
########################################################################################

#make dictionary of fasta  headers (values) by transcript ID and orf location (key)
fasta_line_dict = {}
#make dictionary of transcript ID (keys) and list of orf positions (value) 
orf_dict = {}
#get info from <ORF file>
with open(sys.argv[1], 'r') as fasta:
	for line in fasta:
		if line.startswith('>'):
			linesplit = line.strip().split()
			name = str(linesplit[0]).replace('>','')   
			if name in strand_dict.keys():      
				if strand_dict[name] == '+':
					if name not in orf_dict.keys():
						orf_dict[name] = [[int(linesplit[1]),int(linesplit[3])]]
						fasta_line_dict[name+str(int(linesplit[1]))+str(int(linesplit[3]))] = line
					elif name in orf_dict.keys():
						orf_dict[name] += [[int(linesplit[1]),int(linesplit[3])]]
						fasta_line_dict[name+str(int(linesplit[1]))+str(int(linesplit[3]))] = line
				#reverse sense - can change if needed - ie if six frame translation was used
				elif strand_dict[name] == '-':
					if name not in orf_dict.keys():
						orf_dict[name] = [[int(linesplit[1]),int(linesplit[3])]]
						fasta_line_dict[name+str(int(linesplit[1]))+str(int(linesplit[3]))] = line
					elif name in orf_dict.keys():
						orf_dict[name] += [[int(linesplit[1]),int(linesplit[3])]]
						fasta_line_dict[name+str(int(linesplit[1]))+str(int(linesplit[3]))] = line
				else:
					continue
			else:
				continue

#name and open output file
outfile = open('new_total_'+basename+'.bed', 'w')

#start getting output
for key, value in orf_dict.items():
	if key in exon_dict.keys():
		for i in value:
			#get orf locations
			find_reg_out = find_regions(exon_dict[key], i, strand_dict[key])
			#change negative strand transcripts back to correct orientation 
			if strand_dict[key] == '-':
				find_reg_out = find_reg_out[::-1]
			#set the beginning of the orf 
			start = find_reg_out[0]
			#separates each continuous part of reading frame and writes to bed file
			for m in range(len(find_reg_out)-1):
				if find_reg_out[m+1] - find_reg_out[m] == 1: 
					continue
				else: 
					end = find_reg_out[m] 
					#output numbers should be compatible with UCSC browser
					outfile.write(chrom_dict[key]+'\t'+str(start)+'\t'+str(end)+'\t'+key+'\t'+str(geneid_dict[key])+'\t'+strand_dict[key]+'\t'+fasta_line_dict[str(key)+str(i[0])+str(i[1])])
					#resets start for next section of orf to be written 
					start = find_reg_out[m+1]
			#write current start and last base(-1) of orf
			outfile.write(chrom_dict[key]+'\t'+str(start)+'\t'+str(find_reg_out[-1])+'\t'+key+'\t'+str(geneid_dict[key])+'\t'+strand_dict[key]+'\t'+fasta_line_dict[str(key)+str(i[0])+str(i[1])])
outfile.close()

