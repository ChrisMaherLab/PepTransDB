#!/bin/bash
#the programs would only run from within the input directory with mzML files 
cd $1
#temporarily copy the database file to input directory
cp $2 .
mv *.fasta database.fasta
#temporarily copy the parameter file
cp $3 .
mv *.params msfragger.params
#clean any previous workspace 
philosopher workspace --clean
#initiate new workspace
philosopher workspace --init
#add the database created and indicate that decoy sequences as prefixed by XXX_
philosopher database --custom database.fasta  --prefix XXX_
#remove temporary database
rm -f database.fasta
#run msFragger
java  -jar MSFragger-20190628.jar  *.params *.mzML
#run peptideProphet using Accurate Mass model binning and semi-parametric modeling
philosopher peptideprophet --combine --nonparam  --ppm --accmass --decoy XXX_ --database database.fasta *.pepXML
#run proteinProphet 
ProteinProphet  interact.pep.xml output.protxml
#filter using an FDR or 0.01 at the peptide and protein levels
philosopher filter  --tag XXX_ --pepxml interact.pep.xml --protxml output.protxml --pep 0.01  --prot 0.01
#generate reports 
philosopher report
