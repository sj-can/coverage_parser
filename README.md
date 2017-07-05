# coverage_parser
parses coverage metric file and plots with gnuplot

usage:
python coverage_parser.py -s sample_file -g gene_target_file

#Inputs
a sample list for all the samples requrign coverage of genes plotted
a gene list of all the genes requrieing coveage plotteing for each sample

#outputs
.png file showing coverage of each exon +/- 50 for each gene of each sample specified. 
Intemediate files retained for additional plotting as required

prequisites
python 2.7
gnuplot-py (1.8)

Detailed description.
1. Satisfies input files are correctly formatted
2. initalises Coverage_parser to read target samples and genes into memory and target alamut gene database file
3. Calls exome_coverage_finder_bash to return the coverage files associated with each exome identifier. 
	This is based on location of the coverage file relative to the control file. 
	The control file specifies each sample within a batch
4. Identifies the intervals for each gene from the alamut field, initalises and appends the transcript instances
5. Identifies the longest transcript for every gene
6. For each remaining transcript create an exon +/-50 interval file
7. Sort the interval file based on the transcript strand (forward\reverse)
8. Additional logic to create the complementart exon/interval data files for plotting
9. Finds coverage for each gene in each transcripts
10. plots with gnuplot


