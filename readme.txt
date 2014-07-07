WGMotifanalysis Scripts:

#The following scripts assume you are starting with a FASTA file that contains sequences consisting of a 3025bp promoter region #directly upstream up the transcriptional start site, along with the first 525 bp of the gene for a total sequence size of 3550bp.

                
#WGMotifanalysis-2.7.pl
	# This script takes sequences from a FASTA file and calculates a running average for %GC content, %AT content, CpG, GpC, ApT, 	# TpA, GGG #motifs. The motif frequencies are recorded as both raw data and also normalized to either %GC or %AT content. 
	# This script also catalogues the occurrence of N values in each sequence for downstream Quality Control.
	# Output is placed into the 01-Data folder

#WGMotifanalysis-3.6.pl
        # For each file in 01-Data, this script averages each row to give a mean number
        # for each position within the promoter/gene.
        # This script also accounts for N, removing bias due to lack of data within a specific gene.
        # Resulting data is placed in the 02-Data folder




KEGGmotifanalysis Scripts:

#this set of scripts calculates the same parameters as the WGMotifanalysis scripts. However, the initial FASTA file is separated into #different KEGG annotation groups for analysis of genes within different functional groups.


#KEGGmotifanalysis1.1.pl
       #Re-format KEGG files
       #this script reads in a kegg file and matches kegg annotated terms with a fasta file.
       #resulting data set includes only sequences from the original FASTA file with a KEGG annotation. Data is placed within the 
       #00-KEGG folder.
            
#Use Either 2.1 or 2.2 depending KEGG tier of interest
        
                #KEGGmotifanalysis2.1.pl
                    #This file takes a kegg annotated FASTA file and splits it into separate files based on gene
                    #grouping heirarchies within the file kegger.txt
                    #output is placed within the folder 01-KEGG
                    #Output is for tier 1 kegg annotations
                    
                #KEGGmotifanalysis2.2.pl
                    #This file takes a kegg annotated FASTA file and splits it into separate files based on gene
                    #grouping heirarchies within the file kegger.txt
                    #output is placed within the folder 01-KEGG
                    #Output is tier 2 annotations.
        
#KEGGmotifanalysis3.8.pl
	#This script takes sequences from a FASTA file and runs a running average for GC content, CG and GGG motifs for each nt within 	#the sequence.
	#This script also catalgoues the occurence of N values in each sequence for downstream QC. Output is placed into the 01-Data 	#folder
        
#KEGGmotifanalysis4.5.pl
	# This script transposes each file in 01-DATA so that each column equals a single gene. T
	# The script also averages each position and outputs it as the last column.
	# All data is placed into the 03-Data file
            
#KEGGmotifanalysis5.1.pl
	# This averages across a row and outputs the resulting means into the 04-kegg folder


SynGene Scripts:

# these scripts produce a randomized genome based on the nucleotide composition and codon frequency of the actual genome. This enables 
# the researcher to compare their data set in terms of selection against or for motif structures. The scripts require a FASTA file of # gene transcripts (for calculating codon frequency for the 525 bp region 3’ of the transcriptional start site) and nucleotide 
# composition (fraction of A, T, C, G nucleotides within the 3025 bp region 5’ of the transcriptional start site). The resulting FASTA # file, containing 20000 randomized genes is meant to then be run through the WGMotifanalysis scripts to produce the randomized genome data set for comparison with the actual data set.
