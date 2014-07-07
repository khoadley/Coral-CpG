#!/usr/bin/perl
use strict;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
# 
#This file takes a kegg annotated FASTA file and splits it into separate files based on gene grouping heirarchies within the file kegger.txt
#output is placed within the folder 01-KEGG
#Output is for tier 1 kegg annotations



# - - - - - U S E R    V A R I A B L E S - - - - - - - -


# File I/O parameters . . . . . . . . 
my $folder       = "00-KEGG";           # data folder
my $outfolder    = "01-KEGG";

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %Genes;
my @KOcodes;
my @KeggPath	= ("Metabolism", "Genetic-Information-Processing", "Environmental-Information-Processing", "Cellular-Processes", "Organismal-Systems", "Human-Diseases");
# Kclass

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

foreach my $kpath (@KeggPath)
{	my $pathofinterest = "$kpath";  #input path name of interest.# files extension of target file for processing
	my $outfiletag   = "$pathofinterest-.txt";
	my $pathfolder = "$kpath";
	system("mkdir $outfolder/$pathfolder");
	#populate array with KO codes from specific pathway
	
	open (my $fh, "<", "kegger.txt");
	my @file_array;
	while (my $line = <$fh>)
	    {   chomp $line;
		my @values = split(/[\t]+/, $line);
		if ($values[2] eq $pathofinterest)   # Input heirarchy of interest $values[?]; 2=Kclass, 3=Kpath, 4=Kpathway
		    {my $cocode = $values[0];
		     push (@KOcodes, $cocode);
		     #print "$cocode\n";
		     #print "$values[0]\n";
		    }
	    }
		
	print "\n\nSearching for kegg annotated genomes in: $folder\n\n   Looking For $pathofinterest genes:\n";
	
	opendir(DH, "$folder");
	my @files = readdir(DH);
	#closedir(DH);
	open(OUT, ">".$outfolder."/".$pathfolder."/".$outfiletag);
	foreach my $file (@files)
	{       if ($file =~ m/^\w/)
		{	&FASTAread($folder."/".$file);
			foreach my $genetag (@KOcodes)
				{foreach my $id (keys %Genes)
					{	if ($Genes{$id}{'HEAD'} =~ m/$genetag/)
						{	print OUT ">$genetag\n$Genes{$id}{'ntseq'}\n";
						}
					}
				}
                }	
	undef %Genes;	
	}
	print "\n";
	closedir(DH);
	close(OUT);
        @KOcodes = ();
}
print "\n\n*  *  *  *  *    D  O  N  E     *  *  *  * \n\n";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub FASTAread
{	# 1. Load FIlE . . . . . . . . . .
	$/=">";                                     # set input break string
	my $infile = $_[0];
	open(IN, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
	my @DATA = <IN>; close(IN); shift(@DATA);	
	# 2. Parse sequence data . . . . . . . . . . . . .
	my $unid = 10000;                           # string to generate unique ids
	foreach my $entry (@DATA)
	{	my @data = split('\n', $entry);
		my $seq = '';
		foreach my $i (1..$#data)
		{	$seq .= $data[$i];  }
		$seq =~ s/>//;
		#unless ($seq =~ /[RWYSKMBDHVN\.\-]/i)          # filter for non ATGC base pairs
		{	$Genes{$unid}{'HEAD'}    = $data[0];       # store header
			$Genes{$unid}{'ntseq'}   = uc($seq);       # store sequence
			$Genes{$unid}{'SIZE'}    = length($seq);   # store length
			$unid += 1;
		}
	}
	$/="\n";
}

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - -
#
#
