#!/usr/bin/perl
use strict;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Re-format KEGG files
#this script reads in a kegg file and matches kegg annotated terms with a fasta file.
#resulting data is placed within the 00-KEGG folder.
# - - - - - U S E R    V A R I A B L E S - - - - - - - -

my $path    = "000-Capitella";
my $kegginfile  = "capitella-KEGG.txt"; #KEGG codes for each header annotated
my $fnainfile   = "CapitellaGenes-v4-NR-seq-3000.fasta";
my $outfolder = "00-KEGG";
my $outfile  = "kegg-$fnainfile"; #resulting file with only KEGG annotated sequences

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

my %Genes;
my $kegganotcounter = 0;
my $kegggenecounter = 0;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
print "\n\nWorking for the man every night and day . . .\n\n";

# 1. Load KEGG file . . . . . . . . . .
my $keggpathfile = $path . "/" . $kegginfile;
#$/='01|'; #shows where to break input string
open(IN, "<$keggpathfile") or die "\n\nNADA $keggpathfile you FOOL!!!\n\n";
my @keggDATA = <IN>;
close(IN);

print "$keggpathfile opened with ", $#keggDATA + 1," lines input to memory.\n\n";
#shift(@keggDATA);


# 2. Remove tabs and create array of KEGG headers (@NEWkegg)
my @NEWkegg;
foreach my $keggline (@keggDATA)
{   #print "keggline... $keggline\n";
    $keggline =~ s/\t//g;
    $keggline =~s/\n//g;
    #print "keggline is $keggline \n\n";
    push(@NEWkegg, $keggline);
}

open(OUT, ">".$outfolder."/".$outfile);
&FASTAread($path."/".$fnainfile);
foreach my $keggline (@NEWkegg)
{   my $newkegghdrmatch = $keggline;
    if ($newkegghdrmatch =~ m/(K\d+)$/)
    {   $kegganotcounter +=1;
	my $searchterm = $newkegghdrmatch;
        #print "$searchterm\n";
	foreach my $id (keys %Genes)
	{   my $genehead = $Genes{$id}{'HEAD'};
	    #my $head = substr($genehead, 8, 5);
	    #print "$head\n";
	    if ($searchterm =~ m/$genehead/)
	    {	$kegggenecounter += 1;
		print OUT ">$searchterm\n$Genes{$id}{'ntseq'}\n";
		print "$genehead\t$searchterm\n";
		last;
	    }
	}
    }
}
print "there are $kegganotcounter kegg annotations and $kegggenecounter matches in the file $fnainfile\n\n";
close(OUT);

print "\n\n\n* * * *   D O N E   * * * * \n\n\n";
# - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - -
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




