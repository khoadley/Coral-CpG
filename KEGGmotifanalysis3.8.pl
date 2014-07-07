#!/usr/bin/perl
use strict;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -

#This script takes sequences from a FASTA file and runs a running average for GC content, CG and GGG motifs for each nt within the sequence.
#This script also catalgoues the occurence of N values in each sequence for downstream QC. Output is placed into the 01-Data folder

# - - - - - U S E R    V A R I A B L E S - - - - - - - -

my $infolder = "01-KEGG";
my $outfolder = "02-KEGG";
my $outfile1  = "GCcontent.txt";
my $outfile2  = "CpGmotif.txt";
my $outfile3  = "GGGmotif.txt";
my $outfile4 = "Ndata.txt";
my $outfile5 = "GCmotif.txt";
my $outfile6 = "CpGnorm.txt";
my $outfile7 = "gpcnorm.txt";

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

my %Genes;
my @GCdata;
my @CGdata;
my @GCmdata;
my @GGGdata;
my @Ndata;
my @CpGnorm;
my @gpcnorm;
my $sequencecount = 0;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\n . . . Making data up . . .\n\n";

opendir(DATA, $infolder) or die "\n\nNADA $infolder, idiot . . . . \n\n";
my @filelist = readdir(DATA);
foreach my $file(@filelist)
{	if ($file =~ m/^\w/)
	{	print "processing $file\n";
		#my $pathfolder = "$file";
		system("mkdir $outfolder/$file");
		open(OUT1, ">".$outfolder."/".$file."/".$outfile1);
		open(OUT2, ">".$outfolder."/".$file."/".$outfile2);
		open(OUT3, ">".$outfolder."/".$file."/".$outfile3);
		open(OUT5, ">".$outfolder."/".$file."/".$outfile5);
		open(OUT6, ">".$outfolder."/".$file."/".$outfile6);
		open(OUT7, ">".$outfolder."/".$file."/".$outfile7);
		
		#my $promoter = "$file-.txt";
		&FASTAread($infolder."/".$file."/"."$file-.txt");
		foreach my $id (keys %Genes)
		{   if ($Genes{$id}{'SIZE'} > 3550)
		    {   $sequencecount += 1;
			my $position = 25;
			while ($position < 3525)
			{   my $window = substr($Genes{$id}{'ntseq'}, $position-25, 51); #sets window size for runnng average (20)
			    my $winsize = length($window);
			    my $screen = 'N';
			    if ($window =~ m/$screen/i)
			    {	    push (@GCdata, '-');
				    push (@CGdata, '-');
				    push (@GCmdata, '-');
				    push (@GGGdata, '-');
				    push (@CpGnorm, '-');
				    push (@gpcnorm, '-');
				    $position += 1;
			    }
			    else
			    {
				    # %GC calculator
				    my $count = 0;
				    $count = $window =~ tr/GC/GC/;
				    my $GC = $count/$winsize;
				    push (@GCdata, $GC);
				    
				    # CpG motif calculator
				    my $CGcount = 0;
				    my $char = 'CG';
				    my $pos = index($window, $char, 0);
				    while ($pos > 0)
				    {   $CGcount += (1/(1 + abs(($winsize/2)-$pos))); #allows for values to be weighed depending on distance from position within window
					$pos = index($window, $char, $pos+1);
				    }
				    push (@CGdata, &RoundOff(($CGcount/$winsize), 20));
				    
				    #GC motif calculator
				    my $GCcount = 0;
				    my $chart = 'GC';
				    my $posit = index($window, $chart, 0);
				    while ($posit > 0)
				    {   $GCcount += (1/(1 + abs(($winsize/2)-$posit))); #allows for values to be weighed depending on distance from position within window
					$posit = index($window, $chart, $posit+1);
				    }
				    push (@GCmdata, &RoundOff(($GCcount/$winsize), 20));
				    
				    #CpGnorm
				    my $CpGnormalized = (&RoundOff(($CGcount/$winsize), 10)/(0.00001 + &RoundOff($GC, 10))); 
				    push (@CpGnorm, &RoundOff($CpGnormalized, 20));
				
				    #gpcnorm
				    my $gpcnormalized = (&RoundOff(($GCcount/$winsize), 10)/(0.00001 + &RoundOff($GC, 10)));
				    push (@gpcnorm, &RoundOff($gpcnormalized, 20));
				    
				    # GGG motif calculator
				    my $GGGcount = 0;
				    my $charac = 'GGG';
				    my $posi = index($window, $charac, 0);
				    while ($posi > 0)
				    {   $GGGcount += (1/(1 + abs(($winsize/2)-$posi))); #allows for values to be weighed depending on distance from position within window;
					$posi = index($window, $charac, $posi+1);
				    }
				    push (@GGGdata, &RoundOff(($GGGcount/$winsize), 20)); #normalizes $GGGcount with respect to the size of $winsize
				    $position += 1;
			    }       
			}
	    
			 print OUT1 ">$Genes{$id}{'HEAD'} @GCdata\n";
			 print OUT2 ">$Genes{$id}{'HEAD'} @CGdata\n";
			 print OUT3 ">$Genes{$id}{'HEAD'} @GGGdata\n";
			 print OUT5 ">$Genes{$id}{'HEAD'} @GCmdata\n";
			 print OUT6 ">$Genes{$id}{'HEAD'} @CpGnorm\n";
			 print OUT7 ">$Genes{$id}{'HEAD'} @gpcnorm\n";
			 @GCdata = ();   #reset for next sequence
			 @CGdata = ();   #reset for next sequence
			 @GCmdata = ();  #reset for next sequence
			 @GGGdata = ();  #reset for next sequence
			 @CpGnorm = ();  #reset for next sequence
			 @gpcnorm = ();  #reset for next sequence
		    }
			 
		}
		close (OUT1);
		close (OUT2);
		close (OUT3);
		close (OUT5);
		close (OUT6);
		close (OUT7);
		print "there are $sequencecount in the file $file\n";
		$sequencecount = 0;
		
		#N-value QC check
		open(OUT4, ">".$outfolder."/".$file."/".$outfile4);
		&FASTAread($infolder."/".$file."/"."$file-.txt");
		foreach my $id (keys %Genes)
		{   if ($Genes{$id}{'SIZE'} == 3550)
			{	my $Nposition = 25;
				while ($Nposition < 3525)
				    {	my $Nvalue = substr($Genes{$id}{'ntseq'}, $Nposition, 1);
					my $Ncount = 0;
					if ($Nvalue =~ m/^N/i)
					{    $Ncount += 1;
					}
					push(@Ndata, $Ncount);
					$Nposition += 1;
				    }
				print OUT4 ">$Genes{$id}{'HEAD'} @Ndata\n";
				@Ndata = ();
			}
		}

		close (OUT4);

	undef %Genes;
	}
}
	 
	 
	 
print "\n\n\n. . . . . .DONE. . . . . . \n\n\n ";


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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

sub RoundOff
{	my $num = shift(@_);
	my $decimals = shift(@_);
	my $roundnum = int(($num * 10**$decimals) + 0.5)/(10**$decimals);
	return $roundnum;
}

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - -
#
#
