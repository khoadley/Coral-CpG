#!/usr/bin/perl
use strict;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -

#This script takes sequences from a FASTA file and runs a running average for GC content, CG and GGG motifs for each nt within the sequence.
#This script also catalgoues the occurence of N values in each sequence for downstream QC. Output is placed into the 01-Data folder

# - - - - - U S E R    V A R I A B L E S - - - - - - - -
my $infolder = "00-Data"; #folder with fasta file
my $promoter  = "Nematostella-KEGGtotal.txt"; #fasta file
my $outfolder = "01-Data"; #folder with resulting data
my $outfile1  = "GCcontent.txt"; #1-16 are all data files with respective motifs and nucleotide content
my $outfile2  = "CpGmotif.txt";
my $outfile3  = "GGGmotif.txt";
my $outfile4 = "Ndata.txt";
my $outfile5 = "GCmotif.txt";
my $outfile6 = "CpGnorm.txt";
my $outfile7 = "gpcnorm.txt";
my $outfile8 = "percentG.txt";
my $outfile9 = "percentC.txt";
my $outfile10 = "ATmotif.txt";
my $outfile11 = "TAmotif.txt";
my $outfile12 = "ATcontent.txt";
my $outfile13 = "TAnorm.txt";
my $outfile14 = "ATnorm.txt";
my $outfile15 = "percentA.txt";
my $outfile16 = "percentT.txt";

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

my %Genes;
my @GCdata;
my @CGdata;
my @GCmdata;
my @GGGdata;
my @Ndata;
my @CpGnorm;
my @gpcnorm;
my @G;
my @C;
my @ATmotif;
my @TAmotif;
my @ATcontent;
my @TAnorm;
my @ATnorm;
my @A;
my @T;

my $count = 0;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\n . . . Making data up . . .\n\n";

open(OUT1, ">".$outfolder."/".$outfile1);
open(OUT2, ">".$outfolder."/".$outfile2);
open(OUT3, ">".$outfolder."/".$outfile3);
open(OUT5, ">".$outfolder."/".$outfile5);
open(OUT6, ">".$outfolder."/".$outfile6);
open(OUT7, ">".$outfolder."/".$outfile7);
open(OUT8, ">".$outfolder."/".$outfile8);
open(OUT9, ">".$outfolder."/".$outfile9);
open(OUT10, ">".$outfolder."/".$outfile10);
open(OUT11, ">".$outfolder."/".$outfile11);
open(OUT12, ">".$outfolder."/".$outfile12);
open(OUT13, ">".$outfolder."/".$outfile13);
open(OUT14, ">".$outfolder."/".$outfile14);
open(OUT15, ">".$outfolder."/".$outfile15);
open(OUT16, ">".$outfolder."/".$outfile16);


&FASTAread($infolder."/".$promoter);
    foreach my $id (keys %Genes)
    {   if ($Genes{$id}{'SIZE'} == 3550)
        {   $count += 1;
	    print "$count gene sequences processed\n";
	    my $position = 25;
            while ($position < 3525)
            {   my $window = substr($Genes{$id}{'ntseq'}, $position-25, 51); #sets window size for runnng average (51)
                my $winsize = length($window);
		my $screen = 'N';
                if ($window =~ m/$screen/i)
		{	push (@GCdata, '-');
                        push (@CGdata, '-');
			push (@GCmdata, '-');
                        push (@GGGdata, '-');
			push (@CpGnorm, '-');
                        push (@gpcnorm, '-');
			push (@C, '-');
			push (@G, '-');
                        push (@ATmotif, '-');
			push (@TAmotif, '-');
                        push (@ATcontent, '-');
			push (@TAnorm, '-');
                        push (@ATnorm, '-');
			push (@A, '-');
			push (@T, '-');
                        $position += 1;
                }
                else
		{
			# %GC calculator
			my $count = 0;
			$count = $window =~ tr/GC/GC/;
			my $GC = $count/$winsize;
			push (@GCdata, &RoundOff($GC, 4));
			
			# %C calculator
			my $Ccount = 0;
			$Ccount = $window =~ tr/C/C/;
			my $C = $Ccount/$winsize;
			push (@C, &RoundOff($C, 4));
			
			# %G calculator
			my $Gcount = 0;
			$Gcount = $window =~ tr/G/G/;
			my $G = $Gcount/$winsize;
			push (@G, &RoundOff($G, 4));
			
			# %A calculator
			my $Acount = 0;
			$Acount = $window =~ tr/A/A/;
			my $A = $Acount/$winsize;
			push (@A, &RoundOff($A, 4));
			
			# %A calculator
			my $Tcount = 0;
			$Tcount = $window =~ tr/T/T/;
			my $T = $Tcount/$winsize;
			push (@T, &RoundOff($T, 4));
			
			# %AT calculator
			my $ATcount = 0;
			$ATcount = $window =~ tr/AT/AT/;
			my $AT = $ATcount/$winsize;
			push (@ATcontent, &RoundOff($AT, 4));
			
			# CpG motif calculator
			my $CGcount = 0;
			my $char = 'CG';
			my $pos = index($window, $char, 0);
			while ($pos > 0)
			{   $CGcount += (1/(1 + abs(($winsize/2)-$pos))); #allows for values to be weighed depending on distance from position within window
			    $pos = index($window, $char, $pos+1);
			}
			push (@CGdata, &RoundOff(($CGcount/$winsize), 4));
			
			#GC motif calculator
			my $GCcount = 0;
			my $chart = 'GC';
			my $posit = index($window, $chart, 0);
			while ($posit > 0)
			{   $GCcount += (1/(1 + abs(($winsize/2)-$posit))); #allows for values to be weighed depending on distance from position within window
			    $posit = index($window, $chart, $posit+1);
			}
			push (@GCmdata, &RoundOff(($GCcount/$winsize), 4));
			
			#AT motif calculator
			my $ATmcount = 0;
			my $ATchart = 'AT';
			my $ATposit = index($window, $ATchart, 0);
			while ($ATposit > 0)
			{   $ATmcount += (1/(1 + abs(($winsize/2)-$ATposit))); #allows for values to be weighed depending on distance from position within window
			    $ATposit = index($window, $ATchart, $ATposit+1);
			}
			push (@ATmotif, &RoundOff(($ATmcount/$winsize), 4));
			
			#TA motif calculator
			my $TAmcount = 0;
			my $TAchart = 'TA';
			my $TAposit = index($window, $TAchart, 0);
			while ($TAposit > 0)
			{   $TAmcount += (1/(1 + abs(($winsize/2)-$TAposit))); #allows for values to be weighed depending on distance from position within window
			    $TAposit = index($window, $TAchart, $TAposit+1);
			}
			push (@TAmotif, &RoundOff(($TAmcount/$winsize), 4));
			
			#CpGnorm
			my $CpGnormalized = (&RoundOff(($CGcount/$winsize), 10)/(0.00001 + &RoundOff($GC, 10))); 
			push (@CpGnorm, &RoundOff($CpGnormalized, 6));
			
			#gpcnorm
			my $gpcnormalized = (&RoundOff(($GCcount/$winsize), 10)/(0.00001 + &RoundOff($GC, 10)));
			push (@gpcnorm, &RoundOff($gpcnormalized, 6));
			
			#ATnorm
			my $ATnormalized = (&RoundOff(($ATmcount/$winsize), 10)/(0.00001 + &RoundOff($AT, 10))); 
			push (@ATnorm, &RoundOff($ATnormalized, 6));
			
			#TAnorm
			my $TAnormalized = (&RoundOff(($TAmcount/$winsize), 10)/(0.00001 + &RoundOff($AT, 10)));
			push (@TAnorm, &RoundOff($TAnormalized, 6));
			
			# GGG motif calculator
			my $GGGcount = 0;
			my $charac = 'GGG';
			my $posi = index($window, $charac, 0);
			while ($posi > 0)
			{   $GGGcount += (1/(1 + abs(($winsize/2)-$posi))); #allows for values to be weighed depending on distance from position within window;
			    $posi = index($window, $charac, $posi+1);
			}
			push (@GGGdata, &RoundOff(($GGGcount/$winsize), 4)); #normalizes $GGGcount with respect to the size of $winsize
			$position += 1;
		}       
            }

             print OUT1 ">$Genes{$id}{'HEAD'} @GCdata\n";
             print OUT2 ">$Genes{$id}{'HEAD'} @CGdata\n";
             print OUT3 ">$Genes{$id}{'HEAD'} @GGGdata\n";
	     print OUT5 ">$Genes{$id}{'HEAD'} @GCmdata\n";
	     print OUT6 ">$Genes{$id}{'HEAD'} @CpGnorm\n";
	     print OUT7 ">$Genes{$id}{'HEAD'} @gpcnorm\n";
	     print OUT8 ">$Genes{$id}{'HEAD'} @C\n";
             print OUT9 ">$Genes{$id}{'HEAD'} @G\n";
	     print OUT10 ">$Genes{$id}{'HEAD'} @ATmotif\n";
	     print OUT11 ">$Genes{$id}{'HEAD'} @TAmotif\n";
	     print OUT12 ">$Genes{$id}{'HEAD'} @ATcontent\n";
	     print OUT13 ">$Genes{$id}{'HEAD'} @TAnorm\n";
	     print OUT14 ">$Genes{$id}{'HEAD'} @ATnorm\n";
	     print OUT15 ">$Genes{$id}{'HEAD'} @A\n";
	     print OUT16 ">$Genes{$id}{'HEAD'} @T\n";
             @GCdata = ();   #reset for next sequence
             @CGdata = ();   #reset for next sequence
	     @GCmdata = ();  #reset for next sequence
             @GGGdata = ();  #reset for next sequence
	     @CpGnorm = ();  #reset for next sequence
             @gpcnorm = ();  #reset for next sequence
	     @C = ();   #reset for next sequence
	     @G = ();  #reset for next sequence
             @ATmotif = ();  #reset for next sequence
	     @TAmotif = ();  #reset for next sequence
             @ATcontent = ();  #reset for next sequence
	     @TAnorm = ();  #reset for next sequence
             @ATnorm = ();  #reset for next sequence
	     @A = ();   #reset for next sequence
	     @T = ();  #reset for next sequence
        }
    }

close (OUT1);
close (OUT2);
close (OUT3);
close (OUT5);
close (OUT6);
close (OUT7);
close (OUT8);
close (OUT9);
close (OUT10);
close (OUT11);
close (OUT12);
close (OUT13);
close (OUT14);
close (OUT15);
close (OUT16);

print "A total of $count sequences were processed\n\n";
print "starting QC check for N fraction\n\n";
    
open(OUT4, ">".$outfolder."/".$outfile4);
&FASTAread($infolder."/".$promoter);
    foreach my $id (keys %Genes)
    {   if ($Genes{$id}{'SIZE'} == 3550)
        {   my $position = 25;
	    while ($position < 3525)
	    {	my $Nvalue = substr($Genes{$id}{'ntseq'}, $position, 1);
		my $Ncount = 0;
		if ($Nvalue =~ m/^N/i)
		{    $Ncount += 1;
		}
		push(@Ndata, $Ncount);
		$position += 1;
	    }
	    print OUT4 ">$Genes{$id}{'HEAD'} @Ndata\n";
	    @Ndata = ();
	}
    }	 
close (OUT4);	 
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

sub RoundOff
{	my $num = shift(@_);
	my $decimals = shift(@_);
	my $roundnum = int(($num * 10**$decimals) + 0.5)/(10**$decimals);
	return $roundnum;
}

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - -