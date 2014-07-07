#!/usr/bin/perl
use strict;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
# MAIN FASTA PIPE
# Compile general metrics on the annotated ORFs in a microbial genome.
# Vars for each genome: header, ntseq, aaseq, size, GC, 
# MAST698-2013

# - - - - - U S E R    V A R I A B L E S - - - - - - - -
my $genome   = "transcripts.FilteredModels1.fasta";
my $infolder = "genome";

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -
my %Genes;                                  # Dicitionary with all profile stats
my %CodonTable;                             # translation table
my $outfolder = $infolder;
my $genomefolder;

if ($#ARGV > 0)
{	$genomefolder = $ARGV[0];
	$genome       = $ARGV[1];
	$infolder  = "genome";
	$outfolder = "01-cdna-data";
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\nRunning Profile Calcs on genome:\n      $genome\n\n";

# 1. Run the profile calculations . . . . . . . . . . .
&FASTAread($infolder."/".$genomefolder."/".$genome);   # Defines head, ntseq, size
&GCcalc(\%Genes);                  # Defines GC
&TranslateFasta(\%Genes);          # Defines aaseq                                > (requires LoadCodonTable)
&AAfreq(\%Genes);                  # Defines AAfreq.pX where X = each amino acid  > (requires TranslateFasta)
&AAmass(\%Genes);                  # Defines AAmass                               > (requires AAfreq)
&AAnitrogen(\%Genes);		       # Defines protein sidechain nitrogen req       > (requires AAfreq)
&CodonFreq(\%Genes);               # Defines CDfreq.pXXX where XXX = each codon   > (requires LoadCodonTable)
&CAI(\%Genes);                     # Defines CAI  = codon adaptation index        > (requires CodonFreq)
&Hentropy(\%Genes);                # Defines Haa, Hcd                             > (requires AAfreq, CodonFreq)
&IsoElectricPoint(\%Genes);        # Defines IEP                                  > (requires AAfreq)
&SKEW(\%Genes);                    # Defines SKaa, SKcd

# 2. Output results . . . . . . .
my @ColID = qw | SIZE GC AAmass AAnitrogen CAI Haa Hcd IEP SKaa SKcd |;
my @format = qw | STD AA CD |;       # < options: STD, AA, CD;  each results in a separate output file
system("mkdir $outfolder/$genomefolder");
my $outfile = $outfolder."/".$genomefolder."/".$genomefolder.".ffn";
foreach my $outform (@format)
{	&DataTable(\%Genes, $outfile, \@ColID, $outform); }   # Inputs to DataTable: Dict, genome, cols, format




# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - S U B R O U T I N E S - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub FASTAread
{	print "   Reading file . . . \n";
	# 1. Load FIlE . . . . . . . . . .
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
		unless ($seq =~ /[RWYSKMBDHVN\.\-]/i)          # filter for non ATGC base pairs
		{	$Genes{$unid}{'HEAD'}    = $data[0];       # store header
			$Genes{$unid}{'ntseq'}   = uc($seq);       # store sequence
			$Genes{$unid}{'SIZE'}    = length($seq);   # store length
			$unid += 1;
		}
	}
	$/="\n";
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub TranslateFasta
{	print "   Amino Acid translation . . . \n";
	&LoadCodonTable();
	my $GenDict = $_[0];
	foreach my $id (keys %{$GenDict})
	{	my $prot = '';
		my $gene = ${$GenDict}{$id}{'ntseq'};
		while ($gene =~ s/^(...)//)
		{	$prot .= $CodonTable{$1}; }
		$prot =~ s/\*$//;   # drop terminal stop codon
		${$GenDict}{$id}{'aaseq'} = $prot;
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub GCcalc
{	print "   Counting GCs . . . \n";
	my $GenDict = $_[0];
	foreach my $id (keys %{$GenDict})
	{	my $count = 0;
		$count = ${$GenDict}{$id}{'ntseq'} =~ tr/GC/GC/;
		$count = $count/${$GenDict}{$id}{'SIZE'};
		${$GenDict}{$id}{'GC'} = &RoundOff($count,3);
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub CodonFreq
{	print "   Codon frequencies . . . \n";
	my $GenDict = $_[0];
	foreach my $id (keys %{$GenDict})
	{	my $seq = ${$GenDict}{$id}{'ntseq'};
		my $N = length($seq)/3;
		if ($N < 1) { print "\n\nError in aaseq length; see sub &CodonFreq\n\n"; die; }
		while ($seq =~ s/^(...)//) 
		{	${$GenDict}{$id}{'CDfreq'}{'n'.$1} += 1/$N;  }
		foreach my $codon (keys %CodonTable)
		{	${$GenDict}{$id}{'CDfreq'}{'n'.$codon} = &RoundOff(${$GenDict}{$id}{'CDfreq'}{'n'.$codon}, 4); }
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub AAfreq
{	print "   Amino Acid frequencies . . . \n";
	my @AAx = qw | A C D E F G H I K L M N P Q R S T V W Y |;
	my $GenDict = $_[0];
	foreach my $id (keys %{$GenDict})
	{	my $seq = ${$GenDict}{$id}{'aaseq'};
		my $N = length($seq);
		if ($N < 1) { print "\n\nError in aaseq length; see sub &AAfreq\n\n"; die; }
		while ($seq =~ s/^(.)//) 
		{	${$GenDict}{$id}{'AAfreq'}{'p'.$1} += 1/$N;  }
		foreach my $aa (@AAx)
		{	${$GenDict}{$id}{'AAfreq'}{'p'.$aa} = &RoundOff(${$GenDict}{$id}{'AAfreq'}{'p'.$aa}, 4); }
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub AAmass
{	print "   Protein mass . . . \n";
	my @AAx = qw | A C D E F G H I K L M N P Q R S T V W Y |;
	my %mass = qw | A	71.0788 R	156.1875	N	114.1038	D	115.0886	C	103.1388	E	129.1155	Q	128.1307	G	57.0519	H	137.1411	I	113.1594	L	113.1594	K	128.1741	M	131.1926	F	147.1766	P	97.1167	S	87.0782	T	101.1051	W	186.2132	Y	163.1760	V	99.1326 |;
	my $GenDict = $_[0];
	foreach my $id (keys %{$GenDict})
	{	my $N = length(${$GenDict}{$id}{'aaseq'});
		${$GenDict}{$id}{'AAmass'} = 0;
		foreach my $aa (@AAx)
		{	${$GenDict}{$id}{'AAmass'} += &RoundOff(0.001 * ${$GenDict}{$id}{'AAfreq'}{'p'.$aa} * $N * $mass{$aa}, 2); }
		#&DUMP(${$GenDict}{$id}{'AAmass'});
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub AAnitrogen
{	print "   Protein nitrogen requirement . . . \n";
	my @AAx = qw | A C D E F G H I K L M N P Q R S T V W Y |;
	my %nitrog = qw | R	 3	N	1	Q	1	H	2	K	1	P	1	W	1  |;
	my $GenDict = $_[0];
	foreach my $id (keys %{$GenDict})
	{	my $W = length(${$GenDict}{$id}{'aaseq'});
		${$GenDict}{$id}{'AAnitrogen'} = 0;
		foreach my $aa (@AAx)
		{	${$GenDict}{$id}{'AAnitrogen'} += &RoundOff(0.001 * ${$GenDict}{$id}{'AAfreq'}{'p'.$aa} * $W * $nitrog{$aa}, 2); }
		#&DUMP(${$GenDict}{$id}{'AAnitrogen'});
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub CAI
{	print "   Codon Adaptation Index . . . \n";
	my $GenDict = $_[0];
	foreach my $id (keys %{$GenDict})
	{	# 1. Find maxfj = max(fi) for any group of synonymous codons (i.e., codons with the same AA)
		my %maxfj;
		my @AAx = qw | A C D E F G H I K L M N P Q R S T V W Y |;
		foreach my $aa (@AAx)                                                     # < Look at each AA . . . . 
		{	my $max = -1;                                                         # < Initialize with min value . . 
			$maxfj{$aa} = 0;                                                      # < Initialize with min value . . 
			foreach my $codon (keys %CodonTable)                                  # < Look at each Codon . . . 
			{	if ($CodonTable{$codon} == $aa)                                   # < When current codon = current AA . . . 
				{	if (${$GenDict}{$id}{'CDfreq'}{'n'.$codon} > $max)            # < Check if codon freq is max . . .
					{	$maxfj{$aa} = ${$GenDict}{$id}{'CDfreq'}{'n'.$codon}; }
				}
			}
		}
		# 2. CAI Calculation . . . . . . .
		my $seq = ${$GenDict}{$id}{'ntseq'};
		my $N = length(${$GenDict}{$id}{'aaseq'});  # total codons
		my $cai = 0;
		while($seq =~ s/^(...)//)
		{	my $fi = ${$GenDict}{$id}{'CDfreq'}{'n'.$1};
			my $aa = $CodonTable{$1};
			if ($maxfj{$aa} > 0 && $fi > 0)
			{	$cai += log($fi/$maxfj{$aa}); }
		}
		$cai = exp($cai/$N);
		${$GenDict}{$id}{'CAI'} = &RoundOff($cai, 3);
		#&DUMP(${$GenDict}{$id}{'CAI'}); 
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub Hentropy
{	print "   Entropy calculation . . . \n";
	my @AAx = qw | A C D E F G H I K L M N P Q R S T V W Y |;
	my $GenDict = $_[0];
	foreach my $id (keys %{$GenDict})
	{	# 1. Calculate Haa . . . . . . .
		my $haa = 0;
		foreach my $aa (@AAx)
		{	my $pi = ${$GenDict}{$id}{'AAfreq'}{'p'.$aa};
			if ($pi > 0)
			{	$haa += -1 * ($pi * log($pi)); }
		}
		${$GenDict}{$id}{'Haa'} = &RoundOff(10**$haa, 2);
		# &DUMP(${$GenDict}{$id}{'Haa'});
		
		# 2. Calculate Hcd . . . . . . .
		my $hcd = 0;
		foreach my $codon (keys %CodonTable)
		{	my $pi = ${$GenDict}{$id}{'CDfreq'}{'n'.$codon};
			if ($pi > 0)
			{	$hcd += -1 * ($pi * log($pi)); }
		}
		${$GenDict}{$id}{'Hcd'} = &RoundOff(10**$hcd, 2);
		# &DUMP(${$GenDict}{$id}{'Hcd'});
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub IsoElectricPoint
{	print "   Protein isoelectric point . . . \n";
	my $GenDict = $_[0];
	foreach my $id (keys %{$GenDict})
	{	my $N      = length(${$GenDict}{$id}{'aaseq'});
		my $pH     = 6.5;           # starting point pI = 6.5 - theoretically it should be 7, but
		my $pHlow  = 0.0;
		my $pHhigh = 14.0;  
		my $std    = 0.01;          # [pI = pH ± std]
		my $temp   = 0.0;
		(my $QN1, my $QN2, my $QN3, my $QN4, my $QN5, my $QP1, my $QP2, my $QP3, my $QP4, my $NQ) = (0,0,0,0,0,0,0,0,0,0);
		my $SOLVED = 0;
		my $looper = 0;
		while($SOLVED == 0)
		{	# 1. Calculate charge balance coefficient $NQ . . . . . .   
			$QN1 = -1/(1+10**(3.65-$pH));                                        
			$QN2 = (-1 * ${$GenDict}{$id}{'AAfreq'}{'pD'} * $N)/(1+10**(3.90-$pH));           
			$QN3 = (-1 * ${$GenDict}{$id}{'AAfreq'}{'pE'} * $N)/(1+10**(4.07-$pH));           
			$QN4 = (-1 * ${$GenDict}{$id}{'AAfreq'}{'pC'} * $N)/(1+10**(8.18-$pH));           
			$QN5 = (-1 * ${$GenDict}{$id}{'AAfreq'}{'pY'} * $N)/(1+10**(10.46-$pH));        
			$QP1 = (+1 * ${$GenDict}{$id}{'AAfreq'}{'pH'} * $N)/(1+10**($pH-6.040));            
			$QP2 = 1/(1+10**($pH-8.2));                
			$QP3 = (+1 * ${$GenDict}{$id}{'AAfreq'}{'pK'} * $N)/(1+10**($pH-10.54));           
			$QP4 = (+1 * ${$GenDict}{$id}{'AAfreq'}{'pR'} * $N)/(1+10**($pH-12.48));            
			$NQ = $QN1+$QN2+$QN3+$QN4+$QN5+$QP1+$QP2+$QP3+$QP4;
			
			# 2. Iteratively adjust pH value . . . . . . . .  .
			if ($NQ < 0)              # out of range, thus the new pH value must be smaller    
			{  $temp = $pH;
				$pH = $pH-(($pH-$pHlow)/2);
				$pHhigh = $temp;
			}
			else                      # pH value is too low, needs to be increased
			{   $temp = $pH;
				$pH = $pH + (($pHhigh-$pH)/2);
				$pHlow = $temp;
			}
			
			# 3. Check/Monitor progress . . . . . . . . . . . . . . . .  .. . .. 
			# Error, so stop trying to calculate, exit for loop
			if ($pH >= 14)
			{	${$GenDict}{$id}{'IEP'} = 'NA';
				$SOLVED = -1;
			} 
			# Iterative Solution: terminal condition, isoelectric point is within the given precision $std
			if (($pH-$pHlow < $std) && ($pHhigh-$pH < $std)) 
			{	${$GenDict}{$id}{'IEP'} = &RoundOff($pH, 3);
				$SOLVED = 1;
			}
			# Deadman switch to stop infinite loop . . . .
			if ($looper == 1000)
			{	${$GenDict}{$id}{'IEP'} = 'NA';
				$SOLVED = -1;
			}
			else
			{	$looper += 1; }
		} # end while loop . . . . . . . 
		# &DUMP(${$GenDict}{$id}{'IEP'});
	} # end foreach gene loop . . . . . . . . .. .  . .
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub SKEW                           # Requires the SKEWcalc subroutines . . 
{	print "   Skew calculation . . . \n";
	my @AAx = qw | A C D E F G H I K L M N P Q R S T V W Y |;
	my $GenDict = $_[0];
	foreach my $id (keys %{$GenDict})
	{	# 1. Calculate SKaa . . . . . . .
		my @pAA;
		foreach my $aa (@AAx)
		{	push(@pAA, ${$GenDict}{$id}{'AAfreq'}{'p'.$aa}); }
		${$GenDict}{$id}{'SKaa'} = &SKEWcalc(\@pAA);
		# &DUMP(${$GenDict}{$id}{'SKaa'});
		
		# 2. Calculate SKcd . . . . . . .
		my @nCD;
		foreach my $codon (keys %CodonTable)
		{	push(@nCD,${$GenDict}{$id}{'CDfreq'}{'n'.$codon}); }
		${$GenDict}{$id}{'SKcd'} = &SKEWcalc(\@nCD);
		# &DUMP(${$GenDict}{$id}{'SKcd'});
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub DataTable
{	my $GenDict = $_[0];
	my $outfile = $_[1];
	my $colid   = $_[2];
	my $job     = $_[3];
	my @AAx = qw | A C D E F G H I K L M N P Q R S T V W Y |;
	print "   Data dump to file . . . $job\n";
	
	# 1. Set the output file name . . . . . . . . . 
	if ($job eq 'STD')
	{	$outfile =~ s/\..{3}$/-data.txt/; }
	elsif ($job eq 'AA')
	{	$outfile =~ s/\..{3}$/-AAfreq.txt/; }
	elsif ($job eq 'CD')
	{	$outfile =~ s/\..{3}$/-CDfreq.txt/; }
	else
	{	print "\n\nFAIL: no file output!!\n&DataTable() takes 4 arguments: Dict, genome, colids, dataset\n\n"; }
	
	# 2. Print out the col headers for this data table . . . . . . . . . . 
	open(OUT,">$outfile");
	print OUT "ID";                       # first col header on line
	for my $var (@{$colid})              # now loop through each remaining col header
	{	print OUT "\t$var"; }    # need to print a tab-char '\t' and then new header
# What is the data output format?
	if ($job eq 'AA')
	{	foreach my $aa (@AAx) {	print OUT "\tp$aa"; } }
	elsif ($job eq 'CD')
	{	foreach my $codon (keys %CodonTable) {	print OUT "\tn$codon"; } }	
	print OUT "\n";                          # once complete, need to end the header line with a line-feed \n
	
	# 3. Now print the data/gene rows . . . . . . . . . 
	foreach my $id (keys %{$GenDict})
	{	print OUT "$id";
		foreach my $var (@{$colid})
		{	print OUT "\t${$GenDict}{$id}{$var}"; }
	# What is the data output format?
		if ($job eq 'AA')
		{	foreach my $aa (@AAx) {	print OUT "\t${$GenDict}{$id}{'AAfreq'}{'p'.$aa}"; } }
		elsif ($job eq 'CD')
		{	foreach my $codon (keys %CodonTable) {	print OUT "\t${$GenDict}{$id}{'CDfreq'}{'n'.$codon}"; } }	
		print OUT "\n";       # once row is complete, need to end the gene line with a line-feed \n		
	}
	close(OUT);
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub SKEWcalc                        # Requires the MEAN & MEDIAN subroutines . . 
{	my $array = $_[0];
	my $mean = &MEAN(\@{$array});
	my $median = &MEDIAN(\@{$array});
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	my $skew = &RoundOff(200*($median - $mean)/($mean + $median), 3);
	# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
	return $skew;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub MEAN                          # called by SKEW . . . . . . . . . . . . . . .
{	my $array = $_[0];
	my $N = scalar(@{$array});
	my $sum = 0;
	foreach my $data (@{$array})
	{	$sum += $data; }
	return &RoundOff($sum/$N, 4);
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub MEDIAN                        # called by SKEW . . . . . . . . . . . . . . .
{	my $array = $_[0];
	my $N = scalar(@{$array});
	my @SORT = sort{$a<=>$b}(@{$array});
	my $mid = int($N/2);
	return &RoundOff($SORT[$mid], 4);
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub LoadCodonTable
{	my @bases = ('T', 'C', 'A', 'G');
	my @aa = split(//,"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG");
	foreach my $a (@bases)
	{	foreach my $b (@bases)
		{	foreach my $c (@bases)
			{	$CodonTable{$a.$b.$c}=shift(@aa); }
		}
	}
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub LOG10
{	my $n = shift;
    return log($n)/log(10);
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub RoundOff
{	my $num = shift(@_);
	my $decimals = shift(@_);
	my $roundnum = int(($num * 10**$decimals) + 0.5)/(10**$decimals);
	return $roundnum;
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
sub DUMP
{	my $var = $_[0];
	print "\n\n>> $var << \n\n";
	die; 
}
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -