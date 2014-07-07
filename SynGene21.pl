#!/usr/bin/perl
use strict;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#this takes a predicted gene positions and finds them within a scaffold file. IF the position begins with and ATG, a 2.5k promoter region along with 500bp of the gene are dropped in to a file within the folder 00-Data
# - - - - - U S E R    V A R I A B L E S - - - - - - - -
my $infolder = "genome";
my $file  = "-CDfreq.txt";
my $outfolder = "00-Data";
my $outfile  = "Nemato-syntheticgenome.txt";

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

my %CODONFREQ = ();
my @CODONS = ();
my @s = ();
my $seqcount = 0;
my @p = ();
my @nuc = ("A", "T", "C", "G");
my %NUCFREQ = ();
my $gene = ();
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\n Making data up . . .\n\n";
print "\ncalculating codon frequencies\n";

$NUCFREQ{$nuc[0]} = 0.3003;
$NUCFREQ{$nuc[1]} = 0.2995;
$NUCFREQ{$nuc[2]} = 0.1981;
$NUCFREQ{$nuc[3]} = 0.2020;


my $infile = $infolder."/".$file;
open(IN, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
my @DATA = <IN>; close(IN);
my $head = shift(@DATA);
chomp $head;
my @headers = split (/\t/, $head);
foreach my $i(11..$#headers)
{     push (@CODONS, substr($headers[$i], 1,3));
}
foreach my $line (@DATA)
{	 my @data = split(/\t/, $line);
         foreach my $j(11..$#data)
         {  $CODONFREQ{$CODONS[$j-11]} += $data[$j];
            
         }
}
my $genes = $#DATA+1;
foreach my $codon (@CODONS)
{     $CODONFREQ{$codon} = $CODONFREQ{$codon}/$genes;
      print "$codon = $CODONFREQ{$codon}\n";
}

#this populates the @p array using nucleotide percent composition data
print "populating the promoter nucleotide table using percent composition\n";

foreach my $nucleotides (@nuc)
{     my $nucnum = 10000*$NUCFREQ{$nucleotides};
      my $nucleonum = int($nucnum);
      foreach (1..$nucleonum)
      {     push(@p, $nucleotides);
      }
}
print "@p\n";
#This populates the @s array using the codon frequency data generated above

print "populating codon table based on frequency information";
print "\n\n Making up straight BS. . .\n\n";

foreach my $codon (@CODONS)
{     my $codonnum = 10000*$CODONFREQ{$codon};
      my $num = int($codonnum);
      foreach (1..$num)
      {     push (@s, $codon);
      }
}
print "@s\n";

#This Generates the actual gene by first randomizing the @s array and then randomly picking elements from it to make the gene

open(OUT, ">".$outfolder."/".$outfile);

print "\n\n Making up straight BS. . .\n\n";

foreach my $sequence (1..20000)
{   $seqcount += 1;
    print OUT ">syntheticseq$seqcount\n";
    print "generating synthetic gene $seqcount\n";
    &FYshuffle;
    #die;
    #print "@s\n";
    foreach my $n (1..3024)
    {   my $rand = int(rand(9998));
        my $nucleotide = $p[$rand];
        #print "$rand ";
        #print OUT "$nucleotide";
        $gene = $gene . $nucleotide;
    }
    #print OUT "ATG";
    $gene = $gene . "ATG";
    foreach my $c (1..175)
    {   my $rand = int(rand(9948));
        my $codon = $s[$rand];
        $gene = $gene . $codon;
    }
    my $newgene = substr($gene, 0, 3550);
    #print "@gene\n\n\n";
    print OUT "$newgene\n";
    print "$newgene\n\n\n\n\n\n";
    #print OUT "\n";
    $gene = ();
    $newgene = ();
    sleep(1);
}
close OUT;
print "your new synthetic genome has $seqcount total genes\n";

print "\n\n. . . . . . DONE. . . . . . \n\n";




#subroutine
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
sub FYshuffle
{
	#Fisher-Yates Randomization . . . . .
	#my @s = split(//,$_[0]);
       # @s = shift;
	#print "@s\n";
	foreach my $k (1..10)
	{	#print "\n\n********************\n$k. Loop starting seq: ", join(", ", @s), "\n";
		for (my $i = $#s; $i>0 ; $i--)
		{	my $j = int(rand($i-1));
			@p[$i,$j] = @p[$j,$i];
			#print "---------------------------\n i= $i; j= $j; s= ", join("  |  ", @s),"\n";
		}
	}
	#return join('',@s);
        #print "@s\n";
}

die; 
print "\n\n\n. . . . .DONE. . . . .\n\n\n";

# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - -
#
#
