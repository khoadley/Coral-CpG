#!/usr/bin/env perl
use strict;
# - - - - - H E A D E R - - - - - - - - - - - - - - - - -

# This script transposes each file in 01-DATA so that each column equals a single gene. 
# All data is placed into the 02-Data file

# - - - - - U S E R    V A R I A B L E S - - - - - - - -
my $infolder = "01-Data";
my $outfolder = "02-Data";

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

my %AVERAGE = ();
my %SAMPLENUM = ();
my @NTdata = ();
my @NTwd = ();
my $sequence = ();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\n . . . Making data up . . .\n\n";

opendir(DATA, $infolder) or die "\n\nNADA $infolder, idiot . . . . \n\n";
my @filelist = readdir(DATA);
foreach my $file(@filelist)
{   if ($file =~ m/^\w/)
    {	open(OUT1, ">".$outfolder."/".$file);
        my $infile = $infolder."/".$file;
        open(IN, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
        my @DATA = <IN>; close(IN);
        foreach my $i(1..3500)
        {       push (@NTdata, $i);
                push (@NTwd, $i);
                
        }
        #print "@NTdata\n\n\n\n\n@NTwd\n\n\n\n";
        foreach my $line (@DATA)
        {	my @datum = split(/ /, $line);
                $sequence += 1;
                foreach my $j(1..$#datum)
                {       if ($datum[$j] =~ m/^\d/)
                        {       #print "$datum[$j]\t";
                                $AVERAGE{$NTdata[$j-1]} += $datum[$j];
                                #print "$AVERAGE{$NTdata[$j-1]}\t";
                                $SAMPLENUM{$NTwd[$j-1]} += 1;
                                #print "$SAMPLENUM{$NTwd[$j-1]}\n";
                        }
                        else
                        {
                        }
                        #die;
                }
        #print "$file sequence number $sequence\n";
        }
        foreach my $c(1..3499)
        {       my $NTaveraged = $AVERAGE{$NTdata[$c]}/$SAMPLENUM{$NTwd[$c]};
                print "$file = $AVERAGE{$NTdata[$c]}\t$SAMPLENUM{$NTwd[$c]}\t$c\n";
                print OUT1 "$NTaveraged\n";
                $AVERAGE{$NTdata[$c]} = 0;
                $SAMPLENUM{$NTwd[$c]} = 0;
        }
    my %AVERAGE = ();
    my %SAMPLENUM = ();
    my @NTdata = ();
    my @NTwd = ();
    $sequence = 0;
    }

}

print "\n\n\n. . . . .DONE. . . . \n\n"