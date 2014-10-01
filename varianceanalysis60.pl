#!/usr/bin/env perl
use strict;
# - - - - - H E A D E R - - - - - - - - - - - - - - - - -

# Transposes all data from 09-KEGG-variance so that each column is a kegg category

# - - - - - U S E R    V A R I A B L E S - - - - - - - -
my $infolder = "09-KEGG-variance";
my $outfolder = "10-KEGG-variance";

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

my @rows = ();
my @transposed = ();
my $count = 0;

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\n . . . Making data up . . .\n\n";

opendir(DATA, $infolder) or die "\n\nNADA $infolder, idiot . . . . \n\n";
my @filelist = readdir(DATA);
foreach my $file(@filelist)
{   if ($file =~ m/^\w/)
    {	my $sequences = 0;
        open(OUT1, ">".$outfolder."/".$file);
        my $infile = $infolder."/".$file;
        open(FILENAME, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
        while (<FILENAME>)
        {	chomp;
                push @rows, [split / /];
                $sequences += 1;
        }
        print "there are $sequences sequences in $file\n";
        for my $row (@rows)
        {	for my $column (0..$#{$row})
                {   	push(@{$transposed[$column]}, $row->[$column]);
                }
        }
        for my $new_row (@transposed)
        {   	for my $new_col (@{$new_row})
                {   	print OUT1 $new_col, " ";
                }
                print OUT1 "\n";
        }
        close (OUT1);
        @rows = ();
        @transposed = ()
    }
}

print "\n\n\n. . . . .DONE. . . . \n\n"