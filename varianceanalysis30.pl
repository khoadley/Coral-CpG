#!/usr/bin/env perl
use strict;
# - - - - - H E A D E R - - - - - - - - - - - - - - - - -

# This script transposes each file in 06-KEGG-variance so that each column equals a single gene. T
# The script also averages each position and outputs it as the last column.
# All data is placed into the 07-KEGG-variance folder

# - - - - - U S E R    V A R I A B L E S - - - - - - - -
my $infolder = "06-KEGG-variance";
my $outfolder = "07-KEGG-variance";

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

my @rows = ();
my @transposed = ();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\n . . . Making data up . . .\n\n";

opendir(TOPIC, $infolder) or die "\n\nNADA $infolder, idiot . . . . \n\n";
my @folderlist = readdir(TOPIC);
foreach my $folder(@folderlist)
{	if ($folder =~ m/^\w/)
	{	print "Transposing sequences in folder $folder\n";
		system("mkdir $outfolder/$folder");
		opendir(DATA, $infolder."/".$folder) or die "\n\nNADA $folder, idiot . . . . \n\n";
		my @filelist = readdir(DATA);
		foreach my $file(@filelist)
		{	if ($file =~ m/^\w/)
			{	my $sequences = 0;
				open(OUT1, ">".$outfolder."/".$folder."/".$file);
				my $infile = $infolder."/".$folder."/".$file;
				open(FILENAME, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
				while (<FILENAME>)
				{	chomp;
					push @rows, [split / /];
                                        #print "@rows\n\n";
					$sequences += 1;
				}
				print "there are $sequences sequences in $file\n";
				for my $row (@rows)
				{       #print "$row\n\n";
                                        for my $column (0..$#{$row})
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
				@transposed = ();
			}
		}
	}
}

print "\n\n\n. . . . .DONE. . . . \n\n";