#!/usr/bin/perl
use strict;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
# Places all data from 08-KEGG-variance into a single file, as a column per kegg category
# - - - - - U S E R    V A R I A B L E S - - - - - - - -
my $infolder = "08-KEGG-variance";
my $outfolder = "09-KEGG-variance";


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

my $categories = 'CpGmotif.txt';
my @datum = ();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


print "\n\n Making data up . . .\n\n";


open(OUT1, ">".$outfolder."/".$categories);
opendir(TOPIC, $infolder) or die "\n\nNADA $infolder, idiot . . . . \n\n";
my @folderlist = readdir(TOPIC);
foreach my $folder(@folderlist)
{	if ($folder =~ m/^\w/)
        {	print OUT1 "$folder ";
                opendir(DATA, $infolder."/".$folder) or die "\n\nNADA $folder, idiot . . . . \n\n";
                my @filelist = readdir(DATA);
                foreach my $file(@filelist)
                {	if ($file =~ m/^($categories)/)
                        {       my $infile = $infolder."/".$folder."/".$file;
                                open(IN, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
                                my @DATA = <IN>; close(IN); shift(@DATA);	
                                foreach my $entry (@DATA)
                                {   chomp $entry;
                                    print $entry;
                                    print OUT1 "$entry ";    
                                }       
                        print OUT1 "\n";
                        }
                }
                
        }
}


print "\n\n. . . . . . . . . . . . .DONE!!!!(.)(.)( o )( o ). . . . . \n\n\n";