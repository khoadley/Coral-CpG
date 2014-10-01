#!/usr/bin/perl
use strict;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#This script compiles a mean average CpG motif frequency for each KEGG category. It then calculates the variance from this average for each section of each gene for each KEGG category.
# - - - - - U S E R    V A R I A B L E S - - - - - - - -
my $infolder = "02-KEGG";
my $outfolder = "06-KEGG-variance";


# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

my $totalsequences = 0;
my @y;
my $means = 0;
my $meanscount = 1;
my $categories = "CpGmotif.txt";
my @KEGGcategories = ();
my %KEGGMEAN =();
my $KEGGsum = ();
my $KEGGcount = ();
my $foldnum = 0;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


print "\n\n Making data up . . .\n\n";
opendir(TOPIC, $infolder) or die "\n\nNADA $infolder, idiot . . . . \n\n";
my @folderlist = readdir(TOPIC);
foreach my $folder(@folderlist)
{	if ($folder =~ m/^\w/)
        {       push (@KEGGcategories, $folder);
                #print "computing means $categories for sequences in folder $folder\n";
                opendir(DATA, $infolder."/".$folder) or die "\n\nNADA $folder, idiot . . . . \n\n";
                my @filelist = readdir(DATA);
                foreach my $file(@filelist)
                {	if ($file =~ m/^($categories)/)
                        {       my $zerocount = 0;
                                my $sequences = 0;
                                my $infile = $infolder."/".$folder."/".$file;
                                open(IN, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
                                my @DATA = <IN>; close(IN); shift(@DATA);	
                                foreach my $entry (@DATA)
                                {       $meanscount = 1;
                                        $means = 0;
                                        my @genename = split(/ /, $entry);
                                        my $genehead = $genename[0];
                                        foreach my $mean (1..87)
                                        {   #print "$genename[$mean]\n";
                                            if ($genename[$mean] =~ m/\d/)
                                            {   #print "$genename[$mean]\t";
                                                $meanscount += 1;
                                                #print "$meanscount\t";
                                                $means += $genename[$mean];
                                                #print "$means\n";
                                            }
                                            else
                                            {
                                            }
                                        }
                                        my $genemean = ($means/$meanscount);
                                        if ($genemean == 0)
                                        {
                                        $zerocount += 1;
                                        }
                                        else
                                        {
                                        $KEGGsum += $genemean;
                                        $KEGGcount += 1;
                                        }
                                }
                        #print "there are $zerocount genes with a 0 $file for the KEGG cat $folder\n";
                        #print "there are $KEGGcount genes in folder $KEGGcategories[$foldnum] and the sum is $KEGGsum\n";
                        $KEGGMEAN{$KEGGcategories[$foldnum]} = $KEGGsum/$KEGGcount;
                        print "$KEGGcategories[$foldnum] equals $KEGGMEAN{$KEGGcategories[$foldnum]}\n\n";
                        $KEGGsum = 0;
                        $KEGGcount = 0;
                        $foldnum += 1;
                        }
                
                
                
                }
        }
}


#die;

opendir(TOPIC, $infolder) or die "\n\nNADA $infolder, idiot . . . . \n\n";
@folderlist = readdir(TOPIC);
foreach my $folder(@folderlist)
{	if ($folder =~ m/^\w/)
        {       system("mkdir $outfolder/$folder");
                open(OUT1, ">".$outfolder."/".$folder."/".$categories);
                push (@KEGGcategories, $folder);
                print "computing means $categories for sequences in folder $folder\n";
                #print OUT1 "$folder\t";
                opendir(DATA, $infolder."/".$folder) or die "\n\nNADA $folder, idiot . . . . \n\n";
                my @filelist = readdir(DATA);
                foreach my $file(@filelist)
                {	if ($file =~ m/^($categories)/)
                        {       my $zerocount = 0;
                                my $sequences = 0;
                                my $infile = $infolder."/".$folder."/".$file;
                                open(IN, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
                                my @DATA = <IN>; close(IN); shift(@DATA);	
                                foreach my $entry (@DATA)
                                {       $meanscount = 1;
                                        $means = 0;
                                        my @genename = split(/ /, $entry);
                                        my $genehead = $genename[0];
                                        print OUT1 "$genehead ";
                                        foreach my $mean (1..87)
                                        {   #print "$genename[$mean]\n";
                                            if ($genename[$mean] =~ m/\d/)
                                            {   my $meanie = $genename[$mean];
                                                #print "$meanie\t";
                                                foreach my $j(@KEGGcategories)
                                                {       if ($j eq $folder)
                                                        {   my $catmean = $KEGGMEAN{$KEGGcategories[$j]};
                                                            #print "$KEGGcategories[$j] equals $catmean\n";
                                                            my  $vary = ((($meanie-$catmean)*($meanie-$catmean))/$catmean);
                                                            #print "$vary\n";
                                                            print OUT1 "$vary ";
                                                            last;
                                                        }
                                                }
                                            }
                                            else
                                            {
                                            print OUT1 "- ";
                                            }
                                        }
                                        print OUT1 "\n";
                                }
                        }
             
                }
        }
}
print "\n\n. . . . . . . . . . . . .DONE!!!!(.)(.)( o )( o ). . . . . \n\n\n";