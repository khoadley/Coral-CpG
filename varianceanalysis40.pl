#!/usr/bin/perl
use strict;

# - - - - - H E A D E R - - - - - - - - - - - - - - - - -
#Averages data across a column and accounts for any instances of zero of N in the data-set.
# - - - - - U S E R    V A R I A B L E S - - - - - - - -

my $infolder = "07-KEGG-variance";
my $outfolder = "08-KEGG-variance";

# - - - - - G L O B A L  V A R I A B L E S  - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - M A I N - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\n Making data up . . .\n\n";
opendir(TOPIC, $infolder) or die "\n\nNADA $infolder, idiot . . . . \n\n";
my @folderlist = readdir(TOPIC);
foreach my $folder(@folderlist)
{	if ($folder =~ m/^\w/)
	{      print "averaging sequences in $folder\n";
	       system("mkdir $outfolder/$folder");
	       opendir(DATA, $infolder."/".$folder) or die "\n\nNADA $folder, idiot . . . . \n\n";
	       my @filelist = readdir(DATA);
	       foreach my $file(@filelist)
	       {	if ($file =~ m/^\w/)
			{	 my $sequences = 0;
				 open(OUT1, ">".$outfolder."/".$folder."/".$file);
				 my $infile = $infolder."/".$folder."/".$file;
                                 open(IN, "<$infile") or die "\n\nNADA $infile you FOOL!!!\n\n";
                                 my @DATA = <IN>; close(IN); shift(@DATA);	
                                 # 2. Parse sequence data . . . . . . . . . . . . .
                                 my $unid = 1000000000;                           # string to generate unique ids
                                 foreach my $row (@DATA)
                                 {     my $totalgenes = 1;
                                       my $totalsum = 0;
                                       my @rows = split(/ /, $row);
                                       foreach my $rows(@rows)
                                       {     if ($rows =~ m/\d/)
                                             {     $totalsum += $rows;
                                                   $totalgenes += 1;
                                             }
                                       }
                                       my $Nadjustedvalue = ($totalsum/$totalgenes);
                                       print OUT1 "$Nadjustedvalue\n";
                                 }
                                 close (OUT1);
                        }
               }
        }
                
}
#print "there are $totalsequences sequences in $outfile \n\n";
print "\n\n\n. . . . .DONE. . . . .\n\n\n";



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - S U B R O U T I N E S - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -






# - - - - - EOF - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - -
#
#
