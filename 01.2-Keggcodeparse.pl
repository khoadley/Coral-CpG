#!/usr/bin/perl
use strict;


my $primaryheirarchy;
my $secondaryheirarchy;
my $tertiaryheirarchy;
 
 
my $outfile = "kegglist.txt";

open(OUT, ">$outfile");
print OUT "Kegg-code\tKO-Annotation\tTertiaryHeirarchy\tSecondaryHeirarchy\tPrimaryHeirarchy\n";
 
open (FILE, '01-keggtest.txt') or die ("you suck at this");
my @FILE = <FILE>;
        foreach  my $line (@FILE)
            {my @values = split(' ', $line, 3);
                if ($values[0] eq 'A')
                        {$primaryheirarchy = $values[1]}
                if ($values[0] eq 'B')
                        {$secondaryheirarchy = $values[1]}
                if ($values[0] eq 'C')
                        {$tertiaryheirarchy = $values[2]}
                if ($values[0] eq 'D')
                        {print OUT "$values[1]\t$values[2]\t$tertiaryheirarchy\t$secondaryheirarchy\t$primaryheirarchy\n";}
                        
            }

close OUT;


print "\n\n. . . . . . . DONE. . . . . . . \n\n";
                             
                                
                                
                                
 


 
 


 
 
 
