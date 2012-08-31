#!/gsc/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
#my $parent_dir='/gscmnt/sata141/techd/jhundal/Schreiber_mouse/SchreiberProject/';


my $in  = Bio::SeqIO->new(-file => $ARGV[0] ,
                      -format => 'Fasta');
                           
my $out = Bio::SeqIO->new(-file => ">".$ARGV[1] ,
                        -format => 'Fasta');

while ( my $seq = $in->next_seq() ) {
	my $seq_string = $seq->seq; 
	   if ( $seq_string =~ /[^A-Z]/ )
	   {
	   	print $seq_string."\n";
	   	
	   }
       else 
       {
       	$out->write_seq($seq);
       }   
    }
    
