package Genome::Model::SmallRna::Command::FilterNewBam;

use strict;
use warnings;

use Genome;

class Genome::Model::SmallRna::Command::FilterNewBam {
	is        => 'Genome::Model::Tools::BioSamtools',
	has_input => [
		bam_file => {
			is  => 'String',
			doc => 'Input BAM file of alignments.',
		},
		filtered_bam_file => {
			is  => 'String',
			doc => 'Output BAM file of filtered read alignments .',
			is_output =>1
		},
		
		read_size_bin => {
			is => 'Text',
			doc =>'Min_max read length bin separated by \'_\' : eg 17_75',
		},
		
		xa_tag => {
			is  => 'Boolean',
			doc => 'Remove alignments where Map Score is 0 but no XA tag is reported',
            is_optional => 1,
			
		},
	],
};

sub execute {
	my $self = shift;

	my $bam_file 	= $self->bam_file;
	my $new_bam	 	= $self->filtered_bam_file;
	
	
	
    
	my $in_bam = Bio::DB::Bam->open($bam_file);
	unless ($in_bam) {
		die( 'Failed to open BAM file: ' . $bam_file );
	}
	
	my $out_bam = Bio::DB::Bam->open( $new_bam, 'w' );
	unless ($out_bam) {
		die( 'Failed to open output BAM file: ' . $new_bam );
	}
	
	my $header = $in_bam->header;	
	$out_bam->header_write($header);
	
	while ( my $in_read = $in_bam->read1 ) {
	
    	my $read_length = $in_read->l_qseq;
		my $in_flag 	= $in_read->flag;
		my $map_score   = $in_read->qual;	
    	
    	if (defined($self->read_size_bin)) {
    		my $bin 			= $self->read_size_bin;

    		my ($min,$max) = split('_',$bin);
    		
			if ( $read_length > ($min-1) && $read_length < ($max + 1) && $in_flag != 4) 
				{
				if (defined($self->xa_tag)) {				
					if ( $map_score eq '0')
					  {
					  	if ($in_read->aux_get("XA") )  
					  	{
					  		$out_bam->write1($in_read);
					  	} 
					  }
					else 
					{
						$out_bam->write1($in_read);
					}
				}
				else
				{
					$out_bam->write1($in_read);
				}
		}
	    

		}
	    
	}   
	return 1;
}

1;

__END__