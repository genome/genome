package Genome::Model::Tools::BioSamtools::CompareBamFiles;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::BioSamtools::CompareBamFiles {
	is        => 'Genome::Model::Tools::BioSamtools',
	has_input => [
		first_bam_file => {
			is  => 'String',
			doc => 'First queryname sorted BAM file of (ELAND).',
		},
		second_bam_file => {
			is  => 'String',
			doc => 'Second query-name sorted BAM file (BWA) .',
		},
		merged_bam_file => {
			is => 'String',
			doc =>
			  'The path to output the resulting merged, unsorted BAM file.',
		},

	],
};

sub execute {
	my $self = shift;

	my $first_bam_file  = $self->first_bam_file;
	my $second_bam_file = $self->second_bam_file;
	my $merged_bam_file = $self->merged_bam_file;

	unless ( $first_bam_file && $second_bam_file && $merged_bam_file ) {
		die('Wrong input ');
	}

	my $merged_bam = Bio::DB::Bam->open( $merged_bam_file, 'w' );
	unless ($merged_bam) {
		die( 'Failed to open output BAM file: ' . $merged_bam_file );
	}

	my ( $first_bam, $first_header ) = validate_sort_order($first_bam_file);

	my ( $second_bam, $second_header ) = validate_sort_order($second_bam_file);
	$merged_bam->header_write($first_header);

	################### CHECKING READS#######
	while ( my $second_read = $second_bam->read1 ) {
		my $second_flag       = $second_read->flag;
		my $second_read_qname = $second_read->qname;

		#my $second_read_end   = 0;
		my $first_read = $first_bam->read1;
		my $first_flag = $first_read->flag;

		if ( $first_flag & 4 ) {
			if ( $first_flag & 64 && $second_flag & 64 ) {
				if ( $second_flag & 4 ) {

					#unaligned in both...do nothing
				}
				else {
					$merged_bam->write1($first_read);
				}
			}

			elsif ( $first_flag & 128 && $second_flag & 128 ) {
					
				if ( $second_flag & 4 ){
					
					#unaligned in both...do nothing
				}
					
				else{
					$merged_bam->write1($first_read);
				}
			}
			
		}
}
return 1;

sub validate_sort_order {
	my $bam_file = shift;
	my $bam      = Bio::DB::Bam->open($bam_file);
	unless ($bam) {
		die( 'Failed to open BAM file: ' . $bam_file );
	}
	my $header = $bam->header;
	my $text   = $header->text;
	my @lines  = split( "\n", $text );
	my @hds    = grep { $_ =~ /^\@HD/ } @lines;
	unless ( scalar(@hds) == 1 ) {
		die( 'Found multiple HD lines in header: ' . "\n\t"
			  . join( "\n\t", @hds ) )
		  . "\nRefusing to continue parsing BAM file: "
		  . $bam_file;
	}
	my $hd_line = $hds[0];
	if ( $hd_line =~ /SO:(\S+)/ ) {
		my $sort_order = $1;
		unless ( $sort_order eq 'queryname' ) {
			die(
'Input BAM files must be sorted by queryname!  BAM file found to be sorted by \''
				  . $sort_order
				  . '\' in BAM file: '
				  . $bam_file );
		}
	}
	else {
		die(
'Input BAM files must be sorted by queryname!  No sort order found for input BAM file: '
			  . $bam_file );
	}
	return ( $bam, $header );
}

}
1;
