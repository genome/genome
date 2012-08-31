package Genome::Model::Tools::Synthesizer::SizeBam;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Synthesizer::SizeBam {
	is        => 'Genome::Model::Tools::Synthesizer::Base',
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

	],
};

sub help_brief {
"Run the Synthesizer SizeBam module to fractionate and filter bam based on specifed read lengths "
	  ,;
}

sub help_detail {
"SizeBam module of Synthesizer filters and splits the main alignment BAM into smaller size fractions based on read-lengths.Use this module to filter based on size of species. For longer reads upto 75bp, we recommend using the following bins: 17_25,26_35,36_60,61_70,71_75,17_75. This module also calculates samtools flagstat metrics for each of the filtered & sized bams";
}


sub execute {
	my $self = shift;

	my $bam_file 	= $self->bam_file;
	my $filtered_bam	 	= $self->filtered_bam_file;
	

############## FILTERING AND SIZE-SELECTING BAM ########################### 
	
if ( defined( $self->read_size_bin ) ) {
		my $bin = $self->read_size_bin;
		my $sized_bam = filter_bam ( $bam_file, $filtered_bam,$bin);

	
###################### INDEXING AND CALCULATING FLAGSTAT ##################	
	unless (-s $sized_bam) {
		die( 'Failed to open filtered BAM file: ' . $sized_bam );
	}
	
	
	
	my $index_file= $sized_bam.'.bai';
	my $index_cmd = 'samtools index ' . $sized_bam .' > '. $index_file;
	
	Genome::Sys->shellcmd(
        cmd => $index_cmd,
        input_files => [$sized_bam],
        output_files => [$index_file],
    );    
    

	my $flagstat_file= $sized_bam.'.flagstat';	
    my $cmd = 'samtools flagstat '. $sized_bam .' > '. $flagstat_file ;
    
    Genome::Sys->shellcmd(
        cmd => $cmd,
        input_files => [$sized_bam],
        output_files => [$flagstat_file],
    );  
}
    
    else {
				print  "Please define a min_max size bin to filter BAM ; it should be <minsize>_<maxsize>";
	}  
    
    return 1;
}

sub filter_bam {
	
	my ($bam_file, $new_bam,$bin) = @_;
	
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

		if ( $bin =~ m/^(\d+)(_)(\d+)/ ) {
			my ( $min, $max ) = split( '_', $bin );
			while ( my $in_read = $in_bam->read1 ) {
				my $read_length = $in_read->l_qseq;
				my $in_flag     = $in_read->flag;
				my $map_score   = $in_read->qual;

				if (   $read_length > ( $min - 1 ) && $read_length < ( $max + 1 ) && $in_flag != 4 ){
					if ( $map_score eq '0' ) {
						if ( $in_read->aux_get("XA") ) {
							$out_bam->write1($in_read);
						}
					}
					else {
						$out_bam->write1($in_read);}
				}
			}
		return $new_bam;	
		}
		else{
				print  "Please check format of the min_max size bin ; it should be <minsize>_<maxsize>";
			}  	
		
		
		
};


sub parse_file_into_hashref {   # This subroutine is taken from Genome::Model::Tools::Sam::Flagstat
    my ($class, $flag_file) = @_;

    unless ($flag_file and -s $flag_file) {
        warn "Bam flagstat file: $flag_file is not valid";
        return;
    }

    my $flag_fh = Genome::Sys->open_file_for_reading($flag_file);
    unless($flag_fh) {
        warn 'Fail to open ' . $flag_file . ' for reading';
        return;
    }

    my %data;
    my @lines = <$flag_fh>;
    $flag_fh->close;
    my $line_ct = scalar @lines;

    while ($line_ct and $lines[0] =~ /^\[.*\]/){
        my $error = shift @lines;
        if ($error =~ /EOF marker is absent/){
            $line_ct--;
        }
        push @{ $data{errors} }, $error;
    }

    unless ($line_ct =~ /^1[12]$/) {#samtools 0.1.15 (r949), older versions get 12 lines, newer get 11 line with QC passed/failed separate in each line
        warn 'Unexpected output from flagstat. Check ' . $flag_file;
        return;
    }

   for (@lines) {
        chomp;

        if (/^(\d+) \+ (\d+) in total/) {
            $data{reads_marked_passing_qc} = $1;
            $data{reads_marked_failing_qc} = $2;
            $data{total_reads} = $data{reads_marked_passing_qc} + $data{reads_marked_failing_qc};
        }

        if (/^(\d+) in total/) {
            $data{total_reads} = $1;
        }

        #For older samtools (before r949) total reads count in flagstat includes QC failed read count
        if (/^(\d+) QC failure/) { 
            $data{reads_marked_failing_qc} = $1;
            $data{reads_marked_passing_qc} = $data{total_reads} - $data{reads_marked_failing_qc};
        }

        $data{reads_marked_duplicates} = $1 if /^(\d+) (\+\s\d+\s)?duplicates$/;
        ($data{reads_mapped}, $data{reads_mapped_percentage}) = ($1, $3)
            if /^(\d+) (\+\s\d+\s)?mapped \((\d{1,3}\.\d{2}|nan)\%[\:\)]/;
        undef($data{reads_mapped_percentage})
            if $data{reads_mapped_percentage} && $data{reads_mapped_percentage} eq 'nan';

        $data{reads_paired_in_sequencing} = $1 if /^(\d+) (\+\s\d+\s)?paired in sequencing$/;
        $data{reads_marked_as_read1}      = $1 if /^(\d+) (\+\s\d+\s)?read1$/;
        $data{reads_marked_as_read2}      = $1 if /^(\d+) (\+\s\d+\s)?read2$/;

        ($data{reads_mapped_in_proper_pairs}, $data{reads_mapped_in_proper_pairs_percentage}) = ($1, $3)
            if /^(\d+) (\+\s\d+\s)?properly paired \((\d{1,3}\.\d{2}|nan)\%[\:\)]/;
        undef($data{reads_mapped_in_proper_pairs_percentage})
            if $data{reads_mapped_in_proper_pairs_percentage} && $data{reads_mapped_in_proper_pairs_percentage} eq 'nan';

        $data{reads_mapped_in_pair} = $1 if /^(\d+) (\+\s\d+\s)?with itself and mate mapped$/;

        ($data{reads_mapped_as_singleton}, $data{reads_mapped_as_singleton_percentage}) = ($1, $3)
            if /^(\d+) (\+\s\d+\s)?singletons \((\d{1,3}\.\d{2}|nan)\%[\:\)]/;
        undef($data{reads_mapped_as_singleton_percentage})
            if $data{reads_mapped_as_singleton_percentage} && $data{reads_mapped_as_singleton_percentage} eq 'nan';

        $data{reads_mapped_in_interchromosomal_pairs}    = $1 if /^(\d+) (\+\s\d+\s)?with mate mapped to a different chr$/;
        $data{hq_reads_mapped_in_interchromosomal_pairs} = $1 if /^(\d+) (\+\s\d+\s)?with mate mapped to a different chr \(mapQ>=5\)$/;
    }

    return \%data;
}

	
1;

__END__
