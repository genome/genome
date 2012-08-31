package Genome::Model::Tools::Bwa::AlignReadsMulti;

use strict;
use warnings;

use Genome;
use Genome::Sys;

class Genome::Model::Tools::Bwa::AlignReadsMulti {
	is  => 'Genome::Model::Tools::Bwa',
	has => [
		#######################################################
		dna_type => {
			doc =>
'Optional switch which can be "dna" or "rna".  Each choice causes the application to use a specific primer file.  If no dna_type value is provided, an adaptor_file parameter must be used.',
			is          => 'String',
			is_optional => 1,
		},
		align_options => {
			doc =>
'The bwa aln parameters. These should be specified in a quoted string with single dashes, e.g. "-x -y -z"',
			is            => 'String',
			is_optional   => 1,
			default_value => '',
		},
		upper_bound => {
			doc =>
'The bwa sampe option (-a) designating the upper bound on insert size. Defaults to 600.',
			is            => 'Integer',
			is_optional   => 1,
			default_value => 600,
		},

		top_hits => {
			doc =>
'The bwa samse option (-n) designating how many top hits to output.  Default value -1 to disable outputting multiple hits',
			is            => 'Integer',
			is_optional   => 1,
			default_value => -1,
		},
		max_occurrences => {
			doc =>
'The bwa sampe option (-o) describing maximum occurrences of a read for pairing.  A read with more than -o occurrences will be treated as a single-end read.  Reducing helps faster pairing.  Default 100000.',
			is            => 'Integer',
			is_optional   => 1,
			default_value => 100000,
		},
		force_fragments => {
			doc           => 'Optional switch to force fragment processing.',
			is            => 'Integer',
			is_optional   => 1,
			default_value => 0,
		},
		#####################################################
		#input files
		ref_seq_file => {
			doc =>
'Required input file name containing the reference sequence file.',
			is => 'String',
		},
		files_to_align_path => {
			doc =>
'Path to a directory or a file or a pipe separated list of files containing the reads to be aligned.  Must be in fastq format.',
			is => 'String',
		},
		#####################################################
		#output files
		aligner_output_file => {
			doc => 'Optional output log file containing results of the run.',
			is  => 'String',
		},
		concise_file => {
			doc         => 'Output file containing the aligned map data.',
			is          => 'String',
			is_optional => 0,
		},
		unaligned_reads_file => {
			doc => 'Output file containing unaligned data.',
			is  => 'String',
		},
		duplicate_mismatch_file => {
			doc =>
'Output file containing dumped duplicate mismatches specified by the (-H) bwa parameter. There is no default value.  If this file is not specified duplicate mismatches will not be dumped',
			is          => 'String',
			is_optional => 1,
		},
		temp_directory => {
			doc =>
'Optional temp directory where fastq and bfqs will be stored when generated.  If no temp directory is specified, a temporary directory in /tmp will be created and then removed.',
			is          => 'String',
			is_optional => 1,
		},
		#####################################################
		#private variables
		_files_to_align_list => {
			doc         => 'The list of input files to align.',
			is          => 'List',
			is_optional => 1,
		},
		_ref_seq_index_file => {
			doc         => 'Samtools reference sequence index file.',
			is          => 'String',
			is_optional => 1,
		},
	],
};

sub help_synopsis {
	return <<EOS
    A BWA based utility for aligning reads.;
EOS
}

sub help_brief {
	return <<EOS
    A BWA based utility for aligning reads.;
EOS
}

sub help_detail {
	return <<EOS
Provides an interface to the BWA aligner.  Inputs are:

'ref-seq-file' - The reference sequence file which to align reads to.  Specified by a path to a fasta file. 

'files-to-align-path' - The file or set of files which contain the read fragments which are going to be aligned to the reference sequence.  The path can be a single file, a pipe seperated list of two files for paired end reads, or a directory containing one or two files.  These files must be in the fastq format. 


EOS
}

sub create {
	my $class = shift;
	my $self  = $class->SUPER::create(@_);

	unless ($self) {
		return;
	}

	unless ( $self->use_version ) {
		my $msg = 'use_version is a required parameter to ' . $class;
		$self->delete;
		die($msg);
	}
	unless ( $self->bwa_path ) {
		my $msg =
		    'No path found for bwa version '
		  . $self->use_version
		  . ".  Available versions are:\n";
		$msg .= join( "\n", $self->available_bwa_versions );
		$self->delete;
		die($msg);
	}

	$self->_ref_seq_index_file( $self->ref_seq_file . ".fai" );

	unless ( -e $self->_ref_seq_index_file ) {
		$self->error_message( "Samtools refseq list does not exist.  Expected "
			  . $self->_ref_seq_idx_file
			  . ".  Use 'samtools faidx' to go make one!" );
		return;
	}

	my @listing;
	my $dir_flag = 0;
	my @pipe_list;

	#check to see if files to align path is a pipe delimited list of files
	$self->status_message( "Files to align: " . $self->files_to_align_path );
	my $pipe_char_index = index( $self->files_to_align_path, '|' );

	#$self->status_message("Comma index: ".$pipe_char);
	if ( $pipe_char_index > -1 ) {
		@pipe_list = split( /\|/, $self->files_to_align_path );
		for my $pipe_file (@pipe_list) {

			#make sure each file exists
			if ( -f $pipe_file ) {
				push @listing, $pipe_file;
			}
			else {
				$self->error_message( 'File does not exist: ' . $pipe_file );
			}
		}
	}
	else {
		$self->status_message('Checking existence of reads file: >'.$self->files_to_align_path."<");
		$self->error_message('Checking existence of reads file: >'.$self->files_to_align_path."<");
		push @listing, $self->files_to_align_path;
		#not a pipe delimited list
		#check to see if files to align path is a dir or file
		if ( -e $self->files_to_align_path ) {
			#$self->status_message('Path is a file');
			#push @listing, $self->files_to_align_path;
		}
		elsif ( -d $self->files_to_align_path ) {

			#$self->status_message('Path is dir.');
			@listing = glob( $self->files_to_align_path . '/*' );
		}
		else {
			$self->error_message('Input file does not exist.');
		}
	}

	if ( @listing > 2 || @listing == 0 ) {
		$self->error_message(
			"Should provide either 1 or 2 input files for alignment.");
		return;
	}

	$self->_files_to_align_list( \@listing );

	return $self;
}

sub _check_output_files_and_log_files_for_success {
	my ( $self, $src_name, $log_regex, $check_nonzero, @files_to_check ) = @_;

	for (@files_to_check) {
		if ( !-e $_ ) {
			$self->error_message("Can't find expected output filename $_");
			return;
		}
	}

	my $aligner_output_fh = IO::File->new( $self->aligner_output_file );
	unless ($aligner_output_fh) {
		$self->error_message(
			    "Can't open expected aligner log file to verify completion "
			  . $self->aligner_output_file
			  . "$!" );
		return;
	}

	while (<$aligner_output_fh>) {
		if (m/\[($src_name)\].*?$log_regex/) {
			$aligner_output_fh->close();
			if ( !$check_nonzero || $2 > 0 ) {
				return 1;
			}
			return;
		}
	}
}

sub concise_to_sam {
	my $self          = shift;
	my $file_in       = shift;
	my $file_out      = shift;
	my $unaligned_out = shift;

	my $fh_in  = Genome::Sys->open_file_for_reading($file_in);
	my $fh_out = Genome::Sys->open_file_for_writing($file_out);
	my $fh_unaligned_out =
	  Genome::Sys->open_file_for_writing($unaligned_out);
	my $done = 0;
	my $line = $fh_in->getline;
	my $read_name;

#my $unaligned_block = "4\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*";
	my $unaligned_block =
"4\t*\t0\t0\t*\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB";
	my $aligned_block_example =
"c15473140  0   gi|11612676|gb|AE000513.1|  315532  37  100M    *   0   0   TTCNTCCCCCGTTCGGCCAGCGAGCCGGTCGCCAAGGTGCTCAACAGCGCCAAGGCCAACGCCCTGCACAACGACGAGATGCTCGAAGATCGCCTGTTCG    *   XT:A:U  NM:i:";

#my $aligned_block1 = "255\t100M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t*\tXT:A:U\tNM:i:";
	my $aligned_block1 =
"255\t100M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\tBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB\tXT:A:U\tNM:i:";
	my $aligned_block2 = "X0:i:1\tX1:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tMD:Z:0";

	while ( !$done ) {

		my $first_char = substr( $line, 0, 1 );
		if ( $first_char eq '>' ) {

			#$line =~ m/>(\w*)\s.*$/;
			$line =~ m/>(\w*)\s(\d*).*$/;
			$read_name = $1;
			my $unaligned = $2;
			if ( $unaligned == '0' ) {
				print $fh_out "$read_name\t$unaligned_block\n";
				print $fh_unaligned_out "$read_name\t$unaligned_block\n";
			}
		}
		else {

			#$line =~ m/^(\w*)\s*[-+]\d*\s*\d*$/;
			$line =~ m/^(\S*)\s*[-+](\d*)\s*(\d)*$/;
			print $fh_out
			  "$read_name\t0\t$1\t$2\t$aligned_block1$3\t$aligned_block2\n";
		}

		$line = $fh_in->getline;
		if ( !defined($line) ) {
			$done = 1;
		}

	}

	$fh_in->close;
	$fh_out->close;

	return 1;
}

sub execute {

	my $self = shift;
	$self->dump_status_messages(1);

	$self->status_message("\n");
	$self->status_message('Running AlignReadsMulti with parameters');
	$self->status_message('-----------------------------------');
	$self->status_message('Input Files:');
	$self->status_message( 'Reference sequence file:' . $self->ref_seq_file );
	$self->status_message('Files to align path:' . $self->files_to_align_path );
	$self->status_message('Files to align list:' . join( ",", @{$self->_files_to_align_list} ) );
	$self->status_message('');
	$self->status_message('Output Files:');
	$self->status_message( 'Concise file:' . $self->concise_file );
	$self->status_message(
		'Unaligned reads file:' . $self->unaligned_reads_file )
	  if defined( $self->unaligned_reads_file );
	$self->status_message(
		'Aligner output messages:' . $self->aligner_output_file )
	  if defined( $self->aligner_output_file );
	$self->status_message('');
	$self->status_message('Other Parameters:');
	$self->status_message( 'Align options:' . $self->align_options )
	  if defined( $self->align_options );
	$self->status_message( 'Upper bound value:' . $self->upper_bound )
	  if defined( $self->upper_bound );
	$self->status_message( 'BWA version:' . $self->use_version )
	  if defined( $self->use_version );
	$self->status_message("\n");

	my $tmp_dir;

	if ( defined( $self->temp_directory ) ) {
		$tmp_dir = Genome::Sys->create_directory(
			$self->temp_directory );
	}
	else {
		$tmp_dir = File::Temp::tempdir( CLEANUP => 1 );
	}

	#get the files to align ready
	my @input_file_list = @{ $self->_files_to_align_list };

	#if there is more than one input file, cat them together
	if ( scalar(@input_file_list) > 1 ) {
		$self->status_message("Combining reads files.");
		my $combined_file = "$tmp_dir/combined_read_files.txt";
		$self->status_message( "Combined file: " . $combined_file );
		my $rv_cat = Genome::Sys->cat(
			input_files => \@input_file_list,
			output_file => $combined_file
		);
		
		unless ($rv_cat) {
			$self->error_message(
				"Failed AlignReadsMulti on cat of input reads files.");
		}

		@input_file_list = ($combined_file);
		$self->status_message(
			"Listing contains: " . join( " ", @input_file_list ) );
	}

	### STEP 1: Use "bwa aln" to generate the alignment coordinate file for the input reads (.sai file)
	### Must be run once for each of the input files

	my @sai_intermediate_files;
	my $input = $input_file_list[0];

		#testing
		#my $tmpfile = File::Temp->new( DIR => $tmp_dir, SUFFIX => ".sai", CLEANUP=>0 );
		my $tmpfile = $tmp_dir . "/intermediate.sai";

		my $cmdline = $self->bwa_path
		  . sprintf( ' aln %s %s %s 1> ',
			$self->align_options, $self->ref_seq_file, $input )

		  #. $tmpfile->filename . ' 2>>'
		  . $tmpfile . ' 2>>' . $self->aligner_output_file;

		push @sai_intermediate_files, $tmpfile;

		# run the aligner
		Genome::Sys->shellcmd(
			cmd         => $cmdline,
			input_files => [ $self->ref_seq_file, $input ],

			#output_files => [ $tmpfile->filename ],
			output_files              => [$tmpfile],
			skip_if_output_is_present => 0,
		);


	#### STEP 2: Use "bwa samse" or "bwa sampe" to perform single-ended or paired alignments, respectively.
	#### Runs once for ALL input files

	#my $concise_file = $tmp_dir . "/concise.txt";

	my $sam_command_line = "";

	# single ended run
	my $top_hits_option = '-n ' . $self->top_hits;

	$sam_command_line = $self->bwa_path . sprintf(
		' samse %s %s %s %s',
		$top_hits_option, $self->ref_seq_file,

		#$sai_intermediate_files[0]->filename,
		$sai_intermediate_files[0],
		$input_file_list[0]
	  )
	  . " 1>"
	  . $self->concise_file. " 2>>"
	  . $self->aligner_output_file;

	$self->status_message("Running samXe to get the output alignments");
	$self->status_message("samXe command: $sam_command_line");

	my $rv_samse =  Genome::Sys->shellcmd( cmd => $sam_command_line );

	#my $multi_hit_sam_file = $tmp_dir."/multi.txt";
	#my $multi_hit_sam_file = $self->alignment_file;
	#my $unaligned_sam_file = $tmp_dir . "/unaligned.txt";
	#my $rv_concise = $self->concise_to_sam( $concise_file, $multi_hit_sam_file, $unaligned_sam_file );

	return 1;

}

1;
