package Genome::Model::Tools::Bwa::AlignReads;

use strict;
use warnings;

use Genome;
use Workflow;
use Genome::Sys;

class Genome::Model::Tools::Bwa::AlignReads {
    is  => 'Genome::Model::Tools::Bwa',
    has_param => [
            lsf_resource => 
            {
                value => "-R 'select[model!=Opteron250 && type==LINUX64] span[hosts=1] rusage[mem=10000]' -M 10000000",
            },
    ],
    has_input => [
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
            is_input    => 1,
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
            is_input    => 1,
        },
        files_to_align_path => {
            doc =>
'Path to a directory or a file or a pipe separated list of files containing the reads to be aligned.  Must be in fastq format.',
            is => 'String',
            is_input    => 1,
        },
        #####################################################
        #output files
        aligner_output_file => {
            doc => 'Optional output log file containing results of the run.',
            is  => 'String',
            is_input    => 1,
        },
        alignment_file => {
            doc         => 'Output file containing the aligned map data.',
            is          => 'String',
            is_optional => 0,
            is_input    => 1,
        },
        unaligned_reads_file => {
            doc => 'Output file containing unaligned data.',
            is  => 'String',
            is_input    => 1,
            is_output => 1,
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
        picard_conversion => {
            doc => 'Flag specifying whether to use the samtools conversion tools or picard conversion tools from sam to bam.  Default is samtools. ',
            is          => 'Integer',
            is_optional => 1,
            default_value => 0,
        },
        sam_only => {
            doc => 'Flag to leave the sam files in tact.',
            is          => 'Integer',
            is_optional => 1,
            default_value => 0,
        },
        read_group_tag => {
            doc => 'Flag to leave the sam files in tact.',
            is  => 'String',
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
              . $self->_ref_seq_index_file
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

        #not a pipe delimited list
        #check to see if files to align path is a dir or file
        if ( -f $self->files_to_align_path ) {

            #$self->status_message('Path is a file');
            push @listing, $self->files_to_align_path;
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
    }


    $self->_files_to_align_list( \@listing );

    return $self;
}

sub _check_output_files_and_log_files_for_success {
	my ($self, $src_name, $log_regex, $check_nonzero, @files_to_check) = @_;

	for (@files_to_check) {
		if (!-e $_) {
			$self->error_message("Can't find expected output filename $_");
			return;
		}
	}

    my $aligner_output_fh = IO::File->new($self->aligner_output_file);
    unless ($aligner_output_fh) {
        $self->error_message("Can't open expected aligner log file to verify completion " . $self->aligner_output_file . "$!"
        );
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

sub execute {

    my $self = shift;
    $self->dump_status_messages(1);

    $self->status_message("\n");
    $self->status_message('Running AlignReads with parameters');
    $self->status_message('-----------------------------------');
    $self->status_message('Input Files:');
    $self->status_message( 'Reference sequence file:' . $self->ref_seq_file );
    $self->status_message(
        'Files to align path:' . $self->files_to_align_path );
    $self->status_message(
        'Files to align list:' . $self->_files_to_align_list );
    $self->status_message('');
    $self->status_message('Output Files:');
    $self->status_message( 'Alignment file:' . $self->alignment_file );
    $self->status_message(
        'Unaligned reads file:' . $self->unaligned_reads_file )
      if defined( $self->unaligned_reads_file );
    $self->status_message(
        'Aligner output messages:' . $self->aligner_output_file )
      if defined( $self->aligner_output_file );
    $self->status_message('');
    $self->status_message('Other Parameters:');
    $self->status_message('Sam only flag:'. $self->sam_only );
    $self->status_message( 'Align options:' . $self->align_options )
      if defined( $self->align_options );
    $self->status_message( 'Upper bound value:' . $self->upper_bound )
      if defined( $self->upper_bound );
    $self->status_message( 'BWA version:' . $self->use_version )
      if defined( $self->use_version );
    $self->status_message("\n");

    print "Align Reads begun with " . $self->unaligned_reads_file . "\n";

    my $tmp_dir;

    if ( defined( $self->temp_directory ) ) {
        $tmp_dir =
          Genome::Sys->create_directory(
            $self->temp_directory );
    }
    else {
        $tmp_dir = File::Temp::tempdir( CLEANUP => 1 );
    }

    #get the files to align ready
    my @input_file_list = @{ $self->_files_to_align_list };
    my $files_to_align  = join( ' ', @input_file_list );
    
    unless ( -e $self->unaligned_reads_file && -e $self->alignment_file ) {


	## TODO :: clean this mess up
	my $force_frag_file;
        if ($self->force_fragments) {
                $self->status_message("Forcing fragments.");
                $force_frag_file = "$tmp_dir/force-frag";
                $self->status_message("Frag file: ".$force_frag_file);
                my $cmd = "cat ".join(" ",@input_file_list)." > ".$force_frag_file;
                $self->status_message("Cat command: ".$cmd);
                my @result = `$cmd`;
                #clear the listing and replace it with the new combined file name for processing
                @input_file_list = ($force_frag_file);
                $self->status_message("Listing contains: ".join(" ",@input_file_list) );
        }


    ### STEP 1: Use "bwa aln" to generate the alignment coordinate file for the input reads (.sai file)
    ### Must be run once for each of the input files

    my @sai_intermediate_files;

    foreach my $input (@input_file_list) {


        my $tmpfile = File::Temp->new( DIR => $tmp_dir, SUFFIX => ".sai" );

        my $cmdline = $self->bwa_path
          . sprintf( ' aln %s %s %s 1> ',
            $self->align_options, $self->ref_seq_file, $input )
          . $tmpfile->filename . ' 2>>'
          . $self->aligner_output_file;


        push @sai_intermediate_files, $tmpfile;

        # run the aligner
        Genome::Sys->shellcmd(
            cmd          => $cmdline,
            input_files  => [ $self->ref_seq_file, $input ],
            output_files => [ $tmpfile->filename ],
            skip_if_output_is_present => 0,
        );
        unless ($self->_check_output_files_and_log_files_for_success("bwa_aln_core", '(\d+) sequences have been processed', 1, $tmpfile->filename))
        {
            return;
        }
    }

    #### STEP 2: Use "bwa samse" or "bwa sampe" to perform single-ended or paired alignments, respectively.
    #### Runs once for ALL input files


        my $sam_command_line = "";
        if ( @input_file_list == 1 ) {

            # single ended run
            my $top_hits_option = '-n ' . $self->top_hits;

            $sam_command_line = $self->bwa_path
              . sprintf(
                ' samse %s %s %s %s',
                $top_hits_option, $self->ref_seq_file,
                $sai_intermediate_files[0]->filename,
                $input_file_list[0]
              )
              . " 2>>"
              . $self->aligner_output_file;

        }
        else {

            # paired run
            my $upper_bound_option     = '-a ' . $self->upper_bound;
            my $max_occurrences_option = '-o ' . $self->max_occurrences;
            my $paired_options = "$upper_bound_option $max_occurrences_option";



            $sam_command_line = $self->bwa_path
              . sprintf(
                ' sampe %s %s %s %s',
                $paired_options,
                $self->ref_seq_file,
		join (' ', map {$_->filename} @sai_intermediate_files),
                join (' ', @input_file_list)
              )
              . " 2>>"
              . $self->aligner_output_file;
        }

	$self->status_message("Running samXe to get the output alignments");
	$self->status_message("Command: $sam_command_line");
       #BWA is not nice enough to give us an unaligned output file so we need to
       #filter it out on our own

        my $sam_map_output_fh =
          File::Temp->new( DIR => $tmp_dir, SUFFIX => ".sam" );
        my $unaligned_output_fh =
          IO::File->new( ">" . $self->unaligned_reads_file );

        my $sam_run_output_fh = IO::File->new( $sam_command_line . "|" );
        if ( !$sam_run_output_fh ) {
            $self->error_message("Error running $sam_command_line $!");
            return;
        }

        while (<$sam_run_output_fh>) {
            my @line = split /\s+/;

	    # third column with a * indicates an unaligned read in the SAM format per samtools man page
            if ( $line[2] eq "*" ) {
                $unaligned_output_fh->print($_);
            }
            else {#write out the aligned map, excluding the default header- all lines starting with @.
                my $first_char = substr($line[0],0,1);
	        if ($first_char ne "@") {
                    $sam_map_output_fh->print($_);
                }
            }
        }

        $unaligned_output_fh->close();
        $sam_map_output_fh->close();
        

	unless ($self->_check_output_files_and_log_files_for_success("bwa_aln_core|bwa_sai2sam_pe_core", 'print alignments', 0, $self->unaligned_reads_file, $sam_map_output_fh->filename))
        {
	    $self->error_message("Validation of SAM format output failed");
            return;
        }
        

        #### STEP 3: Convert the SAM format output from samse/sampe (a text file) into the binary
        #### BAM file format.  This requires "samtools import"

        #picard merge command
        #java -Xmx8g -XX:MaxPermSize=256m -cp /gsc/scripts/lib/java/samtools/picard-tools-1.07/SamFormatConverter.jar net.sf.picard.sam.SamFormatConverter I=Q6OBHLc6AP.sam O=out2.bam VALIDATION_STRINGENCY=SILENT

        if ($self->sam_only) {

            my $from_file = $sam_map_output_fh->filename;

            if (defined($self->read_group_tag) ) {
                $self->status_message(">>>Adding readgroup tag.");
                my $tagged_sam_file = File::Temp->new( DIR => $tmp_dir, SUFFIX => ".sam" );
                my $cmd_tag = Genome::Model::Tools::Sam::AddReadGroupTag->create(input_file => $sam_map_output_fh->filename, output_file=>$tagged_sam_file->filename, read_group_tag=>$self->read_group_tag);
                unless ($cmd_tag->execute) {
                    $self->error_message("Error adding read group tag.");
                } 
                $self->status_message("<<<Done adding readgroup tag.");
                $from_file = $tagged_sam_file;
            }
 
            $self->status_message("Leaving alignment file in sam format. Moving from temp file to final file.");
            my $cmd_move = "mv ".$from_file." ".$self->alignment_file;
            my $rv_move = Genome::Sys->shellcmd( cmd => $cmd_move ); 
            unless ($rv_move) {
            $self->status_error("Failed to execute $cmd_move");
            return;  
            } 
        } else {

            $self->status_message("Converting from Sam to Bam.");
            my $samtools_import_command_line;
            my @conversion_input_files;

            if ($self->picard_conversion eq 1) {
                $samtools_import_command_line = sprintf(
                    "java -Xmx8g -XX:MaxPermSize=256m -cp /gsc/scripts/lib/java/samtools/picard-tools-1.07/SamFormatConverter.jar net.sf.picard.sam.SamFormatConverter I=%s O=%s VALIDATION_STRINGENCY=SILENT 2>>%s",
                     $sam_map_output_fh->filename, $self->alignment_file, $self->aligner_output_file
                );
                push(@conversion_input_files,$sam_map_output_fh->filename);
            } else  {
                $samtools_import_command_line = sprintf(
                    "samtools import %s %s %s 2>>%s",
                    $self->_ref_seq_index_file, $sam_map_output_fh->filename,
                   $self->alignment_file,      $self->aligner_output_file
                );
                push(@conversion_input_files,$sam_map_output_fh->filename);
                push(@conversion_input_files,$self->_ref_seq_index_file);
            }

            $self->status_message( "Merging with cmd: " . $samtools_import_command_line );
            
            Genome::Sys->shellcmd(
                cmd         => $samtools_import_command_line,
                #input_files => [ $self->_ref_seq_index_file, $sam_map_output_fh->filename ],
                input_files => \@conversion_input_files,
                output_files              => [ $self->alignment_file ],
                skip_if_output_is_present => 1,
            );

            unless (-e $self->alignment_file && -s $self->alignment_file) {
                    $self->error_message("Alignment output " . $self->alignment_file . " not found or zero length.  Something went wrong");
                    return;
            }

        }

    }

    $self->status_message("Align Reads completed.");
    $self->status_message("Unaligned reads: ".$self->unaligned_reads_file );
    $self->status_message("Aligned reads: ".$self->alignment_file );
    
    return 1;

}

1;
