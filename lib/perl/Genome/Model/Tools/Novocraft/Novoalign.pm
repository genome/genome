package Genome::Model::Tools::Novocraft::Novoalign;

use strict;
use warnings;

use Genome;

my $DEFAULT_OUTPUT_FORMAT = 'SAM';
my $DEFAULT_INPUT_FORMAT = 'ILMFQ';
my $DEFAULT_THREADS = 4;

class Genome::Model::Tools::Novocraft::Novoalign {
    is => 'Genome::Model::Tools::Novocraft',
    has_input => [
        novoindex_file => {
            doc => 'Required input file name containing the reference sequence file indexed by the Novocraft Novoindex software.',
            is => 'String',
        },
        fastq_files => {
            doc => 'A file(fragment reads) or two files(paired-end reads) separated by a space contained within single quotes in fastq format. Default input format is '. $DEFAULT_INPUT_FORMAT,
        },
        output_directory => {
            is => 'String',
            doc => 'The alignment directory where output files are written.  When run in parallel, a temp directory is created to run parallel processes within this defined output directory.  Once all parallel processes are complete, the final output files are written to the defined output directory and all temp  directories containing intermediate files are removed.',
        },
        params => {
            doc => 'The novocraft novoalign parameters. These should be specified in a quoted string with single dashes, e.g. "-x 1 -y -z"',
            is => 'String',
            is_optional => 1,
            default_value => '',
        },
        input_format => {
            is_optional => 1,
            default_value => $DEFAULT_INPUT_FORMAT,
            doc =>  'The input fastq file format. default_value='. $DEFAULT_INPUT_FORMAT,
        },
        output_format => {
            is_optional => 1,
            valid_values => ['Native', 'Pairwise','SAM'],
            default_value => $DEFAULT_OUTPUT_FORMAT,
            doc => 'The output format for alignment files. default_value='. $DEFAULT_OUTPUT_FORMAT,
        },
        threads => {
            is_optional => 1,
            is => 'Integer',
            default_value => $DEFAULT_THREADS,
            doc => 'The number of threads or CPUs to use. If run in parallel, this value must be accompanied by the correct lsf_resource as well.  default_value='. $DEFAULT_THREADS,
        }
    ],
    has_output => [
        log_file => {
            doc => 'Output log file containing results and error messages from the Novocraft Novoalign software.',
            is_optional => 1,
            is => 'String',
        },
        output_file => {
            doc => 'Output file containing the aligned reads in the defined output format. The default output format is '. $DEFAULT_OUTPUT_FORMAT,
            is_optional => 1,
            is => 'String',
        },
        full_param_string => {
            doc => 'The full string of parameters passed to Novoalign upon execution',
            is_optional => 1,
            is => 'String',
        },
    ],
    has_optional => [
        _fastq_count => { },
    ],
};

sub default_input_format {
    my $class = shift;
    return $DEFAULT_INPUT_FORMAT;
}

sub default_output_format {
    my $class = shift;
    return $DEFAULT_OUTPUT_FORMAT;
}

sub default_threads {
    my $class = shift;
    return $DEFAULT_THREADS;
}

sub execute {
    my $self = shift;

    if (ref($self->fastq_files) ne 'ARRAY') {
        my @fastq_files = split(/\s+/,$self->fastq_files);
        $self->fastq_files(\@fastq_files);
    }
    my @fastq_files = @{$self->fastq_files};

    my $fastq_count = scalar(@fastq_files);

    if ($fastq_count > 2) {
        die('Too many fastq files passed to novoalign command');
    }
    $self->_fastq_count($fastq_count);
    my @suffix = qw/\.txt \.fastq/;
    my @output_basenames;
    for my $fastq_file (@fastq_files) {
        my ($basename,$dirname,$suffix) = File::Basename::fileparse($fastq_file,@suffix);
        unless ($basename  && $dirname && $suffix) {
            die('Failed to parse fastq file name '. $fastq_file);
        }
        push @output_basenames, $basename;
    }

    my $output_basename;
    if ($fastq_count == 2) {
        my $lane;
        for my $end_output_basename (@output_basenames) {
            unless ($end_output_basename =~ m/((\d)_[12])/) {
                die('Failed to parse lane and end from file basename '. $end_output_basename);
            }
            my $read_id = $1;
            my $read_lane = $2;
            if ($lane) {
                unless ($lane == $read_lane) {
                    die('Fastq files do not contain reads from the same lane');
                }
            } else {
                $lane = $read_lane;
            }
            $end_output_basename =~ s/$read_id/$lane/;
            $output_basename = $end_output_basename;
        }
    } elsif ($fastq_count == 1) {
        $output_basename = $output_basenames[0];
    } else {
        die('Invalid number of fastq files '. $fastq_count);
    }
    unless (-d $self->output_directory) {
        Genome::Sys->create_directory($self->output_directory);
    }
    $self->output_file($self->output_directory .'/'. $output_basename .'.'. lc($self->output_format));
    $self->log_file($self->output_directory .'/'. $output_basename .'.aligner_out');

    $self->full_param_string(' -F '. $self->input_format .' -o '. $self->output_format .' -c '. $self->threads .' '. $self->params);
    my $cmdline = $self->novoalign_path
        . sprintf(' %s -d %s -f %s 1> %s 2>> %s',
                  $self->full_param_string,
                  $self->novoindex_file,
                  join(' ', @fastq_files),
                  $self->output_file,
                  $self->log_file,
              );
    # db disconnect prior to alignment
    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->status_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }
    Genome::Sys->shellcmd(
        cmd                         => $cmdline,
        input_files                 => [$self->novoindex_file, @fastq_files],
        output_files                => [$self->output_file, $self->log_file],
    );
    $self->_verify_output_files;
    return 1;
}

sub _verify_output_files {
    my $self = shift;

    my $log_fh = Genome::Sys->open_file_for_reading($self->log_file);
    unless ($log_fh) {
        $self->error_message('Failed to open log file '. $self->log_file .": $!");
        die($self->error_message);
    }
    my @lines = <$log_fh>;
    $log_fh->close;
    my $complete = 0;
    my $is_PE = 0;
    my $reads = 0;
    for my $line (@lines) {
        chomp($line);
        if ($line =~ m/^# Done\./) {
            $complete = 1;
        }
        if ($line =~ m/^#\s+Read Sequences:\s+(\d+)/) {
            $reads = $1;
        }
        #Not sure what novoalign does when no alignments found....
        #if ($line =~ m/\[match_index_sorted\] no reasonable reads are available. Exit!/) {
        #    $complete = 1;
        #}
        if ($line =~ m/^#\s+Paired Reads:\s+(\d+)/) {
            $is_PE = 1;
        }
    }
    if ( ($self->_fastq_count == 2) ) {
        unless ($is_PE == 1) {
            die ('Failed to align '. $self->_fastq_count .' fastq files as paired end');
        }
    } else {
        unless ($is_PE == 0) {
            die ('Failed to align '. $self->_fastq_count .' fastq file as fragment');
        }
    }
    unless ( $complete ) {
        die('Incomplete log file '. $self->log_file);
    }
    my $output_fh = Genome::Sys->open_file_for_reading($self->output_file);
    unless ($output_fh) {
        $self->error_message('Failed to open output file '. $self->output_file .": $!");
        die($self->error_message);
    }
    my $alignments;
    if ($self->output_format eq 'SAM') {
        while (my $line = $output_fh->getline) {
            chomp($line);
            if ($line =~ m/^@/) { next; }
            $alignments++;
        }
    } else {
        $output_fh->close;
        die('No logic implemented to validate '. $self->output_format .' format output file '. $self->output_file);
    }
    $output_fh->close;
    unless ($alignments == $reads) {
        die('Expected '. $reads .' alignments but found '. $alignments .'!');
    }
    return 1;
}

1;
