package Genome::Model::Tools::Novocraft::ParallelNovoalign;

use strict;
use warnings;

use Genome;
use Workflow;

my $DEFAULT_SEQUENCES = 1_000_000;
my $DEFAULT_LSF_QUEUE = $ENV{GENOME_LSF_QUEUE_BUILD_WORKER};
my $DEFAULT_LSF_RESOURCE = "-M 8000000 -R 'select[type==LINUX64 && mem>8000] rusage[mem=8000] span[hosts=1]'";
## number of threads will be prepended at run time

class Genome::Model::Tools::Novocraft::ParallelNovoalign {
    is  => ['Workflow::Operation::Command'],
    workflow => sub {
        my $rmapper = Workflow::Operation->create(
            name => 'parallel novocraft novoalign',
            operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Novocraft::Novoalign')
        );
        $rmapper->parallel_by('fastq_files');
        return $rmapper;
    },
    has_optional => [
        lsf_queue => {
            doc => 'The lsf queue to use. default_value='. $DEFAULT_LSF_QUEUE,
            default_value => $DEFAULT_LSF_QUEUE,
        },
        lsf_resource => {
            doc => 'The lsf resource request necessary to execute.  default_value='. $DEFAULT_LSF_RESOURCE,
            default_value => $DEFAULT_LSF_RESOURCE,
        },
        sequences => {
            is => 'Number',
            doc => 'The number of sequences(paired or fragment) to include in each instance. default_value='. $DEFAULT_SEQUENCES,
            default_value => $DEFAULT_SEQUENCES,
        },
        _output_basename => { },
        _base_output_directory => { },
    ],
};

sub help_detail {
    'This documentation should be long explaining how this tool runs novoalign in parallel.  Please finish documentation';
}

sub pre_execute {
    my $self = shift;

    # When run through Workflow, the default value is not set for optional parameters.
    unless ($self->input_format) {
        $self->input_format(Genome::Model::Tools::Novocraft::Novoalign->default_input_format);
    }

    unless ($self->output_format) {
        $self->output_format(Genome::Model::Tools::Novocraft::Novoalign->default_output_format);
    }

    unless ($self->threads) {
        $self->threads(Genome::Model::Tools::Novocraft::Novoalign->default_threads);
    }
    $self->lsf_resource('-n ' . $self->threads . ' ' . $self->lsf_resource);

    $self->_operation->operation_type->lsf_resource($self->lsf_resource);
    $self->_operation->operation_type->lsf_queue($self->lsf_queue);

    if (ref($self->fastq_files) ne 'ARRAY') {
        my @fastq_files = split(/\s+/,$self->fastq_files);
        $self->fastq_files(\@fastq_files);
    }
    my @fastq_files = @{$self->fastq_files};
    my $fastq_count = scalar(@fastq_files);

    $self->_base_output_directory($self->output_directory);
    my $template = '/Genome-Model-Tools-Novocraft-ParallelNovoalign-'. $ENV{'USER'} .'-XXXXX';
    my $tempdir = File::Temp->tempdir($template, DIR=>$self->output_directory, CLEANUP=>1);
    $self->output_directory($tempdir);

    $self->_operation->log_dir($self->output_directory);
    my @output_basenames;
    my @sub_fastqs;
    for my $fastq_file (@fastq_files) {
        my $split = Genome::Model::Tools::Fastq::Split->create(
            split_size => $self->sequences,
            fastq_file => $fastq_file,
            output_directory => $tempdir,
        );
        unless ($split) {
            die('Failed to create fastq split command');
        }
        unless ($split->execute) {
            die('Failed to execute fastq split command');
        }
        my @suffix = qw/\.txt \.fastq/;
        my ($basename,$dirname,$suffix) = File::Basename::fileparse($fastq_file,@suffix);
        unless ($basename && $dirname && $suffix) {
            die('Failed to parse fastq file '. $fastq_file);
        }
        push @output_basenames, $basename;
        push @sub_fastqs, $split->fastq_files;
    }
    $self->_output_basename(join('_',@output_basenames));

    my @parallel_fastq_files;
    if ($fastq_count == 2) {
        my @fastq_1s = @{$sub_fastqs[0]};
        my @fastq_2s = @{$sub_fastqs[1]};
        unless (scalar(@fastq_1s) == scalar(@fastq_2s)) {
            die('The number of Paired End Read 1 fastq '. scalar(@fastq_1s)
                    .' does not match the number of Paired End Read 2 fastq '. scalar(@fastq_2s));
        }
        for (my $i = 0; $i < scalar(@fastq_1s); $i++) {
            my @fastqs = ($fastq_1s[$i],$fastq_2s[$i]);
            push @parallel_fastq_files, \@fastqs;
        }
        my $lane;
        for my $output_basename (@output_basenames) {
            unless ($output_basename =~ m/((\d)_[12])/) {
                die('Failed to parse lane and end from file basename '. $output_basename);
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
            $output_basename =~ s/$read_id/$lane/;
            $self->_output_basename($output_basename);
        }
    } elsif ($fastq_count == 1) {
        @parallel_fastq_files = @{$sub_fastqs[0]};
        $self->_output_basename($output_basenames[0]);
    } else {
        die('Invalid number of fastq files '. $fastq_count);
    }
    $self->fastq_files(\@parallel_fastq_files);
    return 1;
}

sub post_execute {
    my $self = shift;

    $self->status_message(Data::Dumper->new([$self])->Dump);

    #Merge OUTPUT
    my @output_files = @{$self->output_file};

    #The input fastq files better have been uniquely named... otherwise this could stomp on existing files
    $self->output_file($self->_base_output_directory .'/'. $self->_output_basename .'.'. lc($self->output_format));
    if ($self->output_format eq 'SAM') {
        my $output_file_fh = Genome::Sys->open_file_for_writing($self->output_file);
        unless ($output_file_fh) {
            die('Failed to open file for writing '. $self->output_file);
        }

        my $header = 0;
        for my $sub_output_file (@output_files) {
            my $sub_output_file_fh = Genome::Sys->open_file_for_reading($sub_output_file);
            unless ($sub_output_file_fh) {
                die('Failed to open sub output file for reading '. $sub_output_file_fh);
            }
            while (my $line = $sub_output_file_fh->getline) {
                if ($line =~ /^@/ && $header) { next; }
                print $output_file_fh $line;
            }
            $header = 1;
        }
        $output_file_fh->close;
    } else {
        die('Implement logic in post_execute method to merge output format '. $self->output_format);
    }

    #LOG
    my @log_files = @{$self->log_file};
    $self->log_file($self->_base_output_directory .'/'. $self->_output_basename .'.aligner_output');
    Genome::Sys->cat(
        input_files => \@log_files,
        output_file => $self->log_file,
    );

    return 1;
}

1;
