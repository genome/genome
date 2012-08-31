package Genome::Model::Tools::Fastx::ParallelFastqClipper;

use strict;
use warnings;

use Genome;
use File::Basename;
use File::Copy;
use File::Temp;
use Workflow;

class Genome::Model::Tools::Fastx::ParallelFastqClipper {
    is => ['Workflow::Operation::Command'],
    workflow => sub {
        my $workflow = Workflow::Operation->create(
                                                   name => 'parallel fastq clipper',
                                                   operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Fastx::Clipper'),
                                               );
        $workflow->parallel_by('input_file');
        return $workflow;
    },
    has => [
        reads => { is => 'Number', default_value => 5_000_000 },
        _final_output_file => { is => 'Text', is_optional => 1, },
        _basename => { is_optional => 1, },
        _suffix => { is_optional => 1, },
    ],
};

sub pre_execute {
    my $self = shift;

    my $input_file = $self->input_file;
    my ($basename, $dirname, $suffix) = File::Basename::fileparse($input_file,qw/\.txt \.fastq \.fq /);
    $self->_basename($basename);
    $self->_suffix($suffix);
    unless ($self->output_directory) {
        $self->output_directory($dirname);
    }
    my $lines = $self->reads * 4;

    require Cwd;
    my $cwd = Cwd::cwd();
    my $tempdir = File::Temp::tempdir( DIR => $self->output_directory, CLEANUP => 1 );
    chdir $tempdir;
    Genome::Sys->shellcmd(
        cmd => "split -l $lines $input_file $basename",
        input_files => [$input_file],
    );
    chdir $cwd;
    my @split_files = glob($tempdir .'/'. $basename .'*');
    my @fastq_files;
    for my $split_file (@split_files) {
        my $fastq_file = $split_file . $self->_suffix;
        unless (File::Copy::move($split_file,$fastq_file)) {
            die('Failed to move split file '. $split_file .' to fastq file '. $fastq_file);
        }
        push @fastq_files,$fastq_file;
    }
    $self->input_file(\@fastq_files);
    $self->_final_output_file($self->output_file);
    return 1;
}

sub post_execute {
    my $self = shift;

    print Data::Dumper->new([$self])->Dump;

    my @failures = grep { $_ ne 1 } @{$self->result};
    if (@failures) {
        $self->error_message('One or more of the parallel commands failed');
        die($self->error_message);
    }

    my $output_file = $self->output_directory .'/'. $self->_basename .'_clipped'. $self->_suffix;
    Genome::Sys->cat(
        input_files => $self->output_file,
        output_file => $output_file,
    );
    my $log_file = $self->output_directory .'/'. $self->_basename .'_clipped.log';
    Genome::Sys->cat(
        input_files => $self->log_file,
        output_file => $log_file,
    );
    for my $file (@{$self->output_file},@{$self->log_file},@{$self->input_file}) {
        unless (unlink $file) {
            $self->error_message('Failed to remove file '. $file .":  $!");
            die($self->error_message);
        }
    }
    return 1;
}


1;
