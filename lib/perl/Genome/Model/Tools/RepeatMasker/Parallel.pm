package Genome::Model::Tools::RepeatMasker::Parallel;

use strict;
use warnings;

use Genome;
use Workflow;

class Genome::Model::Tools::RepeatMasker::Parallel {
    is => ['Workflow::Operation::Command'],
    workflow => sub {
        my $run = Workflow::Operation->create(
            name => 'run',
            operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::RepeatMasker::Run')
        );
        $run->parallel_by('fasta_file');
        return $run;
    },
    has => [
            jobs => {
                is => 'Number',
                doc => 'The number of jobs to run in parallel',
                default_value => 100,
            },
            _fasta_file => { is_optional => 1, },
        ],
};

sub pre_execute {
    my $self = shift;

    $self->_fasta_file($self->fasta_file);
     my $split = Genome::Model::Tools::Fasta::Split->create(
        number_of_files => $self->jobs,
        fasta_file => $self->fasta_file,
        output_directory => $self->output_directory,
    );
    unless ($split) {
        die('Failed to create fasta split command');
    }
    unless ($split->execute) {
        die('Failed to execute fasta split command');
    }
    my $fasta_files = $split->_split_fasta_files;
    $self->fasta_file($fasta_files);
    return 1;
}

sub post_execute {
    my $self = shift;

    print Data::Dumper->new([$self])->Dump;
    my $output_basename;
    my @fasta_files = @{$self->fasta_file};
    my %files;
    for my $fasta_file (@fasta_files) {
        my $fasta_basename = File::Basename::basename($fasta_file);
        $fasta_basename =~ /(.*)_(\d+)$/;
        my $file_counter = $2;
        $output_basename = $1;
        my $file = $self->output_directory .'/'. $fasta_basename .'.out';
        unless (-e $file) {
            `touch $file`;
        }
        $files{$file_counter} =  $file;
    }
    my @output_files;
    for my $key ( sort {$a <=> $b} keys %files) {
        push @output_files, $files{$key};
    }

    my $output_file = $self->output_directory .'/'. $output_basename .'.out';
    my $merge = Genome::Model::Tools::RepeatMasker::MergeOutput->create(
        input_files => \@output_files,
        output_file => $output_file,
    );
    unless ($merge) {
        die ('Failed to create merge tool for RepeatMasker output');
    }
    unless ($merge->execute) {
        die ('Failed to execute the RepeatMasker output merging tool');
    }

    my @masked_files;
    for my $key ( sort {$a <=> $b} keys %files) {
        my $file = $files{$key};
        $file =~ s/\.out/\.masked/;
        push @masked_files, $file;
    }
    my $masked_file = $self->output_directory .'/'. $output_basename .'.masked';
    @masked_files = grep { -f } @masked_files;
    Genome::Sys->cat(
        input_files => \@masked_files,
        output_file => $masked_file,
    );
    for my $fasta_file (@fasta_files) {
        my @rm_files = glob($fasta_file.'*');
        for my $rm_file (@rm_files) {
            unless (unlink($rm_file)) {
                die('Failed to remove file '. $rm_file .":  $!");
            }
        }
    }
    my $table = Genome::Model::Tools::RepeatMasker::GenerateTable->create(
        fasta_file => $self->_fasta_file,
        output_file => $output_file,
        output_table => $self->output_directory .'/'. $output_basename .'.tsv',
    );
    unless ($table) {
        die('Failed to create table generating tool');
    }
    unless ($table->execute) {
        die('Failed to execute table generating tool');
    }
    return 1;
}


1;
