package Genome::Model::Tools::Shrimp::AlignReads;

use strict;
use warnings;

use Genome;
use Workflow;

class Genome::Model::Tools::Shrimp::AlignReads {
    is  => ['Workflow::Operation::Command'],
    workflow => sub {
        my $rmapper = Workflow::Operation->create(
            name => 'rmapper',
            operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Shrimp::Rmapper')
        );
        $rmapper->parallel_by('fasta_file');
        return $rmapper;
    },
    has_optional => [
        jobs => {
            is => 'Number',
            doc => 'The kb of sequence to include in each instance. default_value=100',
            default_value => 100,
        },
        _fasta_file => { },
    ],
};

sub help_synopsis {
    return <<EOS
    A SHRiMP based utility for aligning reads.;
EOS
}

sub help_brief {
    return <<EOS
    A SHRiMP based utility for aligning reads.;
EOS
}

sub help_detail {
    return <<EOS
Provides an interface to the SHRiMP aligner.  Inputs are:
EOS
}

sub pre_execute {
    my $self = shift;

    $self->_fasta_file($self->fasta_file);

    my $split = Genome::Model::Tools::Fasta::Split->create(
        number_of_files => $self->jobs,
        fasta_file => $self->fasta_file,
    );

    unless ($split) {
        die('Failed to create fasta split command');
    }
    unless ($split->execute) {
        die('Failed to execute fasta split command');
    }
    my $fasta_files = $split->fasta_files;
    $self->fasta_file($fasta_files);
    return 1;
}


sub post_execute {
    my $self = shift;

    print Data::Dumper->new([$self])->Dump;
    my $output_basename;
    my @fasta_files = @{$self->fasta_file};
    my %alignment_files;
    my %aligner_output_files;
    for my $fasta_file (@fasta_files) {
        my $fasta_basename = File::Basename::basename($fasta_file);
        $fasta_basename =~ /(.*)_(\d+)$/;
        my $file_counter = $2;
        $output_basename = $1;
        my $file = $self->output_directory .'/'. $fasta_basename .'.shrimp';
        $alignment_files{$file_counter} =  $file;
        $aligner_output_files{$file_counter} = $self->output_directory .'/'. $fasta_basename .'.aligner_output';
    }
    my @alignment_files;
    my @aligner_output_files;
    for my $key ( sort {$a <=> $b} keys %alignment_files) {
        push @alignment_files, $alignment_files{$key};
        push @aligner_output_files, $aligner_output_files{$key};
    }

    my $output_file = $self->output_directory .'/'. $output_basename .'.shrimp';
    Genome::Sys->cat(
        input_files => \@alignment_files,
        output_file => $output_file,
    );
    Genome::Sys->cat(
        input_files => \@aligner_output_files,
        output_file => $self->output_directory .'/'. $output_basename .'.aligner_output',
    );

    return 1;
}



1;
