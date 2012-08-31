package Genome::Model::Tools::WuBlast::Blastx::Parallel;

use strict;
use warnings;

use Genome;
use Workflow;

class Genome::Model::Tools::WuBlast::Blastx::Parallel {
    is  => ['Workflow::Operation::Command'],
    workflow => sub {
        my $blastx = Workflow::Operation->create(
            name => 'parallel blastx',
            operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::WuBlast::Blastx')
        );
        $blastx->parallel_by('query_file');
        return $blastx;
    },
    has_optional => [
        jobs => {
            is => 'Number',
            doc => 'The number of jobs to run in parallel. default_value=100',
            default_value => 500,
        },
        _query_file => { },
        _output_file => { },
    ],
};

sub create {
    my $class = shift;
    my %params = @_;
    my $output_file = delete($params{output_file});
    my $self = $class->SUPER::create(@_);
    $self->_output_file($output_file);
    return $self;
}


sub pre_execute {
    my $self = shift;

    $self->_query_file($self->query_file);
    #Output directory should really be a temp dir inside of the real output_directory
    #Then we can clean up these intermediate fasta upon completion
    my $split = Genome::Model::Tools::Fasta::Split->create(
        number_of_files => $self->jobs,
        fasta_file => $self->query_file,
        output_directory => $self->output_directory,
    );

    unless ($split) {
        die('Failed to create fasta split command');
    }
    unless ($split->execute) {
        die('Failed to execute fasta split command');
    }
    my $fasta_files = $split->_split_fasta_files;
    $self->query_file($fasta_files);
    return 1;
}


sub post_execute {
    my $self = shift;

    print Data::Dumper->new([$self])->Dump;
    my $results = $self->{result};
    for (my $i = 0; $i < scalar(@$results); $i++) {
        my $rv = $results->[$i];
        if ($rv != 1) {
            foreach my $error (@Workflow::Simple::ERROR) {
                $self->error_message($error->error);
            }
            die($self->error_message);
        }
    }
    my $output_basename;
    my @fasta_files = @{$self->query_file};
    my %blast_files;
    my $error = 0;
    for my $fasta_file (@fasta_files) {
        my $fasta_basename = File::Basename::basename($fasta_file);
        unless ($fasta_basename =~ /(.*)_(\d+)$/) {
            $self->error_message('Failed to parse the fasta file name '. $fasta_file);
            $error++;
        }
        my $file_counter = $2;
        $output_basename = $1;
        my $file = $self->output_directory .'/'. $fasta_basename .'.blast';
        $blast_files{$file_counter} =  $file;
    }
    #Get blast files in order
    my @blast_files;
    for my $key ( sort {$a <=> $b} keys %blast_files) {
        push @blast_files, $blast_files{$key};
    }

    my $output_file = $self->_output_file || $self->output_directory .'/'. $output_basename .'.blast';
    Genome::Sys->cat(
        input_files => \@blast_files,
        output_file => $output_file,
    );
    for my $file (@fasta_files,@blast_files) {
        unless (unlink $file) {
            $self->error_message('Failed to unlink file '. $file);
            $error++;
        }
    }
    if ($error) {
        return;
    }
    return 1;
}



1;
