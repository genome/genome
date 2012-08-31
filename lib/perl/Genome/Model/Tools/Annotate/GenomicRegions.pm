package Genome::Model::Tools::Annotate::GenomicRegions;

use strict;
use warnings;

use Genome;
use Workflow;

class Genome::Model::Tools::Annotate::GenomicRegions {
    is => ['Workflow::Operation::Command'],
    workflow => sub {
        my $workflow = Workflow::Operation->create(
                                                   name => 'parallel genomic region annotation',
                                                   operation_type => Workflow::OperationType::Command->get('Genome::Model::Tools::Annotate::ChromosomeRegions'),
                                               );
        $workflow->parallel_by('chromosome');
        return $workflow;
    },
};

sub pre_execute {
    my $self = shift;
    my $bed = Genome::Model::Tools::RefCov::ROI::Bed->create(
        file => $self->bed_file
    );
    my @chromosomes = $bed->chromosomes;
    $self->chromosome(\@chromosomes);
}

sub post_execute {
    my $self = shift;

    my @failures = grep { $_ ne 1 } @{$self->result};
    if (@failures) {
        $self->error_message('One or more of the parallel commands failed');
        die($self->error_message);
    }
    my ($basename,$dirname,$suffix) = File::Basename::fileparse($self->bed_file,qw/.bed/);
    my $final_file = $self->output_directory .'/'. $basename .'_'. $self->anno_db .'_'. $self->version .'.bed';
    Genome::Sys->cat(
        input_files => $self->anno_file,
        output_file => $final_file,
        skip_if_output_is_present => 0,
    );
    for my $file (@{$self->anno_file}) {
        unless (unlink $file) {
            $self->error_message('Failed to remove file '. $file .":  $!");
            die($self->error_message);
        }
    }
    $self->anno_file($final_file);
    return 1;
}
1;
