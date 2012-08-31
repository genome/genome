package Genome::Model::Build::SimpleAlignment;

use strict;
use warnings;

use Genome;

class Genome::Model::Build::SimpleAlignment {
    is => 'Genome::Model::Build',
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
#    $DB::single=1;
    unless ($self) {
        return;
    }
    my $model = $self->model;
    unless ($model) {
        $self->error_message("Failed to get a model for this build!");
        return;
    }

    return $self;
}

sub workflow_log_directory_path {
    my $self = shift;
    return $self->data_directory."/workflow_logs";
}

sub merged_alignment_file {
    my $self = shift;
    return $self->data_directory."/merged.bam";
}


sub calculate_estimated_kb_usage {
    my $self = shift;

    # 1.5 gig... overestimating by 50% or so...
    #return 1536000;
    return 15360;
}

1;
