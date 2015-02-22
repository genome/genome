package Genome::SoftwareResult::StageableSimple;

use strict;
use warnings FATAL => 'all';
use Genome;

class Genome::SoftwareResult::StageableSimple {
    is_abstract => 1,
    is => 'Genome::SoftwareResult::Stageable',
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;

    $self->_prepare_staging_directory;
    $self->_run;

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub _run {
    my $self = shift;
    die "You must define _run in class" . $self->class;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return File::Spec->join('model_data', 'software-result', $self->id);
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}


1;
