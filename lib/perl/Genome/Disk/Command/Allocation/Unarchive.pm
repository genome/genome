package Genome::Disk::Command::Allocation::Unarchive;

use strict;
use warnings;
use Genome;

class Genome::Disk::Command::Allocation::Unarchive {
    is => 'Genome::Disk::Command::Allocation::UnarchiveBase',
    has => [
        allocations => {
            is => 'Genome::Disk::Allocation',
            is_many => 1,
            shell_args_position => 1,
            doc => 'allocations to be unarchived',
        },
    ],
    doc => 'unarchives the given allocations',
};

sub help_detail {
    return 'unarchives the given allocations, moving them from tape back onto the filesystem';
}

sub help_brief {
    return 'unarchives the given allocations';
}

sub _execute {
    my $self = shift;

    $self->status_message("Starting unarchive command...");

    for my $allocation ($self->allocations) {
        $self->status_message("Unarchiving allocation " . $allocation->id);
        $self->_link_allocation_to_analysis_project($allocation);
        my $rv = $allocation->unarchive(reason => $self->reason);
        unless ($rv) {
            Carp::confess "Could not unarchive alloation " . $allocation->id;
        }
        $self->status_message("Finished unarchiving allocation " . $allocation->id);
    }

    $self->status_message("Done unarchiving, exiting...");
    return 1;
}

sub _link_allocation_to_analysis_project {
    my $self = shift;
    my $allocation = shift;

    my $owner = $allocation->owner;
    unless ($owner) {
        $self->fatal_message('This allocation appears to be orphaned: %s', $allocation->id);
    }

    if($owner->isa('Genome::SoftwareResult')) {
        $owner->add_user(label => 'sponsor', user => $self->analysis_project);
    } elsif ($owner->isa('Genome::Model::Build')) {
        unless($owner->model->analysis_project) {
            $self->fatal_message('No analysis project set on model for build %s.  Please use `genome analysis-project add-model` to correct this.', $owner->__display_name__);
        }
    } else {
        $self->fatal_message('Setting the analysis project of %s is currently not handled.  Please open a support request to unarchive this allocation.', $owner->__display_name__);
    }

    return 1;
}

1;
