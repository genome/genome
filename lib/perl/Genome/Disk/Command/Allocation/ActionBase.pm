package Genome::Disk::Command::Allocation::ActionBase;

use Genome;
use JSON 'to_json';

use strict;
use warnings;

class Genome::Disk::Command::Allocation::ActionBase {
    is => 'Command::V2',
    is_abstract => 1,
    has => {
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'The analysis project for which the data is being acted upon.',
        },
        requestor => {
            is => 'Genome::Sys::User',
            is_optional => 1,
            doc => 'The user who is requesting the action.',
        },
    },
};

sub disk_allocation_action {
    die('Failed to override disk_allocation_action!');
}

sub execute {
    my $self = shift;

    unless ($self->requestor) {
        $self->requestor(Genome::Sys->current_user);
    }
    $self->status_message("%s allocation(s) for Analysis Project %s at the request of %s.",
                          $self->disk_allocation_action, $self->analysis_project->__display_name__, $self->requestor->name);

    return $self->_execute();
}

sub reason {
    my $self = shift;

    my %reason = (
        'analysis_project'=>$self->analysis_project->id,
        'requestor'=> $self->requestor->id
    );
    return to_json(\%reason);
}

1;
