package Genome::Disk::Command::Allocation::UnarchiveBase;

use Genome;
use JSON 'to_json';

use strict;
use warnings;

class Genome::Disk::Command::Allocation::UnarchiveBase {
    is => 'Command::V2',
    is_abstract => 1,
    has => {
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'The analysis project for which the data is being unarchived',
        },
        requestor => {
            is => 'Genome::Sys::User',
            is_optional => 1,
            doc => 'The user who is requesting the unarchive.',
        },
    },
};

sub execute {
    my $self = shift;

    my $name = $self->requestor ? $self->requestor->name : Genome::Sys->username;

    $self->status_message("Unarchiving allocation(s) for Analysis Project %s at the request of %s.",
            $self->analysis_project->__display_name__, $name);

    return $self->_execute();
}

sub reason {
    my $self = shift;

    my $req_id = $self->requestor ? $self->requestor->id : Genome::Sys->username;

    my %reason = (
        'analysis_project'=>$self->analysis_project->id,
        'requestor'=> $req_id,
    );
    return to_json(\%reason);
}

1;
