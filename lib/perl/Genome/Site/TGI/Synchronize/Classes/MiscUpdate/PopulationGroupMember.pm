package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::PopulationGroupMember;

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::PopulationGroupMember {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate',
};

sub perform_update {
    my $self = shift;
    $self->error_message('Cannot UPDATE population group member attribute! It must be deleted and then inserted!');
    return $self->failure;
}

1;

