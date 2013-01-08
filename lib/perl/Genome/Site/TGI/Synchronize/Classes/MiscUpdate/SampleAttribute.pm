package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::SampleAttribute;

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::SampleAttribute {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate',
};

sub perform_update {
    my $self = shift;
    $self->error_message('Cannot UPDATE sample attribute! It must be deleted and then inserted!');
    return $self->failure;
}

1;

