package Genome::Site::TGI::Synchronize::Classes::MiscUpdate::Unsupported;

use strict;
use warnings;

class Genome::Site::TGI::Synchronize::Classes::MiscUpdate::Unsupported {
    is => 'Genome::Site::TGI::Synchronize::Classes::MiscUpdate',
};

sub perform_update {
    my $self = shift;
    $self->error_message('Unsupported subject class name! '.$self->subject_class_name);
    return $self->failure;
}

1;

