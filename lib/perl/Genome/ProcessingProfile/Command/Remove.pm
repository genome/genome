package Genome::ProcessingProfile::Command::Remove;

use strict;
use warnings;

use Genome;

use Data::Dumper;

class Genome::ProcessingProfile::Command::Remove {
    is => 'Command::V2',
    has => [
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            shell_args_position => 1,
            doc => 'Processing profile to delete.',
        },
    ],
    doc => 'Delete a processing profile',
};

sub help_brief { 'Delete a processing profile'; }
sub help_detail { help_brief(); }

sub execute {
    my $self = shift;

    my $processing_profile = $self->processing_profile;
    return if not $processing_profile;

    my $display_name = $processing_profile->__display_name__;

    unless ( $processing_profile->delete ) {
        $self->error_message('Could not remove processing profile '.$display_name);
        return;
    }

    $self->status_message('Removed processing profile '.$display_name);

    return 1;
}

1;

