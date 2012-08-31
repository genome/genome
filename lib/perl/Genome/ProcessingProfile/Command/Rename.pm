package Genome::ProcessingProfile::Command::Rename;

use strict;
use warnings;

use Genome;

class Genome::ProcessingProfile::Command::Rename {
    is => 'Command::V2',
    has => [
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            shell_args_position => 1,
            doc => 'The processing profile to rename',
        },
        new_name => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'The new name for the processing profile.',
        },
    ],
};

sub help_brief {
    return 'rename a processing profile';
}

sub execute {
    my $self = shift;

    # Verify processing profile 
    if ( not $self->processing_profile ) {
        $self->error_message('No prcessing profile to rename');
        return;
    }

    # Verify new name
    unless ( $self->new_name =~ /\w+/ ) {
        $self->error_message("Letters are required to be included in the new name");
        return;
    }
    
    if ( $self->new_name eq $self->processing_profile->name ) {
        $self->error_message(
            sprintf('Processing profile (<ID> %s) already is named "%s"', $self->processing_profile->id, $self->new_name)
        );
        return;
    }
    
    my $existing_pp = $self->processing_profile->class->get(name => $self->new_name);
    if ( $existing_pp ) {
        $self->error_message('A '.$self->processing_profile->type_name.' processing profile already exists for name ('.$self->new_name.'): '.$existing_pp->__display_name__);
        return;
    }

    # Rename
    $self->processing_profile->name( $self->new_name );

    # Sanity chack
    unless ( $self->new_name eq $self->processing_profile->name ) {
        $self->error_message(
            sprintf(
                'Could not rename processing profile (<ID> %s) to "%s" for unkown reasons', 
                $self->processing_profile->id, 
                $self->new_name,
            )
        );
        return;
    }

    printf(
        'Renamed processing profile (<ID> %s) to "%s"', 
        $self->processing_profile->id, 
        $self->processing_profile->name, 
    );
    print "\n";
    
    return 1;
}

1;

