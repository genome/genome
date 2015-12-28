package Genome::Role::Notable::Command::AddNote;

use strict;
use warnings;
use Genome;
use UR::Role;

our $notable_type : RoleParam(notable_type);
role Genome::Role::Notable::Command::AddNote {
    has => [
        notable => {
            is => $notable_type,
            shell_args_position => 1,
            doc => 'notable object to which a note will be added',
        },
        header_text => {
            is => 'Text',
            doc => 'Header of note to be created',
        }
    ],
    has_optional => [
        body_text => {
            is => 'Text',
            doc => 'Body text of note to be created',
        }
    ],
};

sub help_detail {
    return <<EOS
This command can be used to add notes to notable objects
EOS
}

sub execute {
    my $self = shift;

    my $rv = $self->notable->add_note(
        header_text => $self->header_text,
        body_text => $self->body_text,
    );
    unless ($rv) {
        $self->error_message('Failed to add note to ' . $self->notable->__display_name__);
        Carp::confess $self->error_message;
    }

    $self->status_message('Added note with header "' . $self->header_text . '" to object "' . 
        $self->notable->__display_name__ . '" of type ' . $self->notable->class);
    return 1;
}

1;

