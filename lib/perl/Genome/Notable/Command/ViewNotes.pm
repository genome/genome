package Genome::Notable::Command::ViewNotes;

use strict;
use warnings;

use Genome;

class Genome::Notable::Command::ViewNotes {
    is => 'Command::V2',
    is_abstract => 1,
    has => [
        notables => {
            is => 'Genome::Notable', #this class won't work with the command-line object resolution, but this command is usable in code
            is_many => 1,
            shell_args_position => 1,
            doc => 'notable objects on which to view the notes',
        },
        note_type => {
            is => 'Text',
            is_optional => 1,
            doc => 'view notes with this type (header value)',
        },
        _notes => {
            is => 'Genome::MiscNote',
            is_many => 1,
            is_optional => 1,
            is_transient => 1,
            doc => 'Notes belonging to provided notable objects that were not filtered out by note_type',
        },
    ],
    doc => 'View notes that have been set on notable objects',
};

sub help_detail {
    return <<EOS
This command can be used to view the notes that have been added to a notable object.
EOS
}

sub execute {
    my $self = shift;

    my @notes_found;
    my %note_params;
    $note_params{header_text} = $self->note_type if $self->note_type;

    my @notables = $self->notables;
    for my $notable (@notables) {
        print "\n" . 'Notes for ' . $notable->__display_name__ . "\n";
        my @notes = $notable->notes(%note_params);
        push @notes_found, @notes if @notes;

        for my $note (@notes) {
            print $note->header_text . ' by ' . $note->editor_id . ' on ' . $note->entry_date;
            my $body_text = $note->body_text;
            if ($body_text) {
                print ":\n";
                my @body_lines = split("\n", $body_text);
                print '> ' . join("\n> ", @body_lines) . "\n";
            }
            else {
                print ".\n";
            }
        }
        print "\n";
    }

    $self->_notes(\@notes_found);
    return 1;
}

1;
