package Genome::Library::Command::Note::Add;

use strict;
use warnings;

use Genome;

class Genome::Library::Command::Note::Add {
    #is => 'UR::Object::Command::Update',
    is => 'Command',
    has_constant => [
        subject_class_name  => {
            value => 'Genome::Library' 
        },
        friendly_class_name => {
            value => 'library',
        },
    ],
    has => [
        library => {
            is => 'Text',
            doc => 'The id or name of the library to which to attach a note',
            shell_args_position => 1,
        },
        header => {
            is => 'Text',
            shell_args_position => 2,
            doc => 'The header for the note.'
        },
        detail => {
            is => 'Text',
            shell_args_position => 3,
            doc => 'The detail for the note.'
        },
        subject => {
            is => 'Genome::Library',
            id_by => 'subject_id',
            is_optional => 1,
            doc => 'The library to be updated.'
        },
    ],
    doc => "add a note",
};

sub sub_command_sort_position { 3 }

sub _resolve_and_set_subject {
    my $self = shift;
    my $spec = shift;

    my $friendly_class_name = $self->friendly_class_name;
    my $subject_class_name = $self->subject_class_name;

    my $subject;

    if ($spec !~ /\D/) {
        # all numbers, try to get by ID
        $subject = $subject_class_name->get($spec);
    }

    unless ($subject) {
        # see if the spec is the start of a name
        $subject = $subject_class_name->get('name like' => $spec);
    }

    unless ($subject) {
        $self->error_message("Failed to find a $friendly_class_name with the ID, name, or partial name '$spec'!");
        return;
    }

    return $subject;
}

sub execute {
    my $self = shift;

    my $spec = $self->library;
    my $subject = $self->subject ||  $self->_resolve_and_set_subject($spec);
   
    return unless $subject;

    print "add note to $subject h $self->{header} b $self->{detail}\n";

    my $n = Genome::MiscNote->create(
        subject_class_name => $subject->class,
        subject_id => $subject->id,
        header_text => $self->header,
        body_text => $self->detail,
        editor_id => 'me@somewhere',
        entry_date => UR::Time->now(),
    );

    my $n2 = $subject->add_note(
        header_text => $self->header,
        body_text   => $self->detail,
    );

    print Data::Dumper::Dumper($n, $n2);
    $n->delete;
    return 1;
    ##### 
    
    my $friendly_class_name = $self->friendly_class_name;
    my $subject_class_name = $self->subject_class_name;
    my $subject_property_name = $self->subject_property_name;
    
    my $name = $subject->name;
    my $id = $subject->id;

    my $new = $self->value;
    my $old = $subject->$subject_property_name;

    if (defined $old) {
        if ($old eq $new) {
            $self->error_message("The $subject_property_name is already set to '$new' on $friendly_class_name $name ($id)!");
            return;
        }
        $self->warning_message("Changing $subject_property_name from '$old' to '$new' on $friendly_class_name $name ($id)");
        
    }

    $subject->$subject_property_name($new);

    if (my @errors = $subject->__errors__) {
        for my $error (@errors) {
            $self->error_message($error);
        }
        $self->error_message("The requested change cannot be made because it makes $friendly_class_name $name ($id) invalid!");
        $subject->$subject_property_name($old);
        return;
    }

    return 1;
}

1;

