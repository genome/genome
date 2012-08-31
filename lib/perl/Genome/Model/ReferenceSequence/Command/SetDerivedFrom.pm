package Genome::Model::ReferenceSequence::Command::SetDerivedFrom;

use strict;
use warnings;

use Genome;
use Carp qw/confess/;


class Genome::Model::ReferenceSequence::Command::SetDerivedFrom {
    is => 'Genome::Command::Base',
    has_input => [
        child_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'The "child" build on which to set the "derived_from" property',
        },
        parent_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'The "parent" build to which the "child" build\'s "derived_from" property should point',
            is_optional => 1,
        },
        clear => {
            is => 'Boolean',
            doc => 'Clear the value of derived_from. This option is mutually exclusive with --parent-build.',
            is_optional => 1,
        },
        force => {
            is => 'Boolean',
            default => 0,
            doc => 'Set this to disable prompting about overriding existing values of derived_from',
        }
    ]
};

sub help_brief {
    "Mark a reference sequence as 'derived from' another."
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 gmt model reference-sequence set-derived-from --child-build my-ref-seq --parent-build NCBI-human-build36
EOS
}

sub help_detail {
    return <<EOS
Marking a reference sequence as derived_from another makes them compatible
for purposes such as selecting annotation or dbsnp builds. This property can
also specified at reference sequence creation with the --derived-from option.
EOS
}

sub execute {
    my $self = shift;

    unless($self->parent_build xor $self->clear) {
        confess "Please specify exactly one of --parent-build or --clear";
    }

    my $cname = $self->child_build->__display_name__;
    if ($self->clear) {
        my @inputs = grep { $_->name eq 'derived_from' or $_->name eq 'coordinates_from' } $self->child_build->inputs;
        map { $_->delete } @inputs;
        $self->status_message("derived_from property cleared for build $cname.");
        return 1;
    }

    my $pname = $self->parent_build->__display_name__;
    if ($self->child_build->derived_from) {
        my $old_parent = $self->child_build->derived_from;
        my $old_name = $old_parent->__display_name__;

        $self->status_message("Build $cname is already derived from $old_name.");
        return 1 if $old_parent->id == $self->parent_build->id;
        
        unless ($self->force) {
            my $result = $self->_ask_user_question("Are you sure you want to override this and set derived_from to $pname?", 60);
            unless ($result eq "yes") {
                $self->error_message("Aborted by user.");
                return;
            }
        }
    }

    $self->child_build->derived_from($self->parent_build);
    $self->child_build->coordinates_from($self->child_build->derived_from_root);

    $self->status_message("Build $cname is now derived_from $pname.");

    return 1;
}

1;
