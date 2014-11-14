package Genome::VariantReporting::Command::ShowDocs;

use strict;
use warnings;
use Genome;

my %type_lookup = (
    filters => 'filters',
    interpreters => 'interpreters',
    experts => 'adaptors',
    reporters => 'reporters',
);

class Genome::VariantReporting::Command::ShowDocs {
    is => 'Command::V2',
    has_input => [
        component_type => {
            is => 'Text',
            shell_args_position => 1,
            valid_values => [qw(filters interpreters reporters experts)],
            doc => "To list all components of this type by name, only specify the type.",
        },
        component_name => {
            is => 'Text',
            shell_args_position => 2,
            is_optional => 1,
            doc => "To get details about a specific component, specify the name"},
    ],
};

sub execute {
    my $self = shift;
    my $factory = Genome::VariantReporting::Framework::Factory->create();

    if (defined $self->component_name) {
        my $component = $factory->get_dummy_object($type_lookup{$self->component_type}, $self->component_name);
        for my $section ($component->vr_doc_sections) {
            _print_section($section);
        }
        $component->delete;
    }
    else {
        _print_section (
            {
                header => "Available ".$self->component_type,
                items => [sort grep {$_ !~ /^__/} $factory->names($type_lookup{$self->component_type})],
            }
        );
    }
    return 1;
}

sub _print_section {
    my $section = shift;
    print sprintf (
        "%s\n%s\n\n",
        Term::ANSIColor::colored($section->{header}, 'underline'),
        join("\n", @{$section->{items}}),
    );
}

1;

