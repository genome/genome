package Genome::Sample::Command::Import;

use strict;
use warnings;
use Genome;

class Genome::Sample::Command::Import {
    is => 'Command::Tree',
    doc => 'commands for importing samples',
};

#TODO use the upcoming config API!
my @configs = (
    {
        namespace => 'Test',
        nomenclature => 'TeSt',
        sample_name_match => '\d+',
        individual_name_match => '\d+',
    }
);

for my $config ( @configs ) { create_command($config); }
my (%import_namespaces, %import_class_names);
sub create_command {
    my $config = shift;

    my $namespace = $config->{namespace};
    my $class_name = 'Genome::Sample::Command::Import::'.$namespace;
    return $import_class_names{$class_name} if exists $import_class_names{$class_name};

    die 'No individual name match!' if not $config->{individual_name_match};
    die 'No sample name match!' if not $config->{sample_name_match};

    my $nomenclature = $config->{nomenclature} // $namespace;
    $import_namespaces{$namespace} = 1;
    $import_class_names{$class_name} = UR::Object::Type->define(
        class_name => $class_name,
        is => 'Genome::Sample::Command::Import::Base',
        has => [
            name => {
                is => 'Text',
                doc => 'Sample name.',
            },
            nomenclature => { is_constant => 1, value => $nomenclature },
            _individual_name_match => { is_constant => 1, value => $config->{individual_name_match} },
            _sample_name_match => { is_constant => 1, value => $config->{sample_name_match} },
        ],
        doc => "import $nomenclature samples",
    );

    return $import_class_names{$class_name};
}

# Overload sub command classes to return these in memory ones, plus the existing ones
sub sub_command_names {
    my $class = shift;
    my %sub_command_names = map { $_ => 1 } ($class->SUPER::sub_command_names, keys %import_namespaces);
    return keys %sub_command_names;
}

sub sub_command_classes {
    my $class = shift;
    my %sub_command_classes = map { $_ => 1 } ($class->SUPER::sub_command_classes, keys %import_class_names);
    return keys %sub_command_classes;
}

1;

