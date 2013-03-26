package Genome::Sample::Command::Import;

use strict;
use warnings;
use Genome;

class Genome::Sample::Command::Import {
    is => 'Command::Tree',
    doc => 'commands for importing samples',
};

my (%import_namespaces, %import_class_names);
_create_import_commands();
sub _create_import_commands {
    my @configs = _load_import_configs();
    for my $config ( @configs ) { 
        _create_import_command_for_config($config); 
    }
    return 1;
}

sub _load_import_configs {
    #TODO use the upcoming config API!
    return (
        {
            namespace => 'DbGap',
            nomenclature => 'dbGaP',
            sample_name_match => '\d+',
            sample_attribute_names => [qw/ tissue /],
            individual_name_match => '\d+',
            individual_attribute_names => [qw/ race gender /],
        }
    );
}

sub _create_import_command_for_config {
    my $config = shift;

    my $nomenclature = $config->{nomenclature};
    die 'No nomenclautre!' if not $nomenclature;
    my $namespace = $config->{namespace} // ucfirst(lc($nomenclature));
    my $class_name = 'Genome::Sample::Command::Import::'.$namespace;
    return $import_class_names{$class_name} if exists $import_class_names{$class_name};
    $import_namespaces{$namespace} = 1;

    die 'No individual name match!' if not $config->{individual_name_match};
    die 'No sample name match!' if not $config->{sample_name_match};

    my %properties;
    for my $type (qw/ sample individual /) {
        my $key_name = $type.'_attribute_names';
        next if not $config->{$key_name};
        my %type_properties = _get_properties_for_import_command_from_entity($type, @{$config->{$key_name}});
        %properties = ( %properties, %type_properties );
        $properties{ '_'.$key_name } = { is => 'ARRAY', is_constant => 1, value => $config->{$key_name}, };
    }

    #my %sample_properties = %{$config->{sample_properties}} if $config->{sample_properties};
    #my %individual_properties = %{$config->{individual_properties}} if $config->{individual_properties};
    $import_class_names{$class_name} = UR::Object::Type->define(
        class_name => $class_name,
        is => 'Genome::Sample::Command::Import::Base',
        has => [
            %properties,
            _individual_name_match => { is_constant => 1, value => $config->{individual_name_match} },
            _sample_name_match => { is_constant => 1, value => $config->{sample_name_match} },
        ],
        doc => "import $nomenclature samples",
    );

    my $name_property = $import_class_names{$class_name}->property_meta_for_name('name');
    #$name_property->doc();

    my $nomenclature_property = $import_class_names{$class_name}->property_meta_for_name('nomenclature');
    $nomenclature_property->default_value($nomenclature);

    if ( $config->{minimum_unique_source_name_parts} ) {
        die 'Invalid minimum_unique_source_name_parts! '.$config->{minimum_unique_source_name_parts} if $config->{minimum_unique_source_name_parts} > 1;
        my $minimum_unique_source_name_parts = $import_class_names{$class_name}->property_meta_for_name('_minimum_unique_source_name_parts');
        $minimum_unique_source_name_parts->default_value( $config->{minimum_unique_source_name_parts} // 2 );
    }

    return $import_class_names{$class_name};
}

sub _get_properties_for_import_command_from_entity {
    my ($entity, @names) = @_;

    my $class_name = 'Genome::'.ucfirst($entity);
    my $meta = $class_name->__meta__;
    my %properties;
    for my $name ( @names ) {
        my $property = $meta->property_meta_for_name($name);
        $properties{$name} = {
            is => 'Text',
            doc => "The value of '".join(' ', split('_', $name))."' for the $entity.",
        };
        if ( $property ) {
            $properties{$name}->{is} = $property->{is} if $property->{is};
            $properties{$name}->{doc} = $property->doc if $property->doc;
            $properties{$name}->{valid_values} = $property->valid_values if $property->valid_values;
        }
    }

    return %properties;
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

