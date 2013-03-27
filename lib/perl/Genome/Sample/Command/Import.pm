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
            nomenclature => 'ATCC',
            name_regexp => '(ATCC\-[\w\d]+\-[\w\d]+)(\-[\w\d]+)?',
            sample_attributes => {
                age => {},
                disease => { is_optional => 1, },
                organ_name => {},
            },
            individual_attributes => [qw/ ethnicity gender /],
            # gender => { valid_values => [qw/ male female /], }
        },
        {
            nomenclature => 'dbGaP',
            name_regexp => '(dbGaP\-\d+)\-\d+',
            sample_attributes => [qw/ tissue /],
            individual_attributes => [qw/ race gender /],
        },
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

    my $name_regexp_string = $config->{name_regexp};
    die 'No name regexp!' if not $name_regexp_string;
    die 'Invalid name regexp! '.$name_regexp_string if $name_regexp_string !~ /^\($nomenclature/;
    my $name_regexp =  qr|^$name_regexp_string$|; 

    my %properties;
    for my $type (qw/ sample individual /) {
        my $key_name = $type.'_attributes';
        next if not $config->{$key_name};
        my %attributes;
        if ( ref $config->{$key_name} eq 'ARRAY' ) {
            %attributes = map { $_ => {} } @{$config->{$key_name}};
        }
        else { # hash
            %attributes = %{$config->{$key_name}};
        }
        my %type_properties = _get_properties_for_import_command_from_entity($type, %attributes);
        %properties = ( %properties, %type_properties );
        $properties{'_'.$type.'_attribute_names'} = { is => 'ARRAY', is_constant => 1, value => [ keys %attributes ], };
    }

    $import_class_names{$class_name} = UR::Object::Type->define(
        class_name => $class_name,
        is => 'Genome::Sample::Command::Import::Base',
        has => [
            %properties,
            name_regexp => { is_constant => 1, value => $name_regexp, },
            nomenclature => { is_constant => 1, },
        ],
        doc => "import $nomenclature samples",
    );

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
    my ($entity, %attributes) = @_;

    my $class_name = 'Genome::'.ucfirst($entity);
    my $meta = $class_name->__meta__;
    my %properties;
    for my $name ( keys %attributes ) {
        my $property = $meta->property_meta_for_name($name);
        $properties{$name} = $attributes{$name};
        $properties{$name}->{is} = 'Text' if not defined $properties{$name}->{is};
        $properties{$name}->{doc} = "The value of '".join(' ', split('_', $name))."' for the $entity." if not defined $properties{$name}->{doc};
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

