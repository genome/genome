package Genome::Sample::Command::Import;

use strict;
use warnings;
use Genome;

class Genome::Sample::Command::Import {
    is => 'Command::Tree',
    doc => 'commands for importing samples',
};

my %import_namespaces;
_create_import_commands();

sub namespace_for_nomenclature {
    my ($self, $nomenclature) = @_;
    Carp::confess('No nomenclature to get namespace!') if not $nomenclature;
    my ($config) = grep { $_->{nomenclature} eq $nomenclature } values %import_namespaces;
    return $config->{namespace} if $config;
}

sub importer_class_name_for_namespace {
    my ($self, $namespace) = @_;

    Carp::confess('No namespace to get property names for importer!') if not $namespace;

    my $config = $import_namespaces{$namespace};
    if ( not $config ) {
        Carp::confess('No config found for namespace to get property names for importer!');
    }

    return $config->{importer_class_name};
}

sub _create_import_commands {
    my @configs = _load_import_configs();
    for my $config ( @configs ) { 
        _create_import_command_for_config($config); 
    }
    return 1;
}

sub _load_import_configs {
    #TODO
    # use the upcoming config API!
    # mv tests to config
    return (
        {
            nomenclature => 'ATCC',
            name_regexp => '(ATCC\-[\w\d]+\-[\w\d]+)(\-[\w\d]+)?',
            sample_attributes => {
                age => {},
                disease => { is_optional => 1, },
                organ_name => {},
                sample_common_name => {
                    calculate_from => [qw/ disease /],
                    calculate => sub{
                        my $disease = shift;
                        return $disease ? 'tumor' : 'normal';
                    },
                },
            },
            individual_attributes => {
                ethnicity => {},
                gender => { valid_values => [qw/ male female /], },
                individual_common_name => {
                    calculate_from => [qw/ _individual_name /],
                    calculate => sub{
                        my $_individual_name = shift;
                        $_individual_name =~ s/^ATCC\-//;
                        return $_individual_name;
                    },
                },
            },
        },
        {
            nomenclature => 'dbGaP',
            name_regexp => '(dbGaP\-\d+)\-\d+',
            sample_attributes => [qw/ tissue /],
            individual_attributes => [qw/ race gender /],
        },
        {
            nomenclature => 'EMBL-EBI',
            name_regexp => '(EMBL\-[\w\d_]+)\-.+',
            sample_attributes => {
                age => {},
                tissue_label => { doc => 'Tissue from where the sample was taken.', },
                tissue_desc => { calculate_from => [qw/ tissue_label /], calculate => sub{ return $_[0]; }, },
            },
            individual_attributes => {
                ethnicity => {},
                gender => { valid_values => [qw/ male female unknown /], },
            },
        },
        {
            nomenclature => 'METAHIT',
            name_regexp => '(METAHIT\-[\w\d]+)\-\d+',
            sample_attributes => [qw/ age body_mass_index tissue_label /],
            individual_attributes => [qw/ gender /],
        },
        {
            nomenclature => 'LSR',
            name_regexp => '(LSR\-ER\-[\w\d]+)\-.+',
            sample_attributes => {
                extraction_type => {
                    default_value => "rna",
                },
            },
            library_attributes => {
                transcript_strand => {
                    default_value => "unstranded",
                },
            },
       },
       {
            nomenclature => 'SRA',
            name_regexp => '(SRA\-\w+\-[\w\d]+\-[\w\d]+)\-.+',
            sample_attributes => {
                extraction_type => {
                    default_value => "genomic dna",
                },
            },
            library_attributes => {
                transcript_strand => {
                    default_value => "unstranded",
                },
            },
       },
       {
            nomenclature => 'BRC46',
            name_regexp => '(BRC46\-(BRC|CSB)[0-9]+)\_.+',
            sample_attributes => {
                extraction_type => {
                    default_value => "rna",
                },
            },
        },
        {
            nomenclature => 'TCGA',
            name_regexp => '(TCGA\-[\w\d]+\-[\w\d]+)\-[\w\d]+\-[\w\d]+\-[\w\d]+\-[\w\d]+',
            sample_attributes => {
                extraction_type => {
                    calculate_from => [qw/ name /],
                    calculate => sub{
                        my $name = shift;
                        my %extraction_types = (
                            D => 'genomic dna',
                            G => 'ipr product',
                            R => 'rna',
                            T => 'total rna',
                            W => 'ipr product',
                            X => 'ipr product',
                        );
                        my @tokens = split('-', $name);
                        my ($extraction_code) = $tokens[4] =~ /(\w)$/;
                        if ( not $extraction_code ) {
                            die "Cannot get extraction code from name part: $tokens[4]";
                            return;
                        }
                        if ( not $extraction_types{$extraction_code} ) {
                            die "Invalid extraction code ($extraction_code) found in name ($name)";
                            return;
                        }
                        return $extraction_types{$extraction_code};
                    },
                },
                common_name => {
                    calculate_from => [qw/ name /],
                    calculate => sub{
                        my $name = shift;
                        my @tokens = split('-', $name);
                        my ($sample_type) = $tokens[3] =~ /^(\d+)/;
                        if ( not $sample_type ) {
                            die "Cannot get extraction code from name part: $tokens[3]";
                            return;
                        }
                        my $common_name;
                        if ( $sample_type >= 1 and $sample_type <= 9 ) {
                            $common_name = 'tumor';
                        }
                        elsif ( $sample_type >= 10 and $sample_type <= 19 ) {
                            $common_name = 'normal';
                        }
                        else {
                            $common_name = 'unknown';
                        }
                        return $common_name;
                    },
                },
            },
        },
        {
            nomenclature => 'DREAM',
            name_regexp => '(DREAM\-[\w\d]+)\-[\w\d]+',
            sample_attributes => {
                extraction_type => {
                    default_value => "genomic dna",
                },
                common_name => {
                },
            }
        },
        {
            nomenclature => 'WHIM',
            name_regexp => '(WHIM[\d]+)\-[\w\d]+(\-[\w\d]+)*',
        },
    );
}

sub _create_import_command_for_config {
    my $config = shift;

    my $nomenclature = $config->{nomenclature};
    die 'No nomenclature in sample import command config!' if not $nomenclature;
    $config->{namespace} //= join('', map { ucfirst(lc($_)) } split('-', $nomenclature));
    my $class_name = 'Genome::Sample::Command::Import::'.$config->{namespace};
    return 1 if exists $import_namespaces{ $config->{namespace} };
    $import_namespaces{ $config->{namespace} } = $config;
    $config->{importer_class_name} = $class_name;

    my $name_regexp_string = $config->{name_regexp};
    die 'No name regexp in sample import command config!' if not $name_regexp_string;
    #die 'Invalid name regexp in sample import command config! '.$name_regexp_string if $name_regexp_string !~ /^\($nomenclature/;
    die 'Invalid name regexp (no parens to capture individual name) in sample import command config! '.$name_regexp_string if $name_regexp_string !~ /^\(.+\)/;
    my $name_regexp =  qr|^$name_regexp_string$|; 

    my @command_properties;
    for my $entity (qw/ sample individual /) {
        my $key_name = $entity.'_attributes';
        next if not $config->{$key_name};
        my %attributes;
        if ( ref $config->{$key_name} eq 'ARRAY' ) { # convert to HASH
            $config->{$key_name} = { map { $_ => {} } @{$config->{$key_name}} };
        }
        my $properties = $config->{$key_name};
        _add_property_meta_from_entity_to_importer_properties($entity, $properties);
        $properties->{'_'.$entity.'_attribute_names'} = { is => 'ARRAY', is_constant => 1, value => [ sort keys %$properties ], };
        push @command_properties, $properties;
    }

    my $importer_class_meta = UR::Object::Type->define(
        class_name => $class_name,
        is => 'Genome::Sample::Command::Import::Base',
        has => {
            map({%$_ } @command_properties),
            name_regexp => { is_constant => 1, value => $name_regexp, },
            nomenclature => { is_constant => 1, },
        },
        doc => "import $nomenclature samples",
    );

    my $nomenclature_property = $importer_class_meta->property_meta_for_name('nomenclature');
    $nomenclature_property->default_value($nomenclature);

    if ( $config->{minimum_unique_source_name_parts} ) {
        die 'Invalid minimum_unique_source_name_parts! '.$config->{minimum_unique_source_name_parts} if $config->{minimum_unique_source_name_parts} > 1;
        my $minimum_unique_source_name_parts = $importer_class_meta->property_meta_for_name('_minimum_unique_source_name_parts');
        $minimum_unique_source_name_parts->default_value( $config->{minimum_unique_source_name_parts} // 2 );
    }

    return 1;
}

sub _add_property_meta_from_entity_to_importer_properties {
    my ($entity, $properties) = @_;

    Carp::confess('No entity!') if not $entity;
    Carp::confess('No property names!') if not $properties;

    my $class_name = 'Genome::'.ucfirst($entity);
    my $meta = $class_name->__meta__;
    for my $name ( keys %$properties ) {
        my $property = $properties->{$name};
        my $property_from_entity = $meta->property_meta_for_name($name);
        if ( $property_from_entity ) {
            for my $from_entity_property_name (qw/ doc is valid_values /) {
                next if not defined $property_from_entity->{$from_entity_property_name};
                $property->{$from_entity_property_name} = $property_from_entity->{$from_entity_property_name};
            }
        }
        else { 
            $property->{is} = 'Text' if not defined $property->{is};
            $property->{doc} = "The value of '".join(' ', split('_', $name))."' for the $entity." if not defined $property->{doc};
        }
        $property->{is_optional} = 1 if $property->{calculate};
    }

    return 1;
}

# Overload sub command classes to return these in memory ones, plus the existing ones
sub sub_command_names {
    my $class = shift;
    my %sub_command_names = map { $_ => 1 } ($class->SUPER::sub_command_names, keys %import_namespaces);
    return keys %sub_command_names;
}

sub sub_command_classes {
    my $class = shift;
    my %sub_command_classes = map { 
        $_ => 1 
    } ($class->SUPER::sub_command_classes, map({ $_->{importer_class_name} } values %import_namespaces));
    return keys %sub_command_classes;
}

1;

