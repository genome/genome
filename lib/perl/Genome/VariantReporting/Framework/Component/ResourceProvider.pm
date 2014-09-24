package Genome::VariantReporting::Framework::Component::ResourceProvider;

use strict;
use warnings FATAL => 'all';
use UR;
use YAML;
use JSON;
use Data::Dump qw(pp);
use Params::Validate qw(validate validate_pos :types);

my $_JSON_CODEC = new JSON->allow_nonref;

class Genome::VariantReporting::Framework::Component::ResourceProvider {
    has => [
        attributes => {
            is => 'HASH',
            default => {},
        },
    ],
};

sub set_attribute {
    my ($self, $name, $value) = validate_pos(@_, 1, 1, 1);
    $self->attributes->{$name} = $value;
}

sub get_attribute {
    my ($self, $name) = validate_pos(@_, 1, 1);

    if (exists $self->attributes->{$name}) {
        return $self->attributes->{$name};
    } else {
        die "Attempted to get non-existing attribute ($name) from resource-provider, available attributes are ".
          pp(keys %{$self->attributes});
    }
}

sub as_json {
    my $self = shift;

    return $_JSON_CODEC->canonical->encode($self->attributes);
}

sub create_from_json {
    my $class = shift;
    my $json = shift;

    my $hashref = $_JSON_CODEC->decode($json);
    return $class->create(attributes => $hashref);
}

sub write_to_file {
    my $self = shift;
    my $filename = shift;

    YAML::DumpFile($filename, $self->attributes);
}

sub create_from_file {
    my $class = shift;
    my $file = shift;

    Genome::Sys->validate_file_for_reading($file);
    my ($hashref, undef, undef) = YAML::LoadFile($file);

    return $class->create(attributes => $hashref);
}


1;
