package Genome::VariantReporting::Framework::Component::RuntimeTranslations;

use strict;
use warnings FATAL => 'all';
use UR;
use YAML;
use JSON;
use Params::Validate qw(validate validate_pos :types);

my $_JSON_CODEC = new JSON->allow_nonref;

class Genome::VariantReporting::Framework::Component::RuntimeTranslations {
    has => [
        translations => {
            is => 'HASH',
            default => {},
        },
    ],
};

sub as_json {
    my $self = shift;

    return $_JSON_CODEC->canonical->encode($self->translations);
}

sub create_from_json {
    my $class = shift;
    my $json = shift;

    my $hashref = $_JSON_CODEC->decode($json);
    return $class->create(translations => $hashref);
}

sub write_to_file {
    my $self = shift;
    my $filename = shift;

    YAML::DumpFile($filename, $self->translations);
}

sub create_from_file {
    my $class = shift;
    my $file = shift;

    Genome::Sys->validate_file_for_reading($file);
    my ($hashref, undef, undef) = YAML::LoadFile($file);

    return $class->create(translations => $hashref);
}


1;
