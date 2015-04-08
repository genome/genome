package Genome::Utility::Inputs;

use strict;
use warnings FATAL => 'all';
use Genome;
use Params::Validate qw(validate_pos :types);
use Scalar::Util qw();

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    encode
    decode
);


# Given a hashref of inputs, some of which may be UR Datasource backed objects,
# encode them for serialization
sub encode {
    my ($inputs) = validate_pos(@_, {type=>HASHREF});

    my %encoded_inputs = %$inputs;
    while (my ($name, $value) = each %$inputs) {
        if (Scalar::Util::blessed($value)) {
            $encoded_inputs{$name} = convert_obj_to_hash($value);
        } elsif (ref($value) eq 'ARRAY' &&
                 scalar(@{$value}) &&
                 Scalar::Util::blessed($value->[0])) {
            $encoded_inputs{$name} = [map {convert_obj_to_hash($_)} @{$value}];
        }
    }
    return \%encoded_inputs;
}

sub convert_obj_to_hash {
    my $obj = shift;

    return {
        class => $obj->class,
        id => $obj->id,
    };
}


# Given a hashref of inputs, some of which may be encodings for UR Datasource
# backed objects, decode them, instantiating the objects
sub decode {
    my ($inputs) = validate_pos(@_, {type=>HASHREF});

    my %decoded_inputs = %$inputs;
    while (my ($name, $value) = each %$inputs) {
        if (ref($value) eq 'HASH') {
            $decoded_inputs{$name} = convert_hash_to_obj($value);
        } elsif (ref($value) eq 'ARRAY' &&
                 scalar(@{$value}) &&
                 ref($value->[0]) eq 'HASH') {
            $decoded_inputs{$name} = [map {convert_hash_to_obj($_)} @{$value}];
        }
    }

    return \%decoded_inputs;
}

sub convert_hash_to_obj {
    my $hash = shift;

    my $class = $hash->{'class'};
    my $obj = $class->get($hash->{'id'});
    unless (defined($obj)) {
        die sprintf("Couldn't convert hash to class: %s", pp($hash));
    }
    return $obj;
}


1;
