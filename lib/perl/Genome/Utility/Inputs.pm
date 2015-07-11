package Genome::Utility::Inputs;

use strict;
use warnings FATAL => 'all';
use Data::Dump 'pp';
use Try::Tiny qw(try catch);
require Data::Transform::ExplicitMetadata;

use Genome;
use Params::Validate qw(validate_pos :types);
use Scalar::Util qw();

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    encode
    decode
);

# Given a hashref of inputs, some of which may be objects,
# encode them to be safe for JSON serialization
sub encode {
    my ($inputs) = validate_pos(@_, {type=>HASHREF});

    my %encoded_inputs = %$inputs;
    while (my ($name, $value) = each %$inputs) {
        my $encoder = _get_encoder($value);
        if (defined $encoder) {
            $encoded_inputs{$name} = $encoder->($value);
        }
    }

    return \%encoded_inputs;
}

sub _get_encoder {
    my $value = shift;

    if (ref($value) eq 'HASH') {
        return \&encode;
    }

    if (ref($value) eq 'ARRAY') {
        return \&_encode_array;
    }

    if (Scalar::Util::blessed($value)) {
        return \&_encode_object;
    }

    # no encoding to be done
    return;
}

sub _encode_array {
    my $array = shift;

    my @encoded_array;
    for my $element (@$array) {
        my $encoder = _get_encoder($element);
        if (defined $encoder) {
            push @encoded_array, $encoder->($element);
        } else {
            push @encoded_array, $element;
        }
    }

    return \@encoded_array;
}

sub _encode_object {
    my $obj = shift;

    if ($obj->isa("UR::Object")) {
        return {
            __decoder_type__ => 'UR object',
            class => $obj->class,
            id => $obj->id,
        };
    } else {
        my $body = Data::Transform::ExplicitMetadata::encode($obj);
        return {
            __decoder_type__ => 'generic object',
            %$body,
        }
    }
}


# Given a hashref of inputs, some of which may have been encoded,
# decode them, instantiating any encoded objects.
sub decode {
    my ($inputs) = validate_pos(@_, {type=>HASHREF});

    my %decoded_inputs = %$inputs;
    while (my ($name, $value) = each %$inputs) {
        my $decoder = _get_decoder($value);
        if (defined $decoder) {
            $decoded_inputs{$name} = $decoder->($value);
        }
    }

    return \%decoded_inputs;
}

sub _get_decoder {
    my $value = shift;

    if (ref($value) eq 'HASH') {
        if (exists $value->{__decoder_type__}) {
            return _get_decoder_of_type($value->{__decoder_type__});
        } else {
            return \&decode;
        }
    }

    if (ref($value) eq 'ARRAY') {
        return \&_decode_array;
    }

    # no encoding to be done
    return;
}

sub _get_decoder_of_type {
    my $type = shift;

    my %TYPES = (
        'UR object' => \&_decode_ur_object,
        'generic object' => \&_decode_generic_object,
    );

    if (exists $TYPES{$type}) {
        return $TYPES{$type};
    } else {
        die sprintf("Couldn't determine decoder for type (%s)", $type);
    }
}

sub _decode_ur_object {
    my $hash = shift;

    my $class = $hash->{'class'};
    my $obj = try {
        $class->get($hash->{'id'}) || die sprintf(
            "Nothing returned from '%s->get(%s)'", $class, pp($hash->{'id'}));
    } catch {
        my $error = $_ || "Unknown or unspecified error";
        die sprintf("Failed to convert hash to ur object: %s", $error);
    };

    return $obj;
}

sub _decode_generic_object {
    my $hash = shift;

    my %object_metadata = %$hash;
    delete $object_metadata{__decoder_type__};

    my $obj = try {
        Data::Transform::ExplicitMetadata::decode(\%object_metadata) ||
            die sprintf("Nothing returned from " .
                "Data::Transform::ExplicitMetadata::decode called with (%s)",
                pp(\%object_metadata));
    } catch {
        my $error = $_ || "Unknown or unspecified error";
        die sprintf("Failed to convert hash to generic object: %s", $error);
    };

    return $obj;
}

sub _decode_array {
    my $array = shift;

    my @decoded_array;
    for my $element (@$array) {
        my $decoder = _get_decoder($element);
        if (defined $decoder) {
            push @decoded_array, $decoder->($element);
        } else {
            push @decoded_array, $element;
        }
    }

    return \@decoded_array;
}

1;
