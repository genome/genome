package Genome::Command::Copy;

use strict;
use warnings FATAL => 'all';

use Genome;

class Genome::Command::Copy {
    is => 'Command::V2',
    is_abstract => 1,
    has => [
        source => {
            is => 'UR::Object',
        },
        changes => {
            is => 'Text',
            is_many => 1,
            doc => 'comma-separated list of changes',
        },
    ],
};

sub help_detail {
    return join('',
        q(Any non-delegated, non-ID properties may be specified with an operator and a value.),
        q( Valid operators are '=', '+=', '-=', and '.='; function is same as in Perl.),
        qq(\n\nFor example:\n\n),
        qq(   --changes "name.=-RT101912,foo=bar"\n\n),
        q( A value of 'undef' may be used to pass a Perl undef as the value.  Either `foo=` or `foo=''` can be used to set the value to an empty string.),
    );
}

sub execute {
    my $self = shift;

    my $source_type = $self->source->__meta__;
    my @copyable_properties = grep { !$_->is_delegated && !$_->is_id } $source_type->properties;
    my %params = map {
        my $name = $_->property_name;
        my $value = $self->source->$name;
        defined $value ? ($name => $value) : ();
    } @copyable_properties;

    for my $change ($self->changes) {
        my ($key, $op, $value) = $change =~ /^(.+?)(=|\+=|\-=|\.=)(.*)$/;
        unless ($key && $op) {
            $self->error_message("invalid change: $change");
            return;
        }

        if ($value eq 'undef') {
            $value = undef;
        }

        my $property = grep { $_->property_name eq $key } @copyable_properties;

        unless ($property) {
            $self->error_message("unrecognized property: $key");
            return;
        }

        if ($op eq '=') {
            $params{$key} = $value;
        }
        elsif ($op eq '+=') {
            $params{$key} += $value;
        }
        elsif ($op eq '-=') {
            $params{$key} -= $value;
        }
        elsif ($op eq '.=') {
            $params{$key} .= $value;
        }
    }

    my $source_class = $source_type->class_name;
    my $copy = $source_class->create(%params);
    $self->status_message('Created new %s with ID %s', $source_class, $copy->id);

    return 1;
}
