package Genome::ConfigValidator::UR;

use strict;
use warnings;

use Genome::Carp qw(croakf);

sub validate {
    my ($value, $spec) = @_;

    unless ($spec->type->isa('UR::Value')) {
        croakf('%s is not a UR::Value', $spec->type);
    }

    my $value_obj = $spec->type->get($value);
    my @errors = $value_obj->__errors__;
    return unless @errors;
    return sprintf('a %s', $spec->type);
}

1;
