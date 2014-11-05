package Genome::VariantReporting::Framework::Utility;

use strict;
use warnings;
use Genome;

sub get_missing_errors {
    my ($name, $params, $needed, $type_of_thing, $thing_that_needs) = @_;
    my $have = Set::Scalar->new(keys %{$params});
    my @errors;
    unless($needed->is_equal($have)) {
        if (my $still_needed = $needed - $have) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [$still_needed->members],
                desc => sprintf("$type_of_thing required by $thing_that_needs (%s) but not provided: (%s)", 
                    $name, join(",", $still_needed->members)),
            );
        }
    }

    return @errors;
}


1;

