package Genome::Role::Comparable;

use strict;
use warnings;
use UR::Role;

role Genome::Role::Comparable {
    requires => ['compare_output'],
};

sub special_compare_functions {
    return ();
}

my %comparators;
sub file_comparer {
    my $self = shift;
    unless (defined $comparators{$self->id}) {
        my $comparer = Genome::Utility::File::Comparer->create(
            special_compare_functions => [$self->special_compare_functions],
        );
        $comparators{$self->id} = $comparer;
    }
    return $comparators{$self->id};
}
1;

