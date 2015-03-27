package Genome::Utility::File::Comparer;

use strict;
use warnings;
use Genome;
use File::Compare;

class Genome::Utility::File::Comparer {
    has => [
        special_compare_functions => {
            is => 'ARRAY', #Does this work to make sure the pairs come back in order?
            doc => 'Keys are regex that match file names.  Values are code refs',
        },
    ],
};

sub compare {
    my $self = shift;
    my ($file1, $file2) = @_;
    my %special_functions = @{$self->special_compare_functions};
    while (my ($regex, $function) = each %special_functions) {
        if ($file1 =~ /$regex/) {
            return $function->($file1, $file2);
        }
    }
    return regular_compare($file1, $file2);
}

sub regular_compare {
    my ($target, $other_target) = @_;
    return File::Compare::compare($target, $other_target);
}

1;

