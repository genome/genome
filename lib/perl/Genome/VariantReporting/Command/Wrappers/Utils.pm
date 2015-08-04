package Genome::VariantReporting::Command::Wrappers::Utils;

use warnings;
use strict;

use Genome;

class Genome::VariantReporting::Command::Wrappers::Utils { };

sub get_library_name_labels {
    my ($self, $category) = @_;

    my %labels;
    my $counter = 1;

    my $accessor = sprintf('%s_sample', $category);
    for my $library ($self->$accessor->libraries) {
        $labels{$library->name} = sprintf('%s-Library%d(%s)',
            ucfirst($category),
            $counter++,
            $library->name,
        );
    }
    return %labels;
}

1;
