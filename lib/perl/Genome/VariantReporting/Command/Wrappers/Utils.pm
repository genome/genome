package Genome::VariantReporting::Command::Wrappers::Utils;

use warnings;
use strict;

use Genome;
use List::MoreUtils qw(uniq);

class Genome::VariantReporting::Command::Wrappers::Utils { };

sub get_library_name_labels {
    my ($self, $category, $sample, $builds) = @_;
    my %labels;
    my $counter = 1;
    my @libraries = uniq(map {$_->library}
                            grep {$_->library->sample eq $sample}
                                map {$_->instrument_data} @$builds);
    for my $library (sort {$a->name cmp $b->name} @libraries) {
        $labels{$library->name} = sprintf('%s-Library%d(%s)',
            ucfirst($category),
            $counter++,
            $library->name,
        );
    }
    return %labels;
}

1;
