package Genome::DruggableGene::Command::GeneNameReport;

use strict;
use warnings;

use Genome;

class Genome::DruggableGene::Command::GeneNameReport {
    is => 'Command::Tree',
};

sub help_brief {
    "work with gene-name-reports"
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
 genome druggable-gene gene-name-report ...
EOS
}

sub help_detail {
    return <<EOS
A collection of commands to interact with gene-names.
EOS
}

1;
