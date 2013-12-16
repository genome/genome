package Genome::File::Vep::Command::ListGenes;

use Genome;
use Genome::File::Vep::Reader;

use strict;
use warnings;

class Genome::File::Vep::Command::ListGenes {
    is => "Command::V2",
    has_input => [
        input_file => {
            is => "FilePath",
            doc => "A vep annotation input file to look for gene names in",
        }
    ],

    has_output => [
        genes => {
            is_many => 1,
            is => "Text",
            doc => "The list of gene names found in the input file",
        }
    ],
};

sub execute {
    my $self = shift;
    my $reader = new Genome::File::Vep::Reader($self->input_file);
    my %genes;
    while (my $entry = $reader->next) {
        $genes{$entry->{gene}} = 1 if exists $entry->{gene} and defined $entry->{gene};
    }

    $self->genes([keys %genes]);
    return 1;
};

1;
