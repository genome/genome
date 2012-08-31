package Genome::Model::Tools::Sx::Rename::IlluminaToPcap;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Tools::Sx::Rename::IlluminaToPcap {
    is  => 'Genome::Model::Tools::Sx::Base',
};

sub help_brief {
    return <<HELP 
    Rename sequences to pcap style b1/g1
HELP
}

sub help_detail {
    return <<HELP
HELP
}

my @match_and_replace = (
    # use qr{} for speed boost
    [ qr{#.*/1$}, '.b1' ],
    [ qr{#.*/2$}, '.g1' ],
);
sub _eval_seqs {
    my ($self, $seqs) = @_;

    for my $seq ( @$seqs ) { 
        MnR: for my $match_and_replace ( @match_and_replace ) {
            $seq->{id} =~ s/$match_and_replace->[0]/$match_and_replace->[1]/g 
                and last MnR;
        }
    }

    return 1;
}

1;

