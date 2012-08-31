package Genome::Model::Tools::Sx::Dedup::ByName;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Dedup::ByName {
    is  => 'Genome::Model::Tools::Sx::Dedup::Base',
};

sub help_brief { return 'Remove deplicate sequences by name'; }
sub help_detail { return 'Remove deplucicates, leaving only one sequence per name.'; }

sub execute {
    my $self = shift;

    my $init = $self->_init;
    return if not $init;

    my $reader = $self->_reader;
    my $writer = $self->_writer;

    my %seen_names;
    while ( my $seqs = $reader->read ) {
        my @seqs = grep { not $seen_names{ $_->{id} }++ } @$seqs;
        next if not @seqs;
        $writer->write(\@seqs);
    }

    return 1;
}

1;

