package Genome::Model::Tools::Sx::Dedup::BySequence;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::Dedup::BySequence {
    is  => 'Genome::Model::Tools::Sx::Dedup::Base',
};

sub help_brief { return 'Remove deplicate identical'; }
sub help_detail { return 'Remove deplucicates, leaving only one identical sequence.'; }

sub execute {
    my $self = shift;

    my $init = $self->_init;
    return if not $init;

    my $reader = $self->_reader;
    my $writer = $self->_writer;

    my %seen_sequences;
    while ( my $seqs = $reader->read ) {
	my $combined_sequence = join('', map { $_->{seq} } @$seqs);
	next if $seen_sequences{$combined_sequence};
	$seen_sequences{$combined_sequence} = $combined_sequence;
	$writer->write( $seqs );
    }

    return 1;
}

1;

