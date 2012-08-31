package Genome::Model::Tools::Sx::BedWriter;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::Sx::BedWriter {
    is => 'Genome::Model::Tools::Sx::SeqWriter',
};

my $max = 10_000_000;
sub write {
    my ($self, $seq) = @_;

    Carp::confess('No seq to write in bed format!') if not $seq;

    my $length = length $seq->{seq};
    my $cnt = $length / $max;
    for ( my $i = 0; $i <= $cnt; $i++ ) {
        my $end =  ($i + 1) * $max;
        if ( $end > $length ) { $end = $length; }
        $self->{_file}->print( join( "\t", $seq->{id}, ( $i * $max), $end, ( $length > $max ? $seq->{id}.'part'.($i + 1): $seq->{id} ))."\n" );
    }

    return 1;
}

1;

