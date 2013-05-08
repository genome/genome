package Genome::Model::Tools::Sx::Limit::Base;

use strict;
use warnings;

use Genome;

use Regexp::Common;

class Genome::Model::Tools::Sx::Limit::Base {
    is  => 'Genome::Model::Tools::Sx::Base',
    is_abstract => 1,
    has_optional => [
        incoming_sequences => {
            is => 'Text',
            doc => 'To select random sequences, the incoming value [depening on how sequences are being limited] needs to be known. This value can be given as an integer or SX metrics file.',
        },
    ],
};

sub execute {
    my $self = shift;

    if ( $self->incoming_sequences ) {
        return $self->_execute_selecting_random_sequences;
    }
    else {
        return $self->_execute;
    }
}

sub _execute {
    my $self = shift;

    my $init = $self->_init;
    return if not $init;

    my $reader = $self->_reader;
    my $writer = $self->_writer;

    my @limiters = $self->_create_limiters;
    return if not @limiters;

    READER: while ( my $seqs = $reader->read ) {
        $writer->write($seqs);
        for my $limiter ( @limiters ) {
            last READER unless $limiter->($seqs); # returns 0 when done
        }
    }

    return 1;
}

sub _execute_selecting_random_sequences {
    my $self = shift;

    my ($keep_pct, $randomizier) = $self->_init_randomizer;
    return if not $randomizier; # error
    return $self->_execute if $keep_pct >= 1; # randomizing not needed

    my $init = $self->_init;
    return if not $init;

    my $reader = $self->_reader;
    my $writer = $self->_writer;

    my @limiters = $self->_create_limiters;
    return if not @limiters;

    READER: while ( my $seqs = $reader->read ) {
        next if not $randomizier->($seqs);
        $writer->write($seqs);
        for my $limiter ( @limiters ) {
            last READER unless $limiter->($seqs); # returns 0 when done
        }
    }

    return 1;
}

sub _init_randomizer {
    my $self = shift;

    my $metric_name = lc( $self->class );
    $metric_name =~ s/genome::model::tools::sx::limit::by//; # count or bases
    my $limit = $self->$metric_name;

    my $incoming_sequences = $self->incoming_sequences;
    if ( -f $incoming_sequences ) {
        my $metrics = Genome::Model::Tools::Sx::Metrics::Basic->from_file($incoming_sequences);
        return if not $metrics;
        if ( not defined $metrics->$metric_name ) {
            $self->error_message("Failed to find metric ($metric_name) in metrics file: $incoming_sequences");
            return;
        }
        $incoming_sequences = $metrics->$metric_name;
    }

    if ( $incoming_sequences !~ /^$RE{num}{int}$/ ) {
        $self->error_message('Unknown value for incoming sequences. Expected an integer [as param or from metrics file], but was given: '.$incoming_sequences);
        return;
    }

    my $keep_pct = sprintf('%0.2f', ( $limit / $incoming_sequences ));
    return (
        $keep_pct,
        sub{ return ( sprintf('%0.2f', rand()) <= $keep_pct ) ? 1 : 0; },
    );
}

1;

