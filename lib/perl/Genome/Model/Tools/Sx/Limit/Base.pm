package Genome::Model::Tools::Sx::Limit::Base;

use strict;
use warnings;

use Genome;

use Regexp::Common;

class Genome::Model::Tools::Sx::Limit::Base {
    is  => 'Genome::Model::Tools::Sx::Base',
    is_abstract => 1,
    has => [
        select_random_sequences => {
            is_optional => 1,
            doc => 'Select random sequences when limiting. Give an SX metrics file or value. Depending on how the sequences are being limited (bases, count), the incoming value is required. The appropriate value (bases, count) can be loaded from a given metrics file. If no metrics file is available, give the incoming bases when limiting by bases. Or give the count if limiting by count.',
        },
    ],
};

sub execute {
    my $self = shift;

    if ( $self->select_random_sequences ) {
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

    my ($randomizier, $keep_pct) = $self->_init_randomizer;
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

    my $select_random_sequences = $self->select_random_sequences;
    my $metric_value;
    if ( -f $select_random_sequences ) {
        my $metrics = Genome::Model::Tools::Sx::Metrics::Basic->from_file($select_random_sequences);
        return if not $metrics;
        $metric_value = $metrics->$metric_name;
        if ( not $metric_value ) {
            $self->error_message("Failed to find metric ($metric_name) in metrics file: $select_random_sequences");
            return;
        }
    }
    elsif ( $select_random_sequences =~ /^$RE{num}{int}$/ ) {
        $metric_value = $select_random_sequences;
    }
    else {
        $self->error_message('Unknown value to select random sequences. Expect a metrics file or an integer but was given: '.$select_random_sequences);
        return;
    }

    my $keep_pct = sprintf('%0.2f', ( $metric_value / $self->$metric_name ));
    return (
        sub{
            my $rand = sprintf('%0.2f', rand());
            return ( $rand <= $keep_pct ) ? 1 : 0;
        },
        $keep_pct,
    );
}

1;

