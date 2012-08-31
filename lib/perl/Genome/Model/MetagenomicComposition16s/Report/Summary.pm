package Genome::Model::MetagenomicComposition16s::Report::Summary;

use strict;
use warnings;

use Genome;

use Carp 'confess';
use Data::Dumper 'Dumper';
use Regexp::Common;

class Genome::Model::MetagenomicComposition16s::Report::Summary {
    is => 'Genome::Model::MetagenomicComposition16s::Report',
};

#< Generator >#
sub description {
    return 'Summary Report for '.
    Genome::Utility::Text::capitalize_words( $_[0]->build->description );
}

sub _add_to_report_xml {
    my $self = shift;

    $self->_create_metrics;

    my @amplicon_sets = $self->build->amplicon_sets;
    Carp::confess('No amplicon sets for '.$self->build) if not @amplicon_sets;

    for my $amplicon_set ( @amplicon_sets ) {
        next if not $amplicon_set->amplicon_iterator;
        while ( my $amplicon = $amplicon_set->next_amplicon ) {
            $self->_add_amplicon($amplicon);
        }
    }

    # Summary Stats 
    my $summary_stats = $self->get_summary_stats
        or return;
    $self->_add_dataset(
        name => 'stats',
        row_name => 'stat',
        headers => $summary_stats->{headers},
        rows => [ $summary_stats->{stats} ],
    ) or return;

    return 1;
}

sub _create_metrics {
    my $self = shift;

    $self->{_metrix} = {
        lengths => [],
    };

    return 1;
}

sub _add_amplicon {
    my ($self, $amplicon) = @_;

    return 1 if not $amplicon->{classification}; # ok
    push @{$self->{_metrix}->{lengths}}, length $amplicon->{seq}->{seq};

    return 1;
}

sub get_summary_stats {
    my $self = shift;

    my $build = $self->build;
    my $attempted = $build->amplicons_attempted || 0;
    my $processed = $build->amplicons_processed;
    my $processed_success = $build->amplicons_processed_success;
    if ( not defined $processed or $processed == 0 ) {
        return {
            headers => [qw/ amplicons-processed amplicons-attempted amplicons-success /],
            stats => [ 0, $attempted, 0 ] 
        };
    }

    my $sum = sub{
        my $total = 0;
        for ( @_ ) { $total += $_; }
        return $total;
    };

    my @lengths = sort { $a <=> $b } @{ $self->{_metrix}->{lengths} };
    my $length = $sum->(@lengths);

    my %totals = (
        # Amplicons
        'amplicons-attempted' => $attempted,
        'amplicons-processed' => $processed,
        'amplicons-processed-success' => $processed_success,
        'amplicons-classified' => $build->amplicons_classified,
        'amplicons-classified-success' => $build->amplicons_classified_success,
        'amplicons-classification-error' => $build->amplicons_classification_error,
        # Lengths
        'length-minimum' => $lengths[0],
        'length-maximum' => $lengths[$#lengths],
        'length-median' => $lengths[( $#lengths / 2 )],
        'length-average' => sprintf(
            '%.0f',
            $length / $processed,
        ),
        # Reads
        reads => $build->reads_attempted,
        'reads-processed' => $build->reads_processed,
        'reads-processed-success' => $build->reads_processed_success,
    );

    return {
        headers => [ sort { $a cmp $b } keys %totals ],
        stats => [ map { $totals{$_} } sort { $a cmp $b } keys %totals ],
    }
}

1;

