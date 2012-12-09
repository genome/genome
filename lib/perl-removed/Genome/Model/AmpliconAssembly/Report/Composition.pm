package Genome::Model::AmpliconAssembly::Report::Composition;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require Genome::Utility::MetagenomicClassifier::PopulationComposition;

class Genome::Model::AmpliconAssembly::Report::Composition {
    is => 'Genome::Model::Report',
    has => [
    description => {
        calculate => q| 
        return sprintf(
            'Metagenomic Compposition Report for Amplicon Assembly (Name <%s> Build Id <%s>)',
            $self->model_name,
            $self->build_id,
        );
        |,
    },
    confidence_threshold => {
        is => 'Number',
        default_value => 0.8,
        doc => 'The threshold to match or exceed (>=) to consider a taxonic rank confident',
    },
    ],
};

sub _add_to_report_xml {
    my $self = shift;

    my $population_composition = Genome::Utility::MetagenomicClassifier::PopulationComposition->new(
        confidence_threshold => $self->confidence_threshold,
    );

    my $amplicons = $self->build->get_amplicons;
    unless ( $amplicons ) {
        $self->error_message( sprintf("No amplicons for build (ID %s)", $self->build_id) );
        return;
    }
    for my $amplicon ( @$amplicons ) {
        my $classification = $amplicon->get_classification
            or next;
        $population_composition->add_classification($classification);
    }

    my @domains = grep { $_ ne 'anamalia' } Genome::Utility::MetagenomicClassifier->domains;
    my @ranks = grep { $_ ne 'kingdom' } Genome::Utility::MetagenomicClassifier->taxonomic_ranks;
    pop @ranks; # remove species
    for my $domain ( @domains ) {
        my @headers = (qw/ taxonomy rank total /);
        my @ranks_involved;
        for my $rank ( @ranks ) {
            push @ranks_involved, $rank;
            my %counts = $population_composition->get_counts_for_domain_down_to_rank($domain, $rank)
                or next;
            #print Dumper(\%counts);
            my @rows;
            for my $tax ( keys %counts ) {
                push @rows, [
                $tax, 
                $rank,
                $counts{$tax}->{total}, 
                ( map { $counts{$tax}->{$_} } @ranks_involved),
                ];
            }
            $self->_add_dataset(
                #print Dumper({
                name => join('-', lc($domain), lc($rank), 'counts'),
                headers => [ @headers, @ranks_involved ],
                row_name => 'count',
                rows => \@rows,
                #});
            ) or return;
        }
    }

    return 1;
}

1;

#$HeadURL$
#$Id$
