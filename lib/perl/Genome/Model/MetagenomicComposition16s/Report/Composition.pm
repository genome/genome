package Genome::Model::MetagenomicComposition16s::Report::Composition;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require Genome::Utility::MetagenomicClassifier::PopulationComposition;
require Genome::Utility::MetagenomicClassifier::SequenceClassification;

class Genome::Model::MetagenomicComposition16s::Report::Composition {
    is => 'Genome::Model::MetagenomicComposition16s::Report',
    has => [
    confidence_threshold => {
        is => 'Number',
        default_value => 0.8,
        doc => 'The threshold to match or exceed (>=) to consider a taxonic rank confident',
    },
    ],
};

#< Generator >#
sub description {
    my $self = shift;

    return sprintf(
        'Metagenomic Composition Report for Model\'s (<%s> <%s>) Build <%s>)',
        $self->model_name,
        $self->model_id,
        $self->build_id,
    );
}

sub _add_to_report_xml {
    my $self = shift;

    my $population_composition = Genome::Utility::MetagenomicClassifier::PopulationComposition->new(
        confidence_threshold => $self->confidence_threshold,
    );

    my @amplicon_sets = $self->build->amplicon_sets;
    Carp::confess('No amplicon sets for '.$self->build->description) if not @amplicon_sets;

    my @domains = grep { $_ ne 'anamalia' } Genome::Utility::MetagenomicClassifier->domains;
    my @ranks = grep { $_ ne 'kingdom' } Genome::Utility::MetagenomicClassifier->taxonomic_ranks;
    pop @ranks; # remove species

    for my $amplicon_set ( @amplicon_sets ) {
        next if not $amplicon_set->amplicon_iterator; # ok
        while ( my $amplicon = $amplicon_set->next_amplicon ) {
            next if not $amplicon->{classification}; # ok
            my $classification = Genome::Utility::MetagenomicClassifier::SequenceClassification->new_from_classification_array(
                name => $amplicon->{name},
                classifier => $self->model->classifier,
                ranks => [ 'root', @ranks ],
                classifications => [ map { $amplicon->{classification}->[$_] } (2..8) ],
            );
            $population_composition->add_classification($classification);
        }

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
    }

    return 1;
}

1;

#$HeadURL: svn+ssh://svn/srv/svn/gscpan/perl_modules/trunk/Genome/Model/AmpliconAssembly/Report/Composition.pm $
#$Id: Composition.pm 53403 2009-11-30 19:08:28Z ebelter $
