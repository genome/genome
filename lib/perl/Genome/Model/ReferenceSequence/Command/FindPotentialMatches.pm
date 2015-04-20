package Genome::Model::ReferenceSequence::Command::FindPotentialMatches;

use strict;
use warnings;

use Genome;
use Lingua::EN::Inflect;
use UR::Util;

class Genome::Model::ReferenceSequence::Command::FindPotentialMatches {
    is => 'Command::V2',
    has => [
        query_file => {
            is => 'FilePath',
            doc => 'any TSV where the first column contains chromosome names',
        },
    ],
    has_optional => [
        taxon => {
            is => 'Genome::Taxon',
            doc => 'limit search to references matching this taxon',
        },
        include_subsets => {
            is => 'Boolean',
            doc => 'include references that contain a subset of the chromosome names in the query',
            default_value => 0,
        },
        include_supersets => {
            is => 'Boolean',
            doc => 'include references that contain a superset of the chromosome names in the query',
            default_value => 1,
        },
        include_partial_matches => {
            is => 'Boolean',
            doc => 'include references which contain any of the chromosome names in the query',
            default_value => 0,
        },
    ],
    has_transient_optional_output => [
        exact_matches => {
            is => 'ARRAY',
            doc => 'Exact matches that were found',
        },
        supersets => {
            is => 'ARRAY',
            doc => 'Supersets that were found (or undef if --include-supersets was not used).',
        },
        subsets => {
            is => 'ARRAY',
            doc => 'Subsets that were found (or undef if --include-subsets was not used).',
        },
        partial_matches => {
            is => 'ARRAY',
            doc => 'Parial matches that were found (or undef if --include-partial-matches was not used).',
        },
    ],
    doc => 'find reference sequences whose chromosome names match a file',
};

sub help_detail {
    return <<EOHELP;
This command searches through all reference sequences to find ones with chromosome names that match those present in a file.  This does not guarantee that the file actually comes from the references that match--only that the chromosome names are the same.  (For example, 1-22,X,Y will find multiple reference versions of the human reference which are supersets, but it cannot distinguish between them.)
EOHELP
}

sub execute {
    my $self = shift;

    my $query_chromosomes = $self->_query_chromosomes;

    my $reference_iterator = $self->_reference_iterator;

    my @matches;
    my @supersets;
    my @subsets;
    my @partials;
    while (my $reference = $reference_iterator->next) {
        my $reference_chromosomes = $reference->chromosome_array_ref(create_seqdict => 0);
        next unless defined $reference_chromosomes;

        my ($matching, $only_query, $only_reference) = UR::Util::intersect_lists($query_chromosomes, $reference_chromosomes);
        if(@$only_query == 0 and @$only_reference == 0) {
            push @matches, $reference;
        } elsif ($self->include_supersets and @$only_query == 0) {
            push @supersets, [$reference, $only_reference];
        } elsif ($self->include_subsets and @$only_reference == 0) {
            push @subsets, [$reference, $only_query];
        } elsif ($self->include_partial_matches and @$matching > 0) {
            push @partials, [$reference, $only_query, $only_reference];
        }
    }

    $self->exact_matches(\@matches);
    $self->status_message(
        'Found %s that exactly match.',
        Lingua::EN::Inflect::NO('reference', scalar(@matches)),
    );
    for my $match (@matches) {
        $self->status_message('  %s', $match->__display_name__);
    }

    if ($self->include_supersets) {
        $self->supersets(\@supersets);
        $self->status_message(
            'Found %s that are supersets.',
            Lingua::EN::Inflect::NO('reference', scalar(@supersets)),
        );
        for my $superset (@supersets) {
            $self->status_message(
                '  %s (extra: %s)',
                $superset->[0]->__display_name__,
                join(' ', @{$superset->[1]}),
            );
        }
    }

    if ($self->include_subsets) {
        $self->subsets(\@subsets);
        $self->status_message(
            'Found %s that are subsets.',
            Lingua::EN::Inflect::NO('reference', scalar(@subsets)),
        );
        for my $subset (@subsets) {
            $self->status_message(
                '  %s (missing: %s)',
                $subset->[0]->__display_name__,
                join(' ', @{$subset->[1]}),
            );
        }
    }

    if ($self->include_partial_matches) {
        $self->partial_matches(\@partials);
        $self->status_message(
            'Found %s that are partial matches.',
            Lingua::EN::Inflect::NO('reference', scalar(@partials)),
        );
        for my $partial (@partials) {
            $self->status_message(
                '  %s (missing: %s) (extra: %s)',
                $partial->[0]->__display_name__,
                join(' ', @{$partial->[1]}),
                join(' ', @{$partial->[2]}),
            );
        }
    }

    return 1;
}

sub _query_chromosomes {
    my $self = shift;

    my $query_file = $self->query_file;
    my $query_fh = Genome::Sys->open_file_for_reading($query_file);

    my %chromosomes_found;
    while (my $line = <$query_fh>) {
        chomp $line;
        my ($chr) = split("\t", $line);
        $chromosomes_found{$chr} = 1;
    }

    my $query_chromosomes = [keys %chromosomes_found];
    return $query_chromosomes;
}

sub _reference_iterator {
    my $self = shift;

    my @params = (
        'is_rederivable !=' => 1,
        status => 'Succeeded',
        subclass_name => ['Genome::Model::Build::ReferenceSequence', 'Genome::Model::Build::ImportedReferenceSequence'],
    );
    if(my $taxon = $self->taxon) {
        push @params, subject => $taxon;
    }

    my $iterator = Genome::Model::Build::ReferenceSequence->create_iterator(@params);
    return $iterator;
}

1;
