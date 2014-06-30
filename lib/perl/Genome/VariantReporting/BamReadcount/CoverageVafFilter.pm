package Genome::VariantReporting::BamReadcount::CoverageVafFilter;

use strict;
use warnings;

use Genome;

use Genome::VariantReporting::BamReadcount::VafCalculator;
require Memoize;

class Genome::VariantReporting::BamReadcount::CoverageVafFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::BamReadcount::ComponentBase'],
    has => {
        coverages_and_vafs => {
            is => 'HASH',
            doc => 'A mapping of coverages and the minimum vafs.'
        },
    },
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    for my $coverage ( keys %{$self->coverages_and_vafs} ) {
        my $error = $self->_validate_is_int(
            name => 'coverages_and_vafs',
            display_name => 'coverage',
            value => $coverage,
        );
        push @errors, $error and next if $error;
        $error = $self->_validate_is_int(
            name => 'coverages_and_vafs',
            display_name => 'vaf',
            value => $self->coverages_and_vafs->{$coverage},
            max => 100,
        );
        push @errors, $error if $error;
    }

    return @errors;
}

sub _validate_is_int {
    my ($self, %property) = @_;

    Carp::confess('No property name given to check if value is whole percent!') if not $property{name};
    my $display_name = $property{display_name} // $property{name};

    if ( not defined $property{value} ) {
        return UR::Object::Tag->create(
            type => 'error',
            properties => [ $property{name} ],
            desc => "No value given!",
        );
    }

    if ( $property{value} !~ /^\d+$/ ) {
        return UR::Object::Tag->create(
            type => 'error',
            properties => [ $property{name} ],
            desc => "Value given ($property{value}) for $display_name is not a whole number between 0 and 100!",
        );
    }

    if ( defined $property{max} and $property{value} > $property{max} ) {
        return UR::Object::Tag->create(
            type => 'error',
            properties => [ $property{name} ],
            desc => "Value given ($property{value}) for $display_name is not an integer below $property{max}!",
        );
    }

    return;
}

sub name {
    return 'coverage-vaf';
}

sub requires_experts {
    return (qw/ bam-readcount /);
}

sub coverages {
    return sort { $b <=> $a } keys %{$_[0]->coverages_and_vafs};
}
Memoize::memoize('coverages');

sub filter_entry {
    my ($self, $entry) = @_;

    my $readcount_entry = $self->get_readcount_entry($entry);
    return map { $_ => 0 } @{$entry->{alternate_alleles}} if not $readcount_entry;

    my @sample_alt_alleles = sort $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    my %vafs = Genome::VariantReporting::BamReadcount::VafCalculator::calculate_vaf_for_all_alts(
        $entry,
        $self->get_readcount_entry($entry),
    );

    my $vaf_for_coverage = $self->_vaf_for_coverage( $readcount_entry->depth );

    my %return_values = map { $_ => 0 } @{$entry->{alternate_alleles}};
    return %return_values if not defined $vaf_for_coverage;

    for my $alt_allele ( @sample_alt_alleles ) {
        if ( $vafs{$alt_allele} >= $vaf_for_coverage ) {
            $return_values{$alt_allele} = 1;
        }
    }

    return %return_values;
}

sub _vaf_for_coverage {
    my ($self, $coverage) = @_;

    for my $coverage_for_vaf ( $self->coverages ) {
        return $self->coverages_and_vafs->{$coverage_for_vaf} if $coverage >= $coverage_for_vaf;
    }

    return;
}

1;

