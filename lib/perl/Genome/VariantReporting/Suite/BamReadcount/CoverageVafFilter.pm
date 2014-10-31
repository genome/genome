package Genome::VariantReporting::Suite::BamReadcount::CoverageVafFilter;

use strict;
use warnings;

use Genome;

use Genome::VariantReporting::Suite::BamReadcount::VafCalculator;

class Genome::VariantReporting::Suite::BamReadcount::CoverageVafFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::Suite::BamReadcount::ComponentBase'],
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

    while (my ($coverage, $value) = each %{$self->coverages_and_vafs}) {
        my $error = $self->_validate_is_int(
            name => 'coverages_and_vafs',
            display_name => 'coverage',
            value => $coverage,
        );
        push @errors, $error and next if $error;
        $error = $self->_validate_is_int(
            name => 'coverages_and_vafs',
            display_name => 'vaf',
            value => $value,
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

sub requires_annotations {
    return (qw/ bam-readcount /);
}

sub coverages {
    return sort { $b <=> $a } keys %{$_[0]->coverages_and_vafs};
}

sub filter_entry {
    my ($self, $entry) = @_;

    my %return_values = map { $_ => 0 } @{$entry->{alternate_alleles}};
    my @sample_alt_alleles = $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    unless (@sample_alt_alleles) {
        return $self->pass_all_sample_alts($entry);
    }

    #Keep positions without readcount information
    my $readcount_entry = $self->get_readcount_entry($entry);
    unless (defined($readcount_entry)) {
        return $self->pass_all_sample_alts($entry);
    }

    my %vafs = Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_vaf_for_all_alts(
        $entry,
        $readcount_entry,
    );
    my $vaf_for_coverage = $self->_vaf_for_coverage( $readcount_entry->depth );
    unless (defined $vaf_for_coverage) {
        return $self->pass_all_sample_alts($entry);
    }

    for my $alt_allele ( @sample_alt_alleles ) {
        #Keep positions with readcount and coverage of 0
        if ($vafs{$alt_allele} == 0) {
            $return_values{$alt_allele} = 1;
        }
        elsif ( $vafs{$alt_allele} >= $vaf_for_coverage ) {
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

sub vcf_id {
    my $self = shift;
    return 'COVERAGE_VAF';
}

sub vcf_description {
    my $self = shift;

    my @criteria = map {
        my $coverage = $_;
        my $vaf = $self->_vaf_for_coverage($coverage);
        $coverage . "X and VAF>$vaf%";
    } $self->coverages;

    return sprintf(
        'Coverage and VAF for sample %s follows criteria: [%s]',
        $self->sample_name,
        join(',', @criteria)
    );
}

1;

