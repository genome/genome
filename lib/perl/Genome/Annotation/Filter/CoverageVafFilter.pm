package Genome::Annotation::Filter::CoverageVafFilter;

use strict;
use warnings;

use Genome;
use Genome::Annotation::Expert::BamReadcount::VafCalculator;

class Genome::Annotation::Filter::CoverageVafFilter {
    is => ['Genome::Annotation::Filter::Base', 'Genome::Annotation::Expert::BamReadcount::ComponentBase'],
    has => {
        coverages_and_vafs => {
            is => 'Text',
            is_many => 1,
        },
    },
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    my %coverages_and_vafs;
    for my $coverage_and_vaf ( $self->coverages_and_vafs ) {
        my ($coverage, $vaf) = split(/:/, $coverage_and_vaf);
        my $error = $self->_validate_is_int(
            name => 'coverages_and_vafs',
            display_name => 'coverage',
            value => $coverage,
        );
        push @errors, $error and next if $error;
        $error = $self->_validate_is_int(
            name => 'coverages_and_vafs',
            display_name => 'vaf',
            value => $vaf,
            max => 100,
        );
        push @errors, $error and next if $error;
        if ( exists $coverages_and_vafs{$coverage} ) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [qw/ coverages_and_vafs /],
                desc => "Vaf for coverage ($coverage) already set!",
            );
        }
        $coverages_and_vafs{$coverage} = $vaf;
    }
    $self->{_coverages_and_vafs} = \%coverages_and_vafs;

    my @coverages = sort { $b <=> $a } keys %coverages_and_vafs;
    if ( $coverages[$#coverages] != 1 ) { #FIXME needed?? default?
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [qw/ coverages_and_vafs /],
            desc => "There is no VAF for coverage of 1!",
        );
    }
    $self->{_coverages} = \@coverages;

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

sub filter_entry {
    my ($self, $entry) = @_;

    my $readcount_entry = $self->get_readcount_entry($entry);
    return map { $_ => 0 } @{$entry->{alternate_alleles}} if not $readcount_entry;

    my @sample_alt_alleles = sort $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    my %vafs = Genome::Annotation::Expert::BamReadcount::VafCalculator::calculate_vaf_for_multiple_alleles(
        $self->get_readcount_entry($entry),
        \@sample_alt_alleles,
    );

    my $coverage = $readcount_entry->depth;
    my $vaf_for_coverage = $self->_vaf_for_coverage($coverage);

    my %return_values = map { $_ => 0 } @{$entry->{alternate_alleles}};
    return %return_values if not $vaf_for_coverage;

    for my $alt_allele ( @sample_alt_alleles ) {
        #if    (cov > 20) { if (vaf < 5 ) {fail} else {pass}
        #elsif (cov > 10) { if (vaf < 10) {fail} else {pass}
        if ( $vafs{$alt_allele} >= $vaf_for_coverage ) {
            $return_values{$alt_allele} = 1;
        }
    }

    return %return_values;
}

sub _vaf_for_coverage {
    my ($self, $coverage) = @_;

    for my $coverage_for_vaf ( @{$self->{_coverages}} ) {
        return $self->{_coverages_and_vafs}->{$coverage_for_vaf} if $coverage >= $coverage_for_vaf;
    }

    return 0; # FIXME 
}

1;

