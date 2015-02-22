package Genome::VariantReporting::Suite::Fpkm::MinFpkmFilter;

use strict;
use warnings;
use Genome;
use Scalar::Util qw(looks_like_number);

class Genome::VariantReporting::Suite::Fpkm::MinFpkmFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::Suite::Fpkm::EntryParser'],
    has => [
        min_fpkm => {
            is => "Number",
            doc => "The minimum FPKM value to pass",
        },
    ],
    doc => q{Filter variants that don't meet minimum FPKM value},
};

sub name {
    return 'min-fpkm';
}

sub requires_annotations {
    return ('fpkm');
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__;

    my $min_fpkm = $self->min_fpkm;
    unless ( ( looks_like_number($min_fpkm) ) and ($min_fpkm >= 0) ){
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [qw/ min_fpkm /],
            desc => "Value given for min_fpkm ($min_fpkm) is not a positive number!",
        );
    }

    return @errors;
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values;
    my %fpkm_for_genotype_allele = $self->fpkm_for_genotype_allele($entry);
    for my $alt_allele ( @{$entry->{alternate_alleles}} ) {
        $return_values{$alt_allele} = $self->value_passes($fpkm_for_genotype_allele{$alt_allele});
    }

    return %return_values;
}

sub value_passes {
    my ($self, $value) = @_;
    if ( (not defined $value ) or ($value eq '.') or ($value < $self->min_fpkm) ) {
        return 0;
    } else {
        return 1;
    }
}

sub vcf_id {
    my $self = shift;
    return sprintf("MIN_FPKM_%s_%s", $self->min_fpkm, $self->sample_name);
}

sub vcf_description {
    my $self = shift;
    return sprintf("The FPKM value for sample %s is at least %s", $self->sample_name, $self->min_fpkm);
}

1;
