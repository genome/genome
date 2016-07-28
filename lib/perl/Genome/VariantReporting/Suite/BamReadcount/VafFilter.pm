package Genome::VariantReporting::Suite::BamReadcount::VafFilter;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Suite::BamReadcount::VafCalculator;

class Genome::VariantReporting::Suite::BamReadcount::VafFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::Suite::BamReadcount::ComponentBase'],
    has_optional => [
        min_vaf => {
            is => 'Number',
            doc => 'The lower bound for VAF values that will pass the filter (will be kept).  Whether the bound is inclusive or exclusive is controlled by the --exclusive option.',
        },
        max_vaf => {
            is => 'Number',
            doc => 'The upper bound for VAF values that will pass the filter (will be kept).  Whether the bound is inclusive or exclusive is controlled by the --exclusive option.',
        },
        exclusive => {
            is => 'Boolean',
            doc => 'Flag to indicate whether the bounds are exclusive or inclusive.  By default they are inclusive.',
            default => 0,
        },
    ],
    doc => q{Filter variants that don't match vaf cutoffs},
};

sub name {
    return 'vaf-cutoff';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;

    unless (defined($self->min_vaf) || defined($self->max_vaf)) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [$self->min_vaf, $self->max_vaf],
            desc => "Must define at least one of min_vaf or max_vaf",
        );
    }

    if (defined($self->min_vaf) && defined($self->max_vaf)
        && $self->min_vaf > $self->max_vaf
    ) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => [$self->min_vaf, $self->max_vaf],
            desc => "Max_vaf must be larger or equal to min_vaf",
        );

    }

    return @errors;
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my %return_values = map { $_ => 0 } @{$entry->{alternate_alleles}};
    my @sample_alt_alleles = $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    unless (@sample_alt_alleles) {
        return $self->pass_all_sample_alts($entry);
    }

    #Keep positions without readcount information
    my $readcount_entries = $self->get_readcount_entries($entry);
    unless (defined($readcount_entries)) {
        return $self->pass_all_sample_alts($entry);
    }

    my %vafs = Genome::VariantReporting::Suite::BamReadcount::VafCalculator::calculate_vaf_for_all_alts(
        $entry,
        $readcount_entries
    );
    unless (%vafs) {
        return $self->pass_all_sample_alts($entry);
    }
    for my $allele (@sample_alt_alleles) {
        #Keep positions with readcount and coverage of 0
        if ($vafs{$allele} == 0) {
            $return_values{$allele} = 1;
        }
        elsif ($self->passes_filter($vafs{$allele})) {
            $return_values{$allele} = 1;
        }
    }

    return %return_values;
}

sub passes_filter {
    my $self = shift;
    my $vaf_value = shift;

    if (defined($self->min_vaf)) {
        if ($self->exclusive and $vaf_value <= $self->min_vaf or
            $self->inclusive and $vaf_value < $self->min_vaf) {
            return 0;
        }
    }

    if (defined($self->max_vaf)) {
        if ($self->exclusive and $vaf_value >= $self->max_vaf or
            $self->inclusive and $vaf_value > $self->max_vaf) {
            return 0;
        }
    }

    return 1;
}

sub inclusive {
    my $self = shift;
    return !$self->exclusive;
}


sub vcf_id {
    my $self = shift;
    return sprintf("VAF_%s_%s_%s", $self->min_vaf, $self->max_vaf, $self->sample_name);
}

sub vcf_description {
    my $self = shift;
    my $min_vaf = defined $self->min_vaf ? $self->min_vaf : 0;
    my $max_vaf = defined $self->max_vaf ? $self->max_vaf : 100;
    return sprintf("VAF value for sample %s is between %s and %s", $self->sample_name, $min_vaf, $max_vaf);
}

1;
