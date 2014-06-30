package Genome::VariantReporting::BamReadcount::GenotypeVafFilter;

use strict;
use warnings;

use Genome;
use Genome::VariantReporting::BamReadcount::VafCalculator;

class Genome::VariantReporting::BamReadcount::GenotypeVafFilter {
    is => ['Genome::VariantReporting::Framework::Component::Filter', 'Genome::VariantReporting::BamReadcount::ComponentBase'],
    has => {
        min_het_vaf => { is => 'Number', doc => 'Minimum VAF value if the genotype if heterozygous. Expressed as a whole percentage.', },
        max_het_vaf => { is => 'Number', doc => 'Maximum VAF value if the genotype if heterozygous. Expressed as a whole percentage.', },
        min_hom_vaf => { is => 'Number', doc => 'Minimum VAF value if the genotype if homozygous. Expressed as a whole percentage.', },
        max_hom_vaf => { is => 'Number', doc => 'Maximum VAF value if the genotype if homozygous. Expressed as a whole percentage.', },
    },
};

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    PARAM_BASE_NAME: for my $param_base_name (qw/ het_vaf hom_vaf /) {
        my $min_param_name = 'min_'.$param_base_name;
        my $max_param_name = 'max_'.$param_base_name;
        for my $param_name ( $min_param_name, $max_param_name ) {
            my $value = $self->$param_name;
            if ( $value !~ /^\d+$/ or $value > 100 ) {
                push @errors, UR::Object::Tag->create(
                    type => 'error',
                    properties => [ $param_name ],
                    desc => "Value given ($value) is not a whole number between 0 and 100!",
                );
            }
        }
        next PARAM_BASE_NAME if @errors;
        if ( $self->$min_param_name >= $self->$max_param_name ) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [ $min_param_name, $max_param_name ],
                desc => sprintf(
                    'Value for %s (%s) needs to be less than %s (%s)!',
                    $min_param_name, $self->$min_param_name, $max_param_name, $self->$max_param_name,
                ),
            );
        }
    }

    return @errors;
}

sub name {
    return 'genotype-vaf';
}

sub requires_experts {
    return (qw/ bam-readcount /);
}

sub filter_entry {
    my $self = shift;
    my $entry = shift;

    my @sample_alt_alleles = sort $entry->alt_bases_for_sample($self->sample_index($entry->{header}));
    my %vafs = Genome::VariantReporting::BamReadcount::VafCalculator::calculate_vaf_for_all_alts(
        $entry,
        $self->get_readcount_entry($entry),
    );
    my $gt = $entry->genotype_for_sample($self->sample_index($entry->{header}));
    my %return_values = map { $_ => 0 } @{$entry->{alternate_alleles}};

    for my $alt_allele ( @sample_alt_alleles ) {
        my $vaf = $vafs{$alt_allele};
        if ( $gt->is_missing ) {
            # fail
            $return_values{$alt_allele} = 0;
        }
        elsif ( $gt->is_homozygous ) {
            if ( $vaf < $self->min_hom_vaf or $vaf > $self->max_hom_vaf ) { 
                # fail
                $return_values{$alt_allele} = 0;
            }
            else { 
                # pass
                $return_values{$alt_allele} = 1;
            }
        }
        else { # $gt->is_heterozygous
            if ( $vaf < $self->min_het_vaf or $vaf > $self->max_het_vaf ) {
                # fail
                $return_values{$alt_allele} = 0;
            } else {
                # pass 
                $return_values{$alt_allele} = 1;
            }
        }
    }

    return %return_values;
}

1;

