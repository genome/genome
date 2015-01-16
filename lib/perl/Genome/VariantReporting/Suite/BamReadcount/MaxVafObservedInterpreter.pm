package Genome::VariantReporting::Suite::BamReadcount::MaxVafObservedInterpreter;

use strict;
use warnings;
use Genome;
use List::Util qw/ max /;
use List::MoreUtils qw( each_array );
use Scalar::Util qw( looks_like_number );

class Genome::VariantReporting::Suite::BamReadcount::MaxVafObservedInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter'],
    has_optional => [
        tumor_sample_names => {
            is => 'Text',
            is_many => 1,
            is_translated => 1,
            doc => 'A list of tumor samples to use for calculation'
        },
        normal_sample_names => {
            is => 'Text',
            is_many => 1,
            is_translated => 1,
            doc => 'A list of normal samples to use for calculation'
        },
    ],
    doc => 'Calculate the maximum vaf value observed between all the specified normal samples, and the maximum vaf value observed between all the tumor samples',
};

sub name {
    return 'max-vaf-observed';
}

sub requires_annotations {
    return ('bam-readcount');
}

sub field_descriptions {
    my $self = shift;
    return (
        max_normal_vaf_observed => sprintf('The maximum vaf value observed between all the normal samples : %s', join(', ', $self->normal_sample_names)),
        max_tumor_vaf_observed => sprintf('The maximum vaf value observed between all the tumor samples : %s', join(', ', $self->tumor_sample_names)),
    );
}

sub _interpret_entry {
    my $self = shift;
    my $entry = shift;
    my $passed_alt_alleles = shift;

    my %normal_vafs;
    my %tumor_vafs;
    my @sample_name_accessors = qw/normal_sample_names tumor_sample_names/;
    my @vaf_hash_names        = (\%normal_vafs,        \%tumor_vafs);
    my $it                    = each_array(@sample_name_accessors, @vaf_hash_names);
    while ( my ($sample_name_accessor, $vaf_hash_ref) = $it->() ) {
        for my $sample_name ($self->$sample_name_accessor) {
            my $interpreter = Genome::VariantReporting::Suite::BamReadcount::VafInterpreter->create(
                sample_name => $sample_name,
            );
            my %result = $interpreter->interpret_entry($entry, $passed_alt_alleles);
            for my $alt_allele (@$passed_alt_alleles) {
                my $vaf = $result{$alt_allele}->{$interpreter->create_sample_specific_field_name('vaf')};
                $vaf_hash_ref->{$alt_allele}->{$sample_name} = $vaf;
            }
        }
    }

    my %return_values;
    for my $alt_allele (@$passed_alt_alleles) {
        my @normal_vafs = grep { looks_like_number($_) } values %{$normal_vafs{$alt_allele}};
        my @tumor_vafs = grep { looks_like_number($_) } values %{$tumor_vafs{$alt_allele}};

        if (@normal_vafs) {
            $return_values{$alt_allele}{max_normal_vaf_observed} = max(@normal_vafs);
        } else {
            $return_values{$alt_allele}{max_normal_vaf_observed} = $self->interpretation_null_character;
        }
        if (@tumor_vafs) {
            $return_values{$alt_allele}{max_tumor_vaf_observed} = max(@tumor_vafs);
        } else {
            $return_values{$alt_allele}{max_tumor_vaf_observed} = $self->interpretation_null_character;
        }
    }
    return %return_values;
}

1;
