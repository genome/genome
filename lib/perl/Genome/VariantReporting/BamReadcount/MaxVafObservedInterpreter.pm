package Genome::VariantReporting::BamReadcount::MaxVafObservedInterpreter;

use strict;
use warnings;
use Genome;
use List::Util qw/ max /;
use List::MoreUtils qw( each_array );

class Genome::VariantReporting::BamReadcount::MaxVafObservedInterpreter {
    is => ['Genome::VariantReporting::Framework::Component::Interpreter'],
    has => [
        tumor_sample_names => {
            is => 'Text',
            is_many => 1,
            is_translated => 1,
        },
        normal_sample_names => {
            is => 'Text',
            is_many => 1,
            is_translated => 1,
        },
    ],
};

sub name {
    return 'max-vaf-observed';
}

sub requires_annotations {
    return ('bam-readcount');
}


sub available_fields {
    return qw/
        max_normal_vaf_observed
        max_tumor_vaf_observed
    /;
}

sub process_interpret_entry {
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
            my $interpreter = Genome::VariantReporting::BamReadcount::VafInterpreter->create(sample_name => $sample_name);
            my %result = $interpreter->interpret_entry($entry, $passed_alt_alleles);
            for my $alt_allele (@$passed_alt_alleles) {
                $vaf_hash_ref->{$alt_allele}->{$sample_name} = $result{$alt_allele}->{vaf};
            }
        }
    }

    my %return_values;
    for my $alt_allele (@$passed_alt_alleles) {
        my @normal_vafs = values %{$normal_vafs{$alt_allele}};
        my @tumor_vafs = values %{$tumor_vafs{$alt_allele}};
        $return_values{$alt_allele} = {
            max_normal_vaf_observed => max(@normal_vafs),
            max_tumor_vaf_observed => max(@tumor_vafs),
        };
    }
    return %return_values;
}

1;
