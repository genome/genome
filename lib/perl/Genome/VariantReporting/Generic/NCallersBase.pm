package Genome::VariantReporting::Generic::NCallersBase;

use strict;
use warnings;
use Genome;
use List::MoreUtils qw/uniq/;

class Genome::VariantReporting::Generic::NCallersBase {
    is => 'Genome::VariantReporting::Framework::Component::WithSampleName',
    has => [
        valid_callers => {
            is => 'String',
            is_many => 1,
            default_value => [qw(VarscanSomatic Sniper Strelka)],
            doc => 'List of variant callers to include in determination for filtering',
        },
    ],
};

sub get_callers {
    my ($self, $entry, $alts) = @_;
    my %callers;
    for my $alt (@$alts) {
        $callers{$alt} = [];
    }
    for my $caller_name (sort $self->valid_callers) {
        my $sample_name = $self->sample_name_with_suffix($caller_name);
        my $sample_index = eval{ $entry->{header}->index_for_sample_name($sample_name) };
        my $error = $@;
        if ($error =~ /^\QSample name $sample_name not found in header\E/) {
            next;
        }
        my $caller_filtered = $entry->sample_field($sample_index, "FT");
        if (defined $caller_filtered and $caller_filtered ne "." and $caller_filtered ne "PASS") {
            next;
        }
        my @sample_alt_alleles = $entry->alt_bases_for_sample($sample_index);
        for my $sample_alt_allele (uniq @sample_alt_alleles) {
            push(@{$callers{$sample_alt_allele}}, $caller_name);
        }
    }
    return %callers;
}
1;

