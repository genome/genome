package Genome::VariantReporting::Framework::Component::Filter;

use strict;
use warnings FATAL => 'all';
use Genome;
use Carp qw(confess);

class Genome::VariantReporting::Framework::Component::Filter {
    is => 'Genome::VariantReporting::Framework::Component::Interpreter',
    is_abstract => 1,
};

sub filter_entry {
    my $self = shift;
    my $class = $self->class;
    confess "Abstract method 'filter_entry' must be defined in class '$class'";
}

sub field_descriptions {
    my $self = shift;
    return (
        filter_status => $self->vcf_description
    );
}

sub interpret_entry {
    my ($self, $entry, $passed_alt_alleles) = @_;

    my %filter_output = $self->filter_entry($entry);

    my %return_values;
    for my $variant_allele (@$passed_alt_alleles) {
        $return_values{$variant_allele} =  {
            'filter_status' => $filter_output{$variant_allele},
        };
    }
    return %return_values;
}

# Auto pass all alts for this sample (fail alts not present in this sample)
# If this sample has no genotype information, pass everything.
sub pass_all_sample_alts {
    my ($self, $entry) = @_;
    my @sample_alt_alleles = $entry->alt_bases_for_sample($self->sample_index($entry->{header}));

    my %return_values;
    if (@sample_alt_alleles) {
        %return_values = map { $_ => 0 } @{$entry->{alternate_alleles}};
        for my $alt_allele (@sample_alt_alleles) {
            $return_values{$alt_allele} = 1;
        }
    } else {
        %return_values = map { $_ => 1 } @{$entry->{alternate_alleles}};
    }

    return %return_values;
}


sub vcf_description {
    my $self = shift;
    return sprintf("Override method 'vcf_description' must be defined in class '%s' with the real description", $self->class);
}

sub vcf_id {
    my $self = shift;
    my $class_name = $self->class;
    $class_name =~ s/::/_/g;
    return $class_name;
}

1;
