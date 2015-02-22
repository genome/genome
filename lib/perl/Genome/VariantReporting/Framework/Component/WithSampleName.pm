package Genome::VariantReporting::Framework::Component::WithSampleName;

use strict;
use warnings;
use Genome;
use Memoize qw();

class Genome::VariantReporting::Framework::Component::WithSampleName {
    is => ['Genome::VariantReporting::Framework::Component::WithTranslatedInputs'],
    has => [
        sample_name => {
            is => 'Text',
            is_translated => 1,
        },
        sample_name_label => {
            is_optional => 1,
            is => 'Text',
            is_translated => 1,
        },
    ],
};

sub sample_index {
    my $self = shift;
    my $header = shift;

    return $header->index_for_sample_name($self->sample_name);
}

Memoize::memoize("sample_index", LIST_CACHE => 'MERGE');

sub sample_name_with_suffix {
    my $self = shift;
    my $suffix = shift;
    return sprintf("%s-[%s]", $self->sample_name, $suffix);
}

sub get_callers {
    my $self = shift;
    my $header = shift;

    my $sample_name = $self->sample_name;
    my @vcf_sample_names = @{$header->{sample_names}};
    my @callers;
    for my $vcf_sample_name (@vcf_sample_names) {
        if ($vcf_sample_name =~ /$sample_name-\[(.+)\]/) {
            push(@callers, $1);
        }
    }
    return @callers;
}

sub create_sample_specific_field_name {
    my ($self, $field) = @_;

    return join("_", $self->get_sample_label, $field);
}

sub get_sample_label {
    my $self = shift;

    if (defined($self->sample_name_label)) {
        return $self->sample_name_label;
    } else {
        return $self->sample_name;
    }
}

sub sample_names {
    my $self = shift;
    return ($self->sample_name);
}
1;

