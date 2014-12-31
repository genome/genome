package Genome::Model::Tools::DetectVariants2::Result::Filter;

use strict;
use warnings;

use Genome;
use Data::Dump qw(pp);

class Genome::Model::Tools::DetectVariants2::Result::Filter {
    is => ['Genome::Model::Tools::DetectVariants2::Result::DetectionBase'],
    has_param => [
        filter_name => {
            is => 'Text',
            doc => 'The name of the filter to use',
        },
        filter_params => {
            is => 'Text',
            is_optional => 1,
            doc => 'Additional parameters to pass to the filter',
        },
        filter_version => {
            is => 'Text',
            is_optional => 1,
            doc => 'Version of the filter to use',
        },
        previous_filter_strategy => {
            is => 'Text',
            is_optional => 1,
            doc => 'the dispatcher string corresponding to the previously run filters for the data being filtered',
        },
    ],
};

sub previous_result {
    my $self = shift;

    my @parents = grep {pp($_->test_name) eq pp($self->test_name)} $self->parents;
    if (scalar(@parents) == 1) {
        return shift @parents;
    } else {
        my $message = sprintf("Number of previous results (%d) is not 1.  Result IDs: (%s)",
            scalar(@parents), join(', ', map {$_->id} @parents));
        die $message;
    }
}

sub vcf_result_params {
    my $self = shift;
    my $aligned_reads_sample = shift;
    my $control_aligned_reads_sample = shift;

    return (
        filter_name => $self->filter_name,
        filter_params => $self->filter_params,
        filter_version => $self->filter_version,
        incoming_vcf_result => $self->previous_result->get_vcf_result($aligned_reads_sample, $control_aligned_reads_sample),
        input_id => $self->id,
        previous_filter_strategy => $self->previous_filter_strategy,
        test_name => $self->test_name,
        vcf_version => Genome::Model::Tools::Vcf->get_vcf_version,
        aligned_reads_sample => $aligned_reads_sample,
        ($control_aligned_reads_sample? (control_aligned_reads_sample => $control_aligned_reads_sample) : ()),
    );
}

sub vcf_result_class {
    'Genome::Model::Tools::DetectVariants2::Result::Vcf::Filter';
}

1;
