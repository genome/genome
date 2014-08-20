package Genome::Model::Tools::DetectVariants2::Result::Filter;

use strict;
use warnings;

use Genome;

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

#Most filter-specific logic is in Detector.pm

sub previous_result {
    my $self = shift;

    my @users = Genome::SoftwareResult::User->get(user=>$self);
    if (scalar(@users) == 1) {
        my $user = shift @users;
        return $user->software_result;
    } else {
        my $message = sprintf("Number of previous results (%d) is not 1.  Result IDs: (%s)",
            scalar(@users), join(', ', map {$_->software_result->id} @users));
        die $message;
    }
}

sub vcf_result_params {
    my $self = shift;

    return (
        filter_name => $self->filter_name,
        filter_params => $self->filter_params,
        filter_version => $self->filter_version,
        incoming_vcf_result => $self->previous_result->get_vcf_result,
        input_id => $self->id,
        previous_filter_strategy => $self->previous_filter_strategy,
        test_name => $self->test_name,
        vcf_version => Genome::Model::Tools::Vcf->get_vcf_version,
    );
}

sub vcf_result_class {
    'Genome::Model::Tools::DetectVariants2::Result::Vcf::Filter';
}

1;
