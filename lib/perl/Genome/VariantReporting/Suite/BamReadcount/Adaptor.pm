package Genome::VariantReporting::Suite::BamReadcount::Adaptor;

use strict;
use warnings;
use Genome;

class Genome::VariantReporting::Suite::BamReadcount::Adaptor {
    is => 'Command::V2',
    has_input => [
        plan_json => {
            is => 'Text',
        }
    ],
    has_transient_output => [
        aligned_bam_result_id => {
            is => 'Text',
            is_many => 1,
            doc => 'The bam result used to calculate read counts',
        },
        version => {
            is  => 'Version',
            example_values => ['0.5'],
            doc => "bam-readcount version to be used.",
        },
    ],
};

sub execute {
    my $self = shift;

    $self->resolve_plan_attributes;

    return 1;
}

sub resolve_plan_attributes {
    my $self = shift;

    my $variant_reporting_plan = Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_json($self->plan_json);
    my $specific_plan = $variant_reporting_plan->get_plan('expert', 'bam-readcount');
    $self->aligned_bam_result_id($specific_plan->run_params->{aligned_bam_result_id});
    $self->version($specific_plan->run_params->{version});
    return;
}

1;
