package Genome::VariantReporting::Framework::ReportResult;

use strict;
use warnings FATAL => 'all';
use Genome::File::Vcf::Reader;
use Genome;

class Genome::VariantReporting::Framework::ReportResult {
    is => 'Genome::SoftwareResult::Stageable',
    has_input => [
        input_vcf_lookup => {
            is => 'Text',
        },
    ],
    has_param => [
        plan_json_lookup => {
            is => 'Text',
        },
        variant_type => {
            is => 'Text',
            valid_values => ['snvs', 'indels'],
        },
    ],
    has_transient_optional => [
        input_vcf => {
            is => 'Path',
        },
        plan_json => {
            is => 'Text',
        },
        plan => {
            is => 'Genome::VariantReporting::Framework::Plan::MasterPlan',
        },
    ],
};

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;

    $self->_prepare_staging_directory;
    $self->_run;

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    return File::Spec->join('/', 'model_data', 'software-result', $self->id);
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

sub plan {
    my $self = shift;

    unless (defined($self->__plan)) {
        my $plan = Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_json($self->plan_json);
        $self->__plan($plan);
    }
    return $self->__plan;
}

sub _run {
    my $self = shift;

    $self->status_message("Reading from: ".$self->input_vcf."\n");

    my @reporters = $self->create_reporters;
    $self->initialize_reporters(@reporters);

    my $vcf_reader = Genome::File::Vcf::Reader->new($self->input_vcf);
    while (my $entry = $vcf_reader->next) {
        for my $reporter (@reporters) {
            $reporter->process_entry($entry);
        }
    }

    $self->finalize_reporters(@reporters);
    return 1;
}

sub create_reporters {
    my $self = shift;

    my @reporters;
    for my $reporter_plan ($self->plan->reporter_plans) {
        push @reporters, $reporter_plan->object;
    }

    return @reporters;
}

sub initialize_reporters {
    my $self = shift;
    map {$_->initialize($self->temp_staging_directory)} @_;
    return;
}

sub finalize_reporters {
    my $self = shift;
    map {$_->finalize()} @_;
    return;
}


1;
