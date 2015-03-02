package Genome::VariantReporting::Suite::BamReadcount::Run;

use strict;
use warnings FATAL => 'all';
use Genome;
use version 0.77;
use Params::Validate qw(validate_pos);

class Genome::VariantReporting::Suite::BamReadcount::Run {
    is => 'Genome::VariantReporting::Framework::Component::Expert::Command',
    has_planned_input => [
        aligned_bam_result_id => {
            is => 'Text',
            is_translated => 1,
            doc => 'The bam result used to calculate read counts',
        },
        version => {
            is  => 'Version',
            example_values => ['0.5'],
            doc => "bam-readcount version to be used.",
        },
    ],
    has_planned_optional => [
        minimum_mapping_quality => {
            is => 'Integer',
            example_values => [0],
            doc => "filter reads with mapping quality less than this. This is the -q parameter.",
        },
        minimum_base_quality => {
            is => 'Integer',
            example_values => [0],
            doc => "don't include reads where the base quality is less than this. This is the -b parameter. This is only available in versions 0.3 and later.",
        },
        max_count => {
            is  => 'Integer',
            example_values => [10_000_000],
            doc => "max depth to avoid excessive memory. This is the -d parameter in version 0.5.",
        },
        per_library => {
            is  => 'Boolean',
            example_values => [0],
            doc => "report results per library. This is the -p parameter in version 0.5.",
        },
        insertion_centric => {
            is  => 'Boolean',
            example_values => [0],
            doc => "do not include reads containing insertions after the current position in per-base counts. This is the -i parameter in version 0.5.",
        },
    ],
    doc => 'Add bam readcount information to a vcf',
};

sub name {
    'bam-readcount';
}

sub result_class {
    'Genome::VariantReporting::Suite::BamReadcount::RunResult';
}

sub __errors__ {
    my $self = shift;
    my @errors = $self->SUPER::__errors__;

    unless ($self->version >= 0.5) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => ['version'],
            desc => sprintf("Version provided (%s) must be greater than or equal to 0.5", $self->version),
        );
    }
    return @errors;
}

my $MIN_VERSION = '0.7';

sub __planned_errors__ {
    my ($self, $params) = validate_pos(@_, 1, 1);

    my @errors = $self->SUPER::__planned_errors__($params);

    if (version->parse("v".$params->{version}) < version->parse("v$MIN_VERSION")) {
        push @errors, UR::Object::Tag->create(
            type => 'error',
            properties => ['version'],
            desc => sprintf("The BamReadcount expert requires version (%s) ".
                    "or higher, but you supplied version (%s).",
                    $MIN_VERSION, $params->{version} || 'undef'),
        );
    }
    return @errors;
}

sub resolve_plan_attributes {
    my $self = shift;

    my $aligned_bam_result_id = $self->aligned_bam_result_id;
    $self->SUPER::resolve_plan_attributes;
    $self->aligned_bam_result_id($aligned_bam_result_id);
}

sub all_translated_inputs {
    my $class = shift;
    return grep {$_->property_name ne 'aligned_bam_result_id'} $class->SUPER::all_translated_inputs;
}

sub all_translated_is_many_inputs {
    my $class = shift;
    my @inputs = $class->SUPER::all_translated_is_many_inputs;
    push @inputs, $class->__meta__->properties(property_name => 'aligned_bam_result_id');
    return @inputs;
}


1;
