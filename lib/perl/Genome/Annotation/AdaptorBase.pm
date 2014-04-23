package Genome::Annotation::AdaptorBase;

use strict;
use warnings;
use Genome;
use Params::Validate qw(validate validate_pos :types);

class Genome::Annotation::AdaptorBase {
    is => ['Command::V2', 'Genome::Annotation::ComponentBase'],
    is_abstract => 1,
    attributes_have => {
        is_planned => {
            is => "Boolean",
            default => 0,
        },
    },
    has_input => [
        build_id => {
            is => 'Text',
        },
        variant_type => {
            is => 'Text',
            is_output => 1,
            valid_values => ['snvs', 'indels'],
        }
    ],
    has_output => [
        bam_results => {
            is_many => 1,
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
        output_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Vcf',
        },
    ],
};

sub shortcut {
    #TODO
}

sub build {
    my $self = shift;
    return Genome::Model::Build->get($self->build_id);
}

sub execute {
    my $self = shift;
    $self->resolve_bam_results;
    $self->resolve_output_result;
    $self->resolve_plan_attributes;
    return 1;
}

sub resolve_bam_results {
    my $self = shift;

    my $results;
    if ($self->build->isa('Genome::Model::Build::SomaticVariation')) {
        $results = $self->_resolve_bam_results_variation;
    } elsif ($self->build->isa('Genome::Model::Build::SomaticValidation')) {
        $results = $self->_resolve_bam_results_validation;
    } else {
        die "This adaptor can only work on SomaticValidation or SomaticVariation type builds";
    }
    $self->bam_results($results);
}

sub _resolve_bam_results_variation {
    my $self = shift;
    my @bam_results;
    for my $type qw(normal_build tumor_build) {
        push @bam_results, $self->build->$type->merged_alignment_result;
    }
    return \@bam_results;
}

sub _resolve_bam_results_validation {
    my $self = shift;
    return [$self->build->control_merged_alignment_result, $self->build->merged_alignment_result];
}

sub resolve_output_result {
    my $self = shift;

    my $accessor = sprintf('get_detailed_%s_vcf_result', $self->variant_type);
    my $result = eval {$self->build->$accessor};
    my $error = $@;
    if ($error) {
        $self->debug_message("No %s result found on build %s:\n%s",
            $self->variant_type, $self->build->id, $error);
    }
    $self->output_result($result);
}

sub resolve_plan_attributes {
    my $self = shift;

    my $annotation_plan = $self->build->annotation_plan;
    my $specific_plan = $annotation_plan->get_plan($self->category, $self->name);
    for my $name (keys %{$specific_plan->params}) {
        $self->$name($specific_plan->params->{$name});
    }
}

sub planned_output_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(
        is_output => 1, is_planned => 1);
    return map {$_->property_name} @properties;
}

sub validate_with_plan_params {
    my ($self, $params) = validate_pos(@_, 1, 1);

    my $needed = Set::Scalar->new($self->planned_output_names);
    my $have = Set::Scalar->new(keys %{$params});

    unless($needed->is_equal($have)) {
        if (my $still_needed = $needed - $have) {
            $self->error_message("Parameters required by adaptor (%s) but not provided: (%s)",
                $self->class, join(",", $still_needed->members));
        }
        if (my $not_needed = $have - $needed) {
            $self->error_message("Parameters provided by plan but not required by adaptor (%s): (%s)",
                $self->class, join(",", $not_needed->members));
        }
        die $self->error_message("Provided parameters and required parameters do not match");
    }

}

1;
