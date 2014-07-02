package Genome::VariantReporting::Framework::Component::Adaptor;

use strict;
use warnings;
use Genome;
use Params::Validate qw(validate validate_pos :types);

class Genome::VariantReporting::Framework::Component::Adaptor {
    is => ['Command::V2', 'Genome::VariantReporting::Framework::Component::Base'],
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
        },
        plan_json => {
            is => 'Text',
        },
    ],
    has_output => [
        bam_results => {
            is_many => 1,
            is => 'Genome::InstrumentData::AlignedBamResult',
        },
    ],
};

sub name {
    die "Abstract";
}

sub resolve_expert_specific_attributes_from_build {
    my $self = shift;
    # This may be defined in subclasses
    return;
}

sub shortcut {
    my $self = shift;
    return $self->execute();
}

sub execute {
    my $self = shift;
    $self->debug_message("Resolving bam results");
    $self->resolve_bam_results;

    $self->debug_message("Resolving plan attributes");
    $self->resolve_plan_attributes;

    $self->debug_message("Resolving any expert-specific attributes from the build");
    $self->resolve_expert_specific_attributes_from_build;
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

sub build {
    my $self = shift;

    my $build = Genome::Model::Build->get($self->build_id);
    if ($build) {
        return $build;
    } else {
        die $self->error_message("Couldn't find a build for id (%s)",
            $self->build_id);
    }
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

sub resolve_plan_attributes {
    my $self = shift;

    my $variant_reporting_plan = $self->plan;
    my $specific_plan = $variant_reporting_plan->get_plan('expert', $self->name);
    for my $name (keys %{$specific_plan->params}) {
        $self->$name($specific_plan->params->{$name});
    }
}

sub plan {
    my $self = shift;

    return Genome::VariantReporting::Framework::Plan::MasterPlan->create_from_json($self->plan_json);
}

sub planned_output_names {
    my $self = shift;

    my @properties = $self->__meta__->properties(
        is_output => 1, is_planned => 1);
    return map {$_->property_name} @properties;
}

# TODO this is not covered by tests
sub validate_with_plan_params {
    my ($self, $params) = validate_pos(@_, 1, 1);

    my @errors = $self->__planned_output_errors__($params);
    if (@errors) {
        $self->print_errors(@errors);
        die $self->error_message("Failed to validate_with_plan_params with params:\n" . Data::Dumper::Dumper $params);
    }
    return;
}

# TODO this is not covered by tests
sub __planned_output_errors__ {
    my ($self, $params) = validate_pos(@_, 1, 1);
    my $needed = Set::Scalar->new($self->planned_output_names);
    return $self->_get_missing_errors($params, $needed),
        $self->_get_extra_errors($params, $needed);
}
sub _get_missing_errors {
    my ($self, $params, $needed) = validate_pos(@_, 1, 1, 1);

    my $have = Set::Scalar->new(keys %{$params});
    my @errors;
    unless($needed->is_equal($have)) {
        if (my $still_needed = $needed - $have) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [$still_needed->members],
                desc => sprintf("Parameters required by adaptor (%s) but not provided: (%s)", 
                    $self->class, join(",", $still_needed->members)),
            );
        }
    }

    return @errors;
}

sub _get_extra_errors {
    my ($self, $params, $needed) = validate_pos(@_, 1, 1, 1);

    my $have = Set::Scalar->new(keys %{$params});
    my @errors;
    unless($needed->is_equal($have)) {
        if (my $not_needed = $have - $needed) {
            push @errors, UR::Object::Tag->create(
                type => 'error',
                properties => [$not_needed->members],
                desc => sprintf("Parameters provided but not required by adaptor (%s): (%s)",
                    $self->class, join(",", $not_needed->members)),
            );
        }
    }

    return @errors;
}

1;
