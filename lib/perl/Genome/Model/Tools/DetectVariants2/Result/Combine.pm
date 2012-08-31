package Genome::Model::Tools::DetectVariants2::Result::Combine;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::DetectVariants2::Result::Combine {
    is  => ['Genome::Model::Tools::DetectVariants2::Result::Base'],
    has_input => [
        input_a_id => {
            is => 'Text',
        },
        input_b_id => {
            is => 'Text',
        },
    ],
    has_param => [
        version => {
            is => 'Integer',
        },
    ],
    has_optional => [
        # these should be DV2 results but there is not yet a common base for all types of results
        # filters and detectors are Result/Base
        # combine are not (because Result/Base has a lot of properties that aren't needed)
        input_a => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            id_by => 'input_a_id',
        },
        input_b => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            id_by => 'input_b_id',
        },
        input_directory_a => {
            is => 'Text',
            via => 'input_a',
            to => 'output_dir',
        },
        input_directory_b => {
            is => 'Text',
            via => 'input_b',
            to => 'output_dir',
        },
    ],
};

sub _variant_type { die 'override _variant_type' };

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return unless ($self);

    unless($self->_validate_inputs) {
        die $self->error_message('Failed to validate inputs.');
    }

    unless($self->_prepare_staging_directory) {
        die $self->error_message('Failed to prepare staging directory.');
    }

    unless($self->_combine_variants){
        die $self->error_message('Failed to combine variants');
    }

    unless($self->_validate_output) {
        die $self->error_message('Failed to validate output.');
    }

    unless ($self->_prepare_output_directory) {
        die $self->error_message('Failed to prepare output directory.');
    }

    unless($self->_promote_data) {
        die $self->error_message('Failed to promote data.');
    }

    unless($self->_reallocate_disk_allocation) {
        die $self->error_message('Failed to reallocate disk allocation.');
    }

    unless($self->_add_as_user_of_inputs) {
        die $self->error_message('Failed to add self as user of inputs.');
    }

    return $self;
}

sub _gather_params_for_get_or_create {
    my $class = shift;

    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }
    }

    my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my %software_result_params = (
        params_id => $params_bx->id,
        inputs_id => $inputs_bx->id,
        subclass_name => $class,
    );

    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs => \%is_input,
        params => \%is_param,
    };
}

sub _needs_symlinks_followed_when_syncing { die "overload this in subclass" };
sub _working_dir_prefix { die "overload this in subclass" };
sub resolve_allocation_disk_group_name { die "overload this in subclass" };
sub allocation_subdir_prefix { die "overload this in subclass" };

sub _combine_variants {
    die "overload this function to do work";
}

sub estimated_kb_usage {
    my $self = shift;

    my $input_a = $self->input_a;
    my $input_b = $self->input_b;

    my $input_a_allocation = $input_a->disk_allocations;
    my $input_b_allocation = $input_b->disk_allocations;

    my $estimated_kb_usage;
    $estimated_kb_usage += $input_a_allocation->kilobytes_requested;
    $estimated_kb_usage += $input_b_allocation->kilobytes_requested;

    return $estimated_kb_usage;
}

sub resolve_allocation_subdirectory {
    my $self = shift;
    my $staged_basename = File::Basename::basename($self->temp_staging_directory);
    return join('/', 'build_merged_alignments', $self->id, $self->allocation_subdir_prefix . $staged_basename);
};

sub _validate_inputs {
    my $self = shift;

    my $input_dir = $self->input_directory_a;
    unless (Genome::Sys->check_for_path_existence($input_dir)) {
        $self->error_message("input_directory_a input $input_dir does not exist");
        return;
    }
    $input_dir = $self->input_directory_b;
    unless (Genome::Sys->check_for_path_existence($input_dir)) {
        $self->error_message("input_directory_b input $input_dir does not exist");
        return;
    }

    return 1;
}

sub line_count {
    my $self = shift;
    my $input = shift;
    unless( -e $input ) {
        die $self->error_message("Could not locate file for line count: $input");
    }
    my $result = `wc -l $input`; 
    my ($answer)  = split /\s/,$result;
    return $answer
}

sub _validate_output {
    my $self = shift;
    my $variant_type = $self->_variant_type;
    my $input_a_file = $self->input_directory_a."/".$variant_type.".hq.bed";
    my $input_b_file = $self->input_directory_b."/".$variant_type.".hq.bed";
    my $hq_output_file = $self->temp_staging_directory."/".$variant_type.".hq.bed";
    my $lq_output_file = $self->temp_staging_directory."/".$variant_type.".lq.bed";
    my $input_total = $self->line_count($input_a_file) + $self->line_count($input_b_file);
    my $output_total = $self->line_count($hq_output_file) + $self->line_count($lq_output_file);
    unless(($input_total - $output_total) == 0){
        die $self->error_message("Combine operation in/out check failed. Input total: $input_total \toutput total: $output_total");
    }
    return 1;
}

sub _add_as_user_of_inputs {
    my $self = shift;

    my $input_a = $self->input_a;
    my $input_b = $self->input_b;

    return (
        $input_a->add_user(user => $self, label => 'uses')
        && $input_b->add_user(user => $self, label => 'uses')
    );
}

1;
