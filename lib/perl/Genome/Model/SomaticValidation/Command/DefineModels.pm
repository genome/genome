package Genome::Model::SomaticValidation::Command::DefineModels;

use strict;
use warnings;

use Genome;

class Genome::Model::SomaticValidation::Command::DefineModels {
    is => 'Command::V2',
    has_input => [
        models => {
            is => 'Genome::Model',
            doc => 'the somatic models from discovery',
            is_many => 1,
        },
        new_model_group => {
            is => 'Genome::ModelGroup',
            is_optional => 1,
            doc => 'the group to which to add the new models',
        },
        design => {
            is => 'Genome::FeatureList',
            doc => 'The designs sent to the vendor',
        },
        target => {
            is => 'Genome::FeatureList',
            doc => 'The target set received back from the vendor',
        },
        region_of_interest_set => {
            is => 'Genome::FeatureList',
            is_optional => 1,
            doc => 'The regions on which to perform refcov analysis (defaults to the target)',
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            is_optional => 1,
            doc => 'A processing profile to use for all models defined (the current default will be used otherwise)',
        },
        variant_file_list => {
            is => 'Text',
            is_optional => 1,
            doc => 'A file listing the variants for each patient',
        },
        variant_file_format => {
            is => 'Text',
            is_optional => 1,
            default_value => 'bed',
            doc => 'format of the files listed in the variant_file_list, if provided',
        },
        generate_variant_lists => {
            is => 'Boolean',
            is_optional => 1,
            default_value => 0,
            doc => 'In lieu of provided variant lists, automatically generate lists based on the intersection of the variants in the discovery models and the target',
        },
    ],
    has_output => [
        result_models => {
            is => 'Genome::Model::SomaticValidation',
            doc => 'the models created',
            is_many => 1,
        },
        result_model_ids => {
            is => 'Text',
            is_many => 1,
        },
    ],
    doc => 'define models based on an existing model-group',
};

sub sub_command_category { 'analyst tools' }

sub execute {
    my $self = shift;

    my @models = $self->models;

    my @params;

    for my $param ('processing_profile', 'region_of_interest_set', 'target', 'design') {
        push @params, $param => $self->$param
            if $self->$param;
    }

    if(!$self->region_of_interest_set and $self->target) {
        $self->region_of_interest_set($self->target);
    }

    my $reference_sequence_build = $self->region_of_interest_set->reference;
    push @params, reference_sequence_build => $reference_sequence_build;

    my %variants = $self->_generate_variant_mapping;

    my @new_m;
    for my $model (@models) {
        my $patient = $model->subject;
        if($patient->isa('Genome::Sample')) { $patient = $patient->source; }

        my @results;
        if($self->variant_file_list) {
            ACCESSOR: for my $accessor ('id', 'name', 'common_name') {
                if(exists $variants{$patient->$accessor}) {
                    for my $variant_type (keys %{ $variants{$patient->$accessor} }) {
                        my $result = $self->_upload_result($model, $variant_type, @{ $variants{$patient->$accessor}{$variant_type} });
                        push @results, $result;
                    }
                    last ACCESSOR; #found it already
                }
            }
        } elsif($self->generate_variant_lists) {
            my $extract_variants_cmd = Genome::Model::SomaticValidation::Command::ExtractVariantsToValidate->create(
                model => $model,
                target => $self->target,
            );
            $extract_variants_cmd->dump_status_messages(1);
            unless($extract_variants_cmd->execute()) {
                die $self->error_message('Failed to extract variants for model ' . $model->__display_name__);
            }

            for my $type ('snv', 'indel', 'sv', 'cnv') {
                my $accessor = $type . '_variant_list';
                my $result = $extract_variants_cmd->$accessor;

                push @results, $result if $result;
            }
        }

        my ($tumor_sample, $normal_sample);
        if($model->can('tumor_model')) {
            $tumor_sample = $model->tumor_model->subject;
        } else {
            $tumor_sample = $model->tumor_sample;
        }
        if($model->can('normal_model')) {
            $normal_sample = $model->normal_model->subject;
        } else {
            $normal_sample = $model->normal_sample;
        }

        my $define_cmd = Genome::Model::Command::Define::SomaticValidation->create(
            @params,
            (scalar @results?
                () :
                (tumor_sample => $tumor_sample,
                $normal_sample?
                    (normal_sample => $normal_sample) :
                    ()
                )
            ),
            variants => \@results,
        );

        $define_cmd->dump_status_messages($self->dump_status_messages);
        unless($define_cmd->execute) {
            die $self->error_message('Failed to create validation model based on ' . $model->__display_name__);
        }
        push @new_m, $define_cmd->result_models;
    }

    if($self->new_model_group) {
        $self->new_model_group->assign_models(@new_m);
        $self->debug_message('Model group updated: ' . $self->new_model_group->__display_name__);
    }

    $self->result_models(\@new_m);
    $self->result_model_ids([map($_->id, @new_m)]);

    return 1;
}

#TODO Just defer to import-variants command
sub _generate_variant_mapping {
    my $self = shift;

    return unless $self->variant_file_list;

    my %data;

    my $variant_file_list_fh = Genome::Sys->open_file_for_reading($self->variant_file_list);

    my $variant_type;
    while(my $line = <$variant_file_list_fh>) {
        chomp $line;
        if($line =~ m{.*/((?:\w+\d+)|(?:\w_\w\w-[^/]+))/[^/]+$}) {
            my $patient = $1;

            $data{$patient}{$variant_type} ||= [];
            push @{ $data{$patient}{$variant_type} }, $line;
        } elsif($line =~ m{.*/([^/.]+)(?:\.[^./]+)+$}) {
            my $patient = $1;

            $data{$patient}{$variant_type} ||= [];
            push @{ $data{$patient}{$variant_type} }, $line;
        } elsif(grep($_ eq $line, 'snvs', 'indels', 'svs')) {
            $variant_type = $line; #header indicating SNVs, indels, or SVs
            $variant_type =~ s/s$//;
        } else {
            die $self->error_message('Could not determine patient for this file: ' . $line);
        }
    }

    return %data;
}

sub _upload_result {
    my $self = shift;
    my $model = shift;
    my $variant_type = shift;
    my @files = @_;

    my $file;
    if(scalar(@files) > 1) {
        $file = Genome::Sys->create_temp_file_path;
        Genome::Sys->cat(
            input_files => \@files,
            output_file => $file,
        );
    } elsif (scalar(@files) == 1) {
        $file = $files[0];
    } else {
        die '_upload_result called with no variant files';
    }

    my $build = $model->last_complete_build;
    unless($build) {
        die $self->error_message('No complete build found for model ' . $model->__display_name__);
    }

    my $result_cmd = Genome::Model::SomaticValidation::Command::ManualResult->get_or_create(
        source_build => $build,
        variant_file => $file,
        variant_type => $variant_type,
        format => $self->variant_file_format,
        description => 'generated from ' . $self->variant_file_list,
    );

    unless($result_cmd->execute) {
        die $self->error_message('Failed to create result for model ' . $model->__display_name__ . ' ' . $variant_type . 's using files: ' . join(', ', @files));
    }

    return $result_cmd->manual_result;
}

1;

