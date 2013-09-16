package Genome::Model::GenePrediction::Eukaryotic;

use strict;
use warnings;
use Genome;
use Carp 'confess';

class Genome::Model::GenePrediction::Eukaryotic {
    is => 'Genome::Model::GenePrediction',
    has_param => [
        max_bases_per_fasta => {
            is => 'Number',
            is_optional => 1,
            default => 5000000,
            doc => 'Maximum allowable base pairs in a fasta file, used for fasta splitting',
        },
        xsmall => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'If set to true, lower-case masking characters are used by repeat masker instead of N',
        },
        rfam_no_big_flag => {
            is => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => 'If set, rfam will skip big trnas',
        },
        rnammer_version => {
            is => 'Text',
            is_optional => 1,
            default => '1.2.1', 
            doc => 'Version of rnammer predictor to use',
        },
        rfamscan_version => {
            is => 'Text',
            is_optional => 1,
            valid_values => ['7.0', '8.0', '8.1', '8.1.skip_introns'],
            default => '8.1.skip_introns',
            doc => 'Version of rfamscan predictor to use',
        },
        snap_version => {
            is => 'Text',
            is_optional => 1,
            valid_values => ['2004-06-17', '2007-12-18', '2010-07-28', '2010-07-28.2'],
            default => '2010-07-28.2',
            doc => 'Version of SNAP predictor to use',
        },
        skip_masking_if_no_rna => {
            is => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => 'If set, rna sequence masking is skipped if no rna files are found. If this is false, not ' .
                   'finding an rna file is a fatal error',
        },
        skip_repeat_masker => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'If set, the repeat masker step is skipped',
        },
        exclude_overly_masked => {
            is => 'Boolean',
            is_optional => 1,
            default => 1,
            doc => 'If set, sequences that have N content greater than maximum_percent_masked are not included in masked fasta file',
        },
        maximum_percent_masked => {
            is => 'Number',
            is_optional => 1,
            default => 80,
            doc => 'If exclude_overly_masked is set, sequences with N percentage greater than this value are excluded from output file',
        },
        skip_rnammer => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'If set, rnammer prediction is skipped',
        },
        skip_trnascan => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'If set, trnascan prediction is skipped',
        },
        skip_rfamscan => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'If set, rfamscan prediction is skipped',
        },
        skip_snap => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'If set, snap prediction is skipped',
        },
        skip_fgenesh => {
            is => 'Boolean',
            is_optional => 1,
            default => 0,
            doc => 'If set, fgenesh prediction is skipped',
        },
    ],
    has_optional_input => [
        repeat_library => {
            is => 'String',
        },
        snap_models => {
            is => 'String',
        },
        fgenesh_model => {
            is => 'String',
        },
    ],
};

sub create {
    my $class = shift;
    my %params = @_;

    # Anything left in the params hash will be made into an input on the model
    my $self = $class->SUPER::create(
        name                             => delete $params{name},
        processing_profile_id            => delete $params{processing_profile_id},
        subject_id                       => delete $params{subject_id},
        subject_class_name               => delete $params{subject_class_name},
        auto_assign_inst_data            => delete $params{auto_assign_inst_data},
        auto_build_alignments            => delete $params{auto_build_alignments},
        create_assembly_model            => delete $params{create_assembly_model},
        assembly_processing_profile_name => delete $params{assembly_processing_profile_name},
        start_assembly_build             => delete $params{start_assembly_build},
        assembly_contigs_file            => delete $params{assembly_contigs_file},
    );
    return unless $self;

    # Add inputs to the model
    my $class_meta = $class->__meta__;

    for my $key (keys %params) {
        my $param_property_meta = $class_meta->property_meta_for_name($key);
        my $data_type = $param_property_meta->data_type;
        $self->add_input(
            value_class_name => 'UR::Value::' . $data_type,
            value_id => $params{$key},
            name => $key,
        );
    }

    # The snap models, fgenesh model, and repeat library all should point to existing files
    $DB::single = 1;
    if (defined $self->fgenesh_model) {
        confess "No fgenesh model file found at " . $self->fgenesh_model unless -e $self->fgenesh_model;
    }
    if (defined $self->snap_models) {
        for my $snap_model (split(",", $self->snap_models)) {
            confess "No snap model found at $snap_model" unless -e $snap_model;
        }
    }
    if (defined $self->repeat_library) {
        confess "No repeat library found at " . $self->repeat_library unless -e $self->repeat_library;
    }

    return $self;
}


sub validate_created_processing_profile {
    my $self = shift;
    my $pp = shift;
    if ($pp->skip_rnammer == 1 and $pp->skip_trnascan == 1 and $pp->skip_rfamscan == 1
            and $pp->skip_masking_if_no_rna == 0) {
        $self->error_message('All RNA predictors are disabled and the rna masking step has been ' .
            'configured to fail in the case that no RNA predictions can be found! This processing ' .
            'profile is invalid!');
        return 0;
    }
    return 1;
}

sub help_detail_for_create {
    return "create a eukaryotic gene prediction processing profile";
}

sub help_synopsis_for_create {
    return "create a eukaryotic gene prediction processing profile";
}

sub _resolve_workflow_for_build { 
    my ($self, $build) = @_;

    my $xml = __FILE__ . '.xml';
    confess "Did not find workflow xml file at $xml!" unless -e $xml;

    my $workflow = Workflow::Operation->create_from_xml($xml);
    confess "Could not create workflow object from $xml!" unless $workflow;

    $workflow->log_dir($build->log_directory);
    $workflow->name($build->workflow_name);

    return $workflow;
}

sub map_workflow_inputs {
    my ($self, $build) = @_;
    my $model = $build->model;
    confess "Could not get model from build " . $build->build_id unless $model;

    my @inputs;

    push @inputs,
        sorted_fasta => $build->sorted_fasta_file,
        domain => $self->domain,
        max_bases_per_fasta => $self->max_bases_per_fasta,
        xsmall => $self->xsmall,
        rfam_no_big_flag => $self->rfam_no_big_flag,
        rnammer_version => $self->rnammer_version,
        rfamscan_version => $self->rfamscan_version,
        snap_version =>  $self->snap_version,
        skip_masking_if_no_rna => $self->skip_masking_if_no_rna,
        repeat_library => (defined $model->repeat_library ? $model->repeat_library : '' ),
        snap_models => $model->snap_models, #FIXME using the model here is bad--use the build inputs
        fgenesh_model => $model->fgenesh_model,
        contig_fasta => $model->assembly_contigs_file,
        split_fastas_output_directory => $build->split_fastas_output_directory,
        raw_output_directory => $build->raw_output_directory,
        prediction_directory => $build->prediction_directory,
        skip_repeat_masker => $self->skip_repeat_masker,
        exclude_overly_masked => $self->exclude_overly_masked,
        maximum_percent_masked => $self->maximum_percent_masked,
        repeat_masker_ace_file => $build->repeat_masker_ace_file,
        repeat_masker_gff_file => $build->repeat_masker_gff_file,
        remove_merged_files => 1, # Don't want to keep the small unmerged files, they're 
                                  # unnecessary and clutter the data directory
        predictions_ace_file => $build->predictions_ace_file,
        overly_masked_sequence_fasta => $build->overly_masked_sequence_fasta_file,
        coding_predictions_only_flag => 1, # Similar to the rna predictions flag below, this just tells
                                           # the prediction ace file generator to only include coding gene
                                           # predictions for one particular point. The same module is used
                                           # to produce the rna gene ace file, but has a different flag set.
        rna_predictions_ace_file => $build->rna_predictions_ace_file,
        rna_predictions_only_flag => 1, # This is just used to tell the step makes the rna predictions ace
                                        # file to only look at rna. There's unfortunately no other way I'm
                                        # aware of that'll do this, and since this same module is used to
                                        # create the coding gene predictions ace file too, I can't change
                                        # the default value of the module itself.
        skip_rnammer => (defined $self->skip_rnammer ? $self->skip_rnammer : 0),
        skip_trnascan => (defined $self->skip_trnascan ? $self->skip_trnascan : 0),
        skip_rfamscan => (defined $self->skip_rfamscan ? $self->skip_rfamscan : 0),
        skip_snap => (defined $self->skip_snap ? $self->skip_snap : 0),
        skip_fgenesh => (defined $self->skip_fgenesh ? $self->skip_fgenesh : 0),
        gff_file => join('/', $build->data_directory, 'pred_gff');

    my $params;
    for (my $i = 0; $i < (scalar @inputs); $i += 2) {
        my $key = $inputs[$i];
        my $value = $inputs[$i + 1];
        $value = 'undef' if not defined $value;
        $params .= "$key : $value\n";
    }
    $self->status_message("Parameters for workflow are: \n$params");

    return @inputs;
}

sub _resolve_resource_requirements_for_build {
    # This is called during '_resolve_workflow_for_build' to specify the lsf resource requirements of the one-step
    # workflow that is generated from '_execute_build'.
    return "-M 8000000 -R 'select[type==LINUX64 && mem>8000] rusage[mem=8000]'";
}

1;

