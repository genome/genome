package Genome::Model::Command::Define::ImportedVariationList;

use strict;
use warnings;

use Data::Dumper;
use Genome;

my $pp_name = "imported-variation-list";

class Genome::Model::Command::Define::ImportedVariationList {
    is => 'Genome::Model::Command::Define::Helper',
    has_input => [
        version => {
            is => 'Text',
            doc => 'The version of the build to create or update',
        },
        prefix => {
            is_optional => 1,
            is => 'Text',
            doc => 'The prefix for the name of the model to create or update (no spaces)',
        },
        model_name => {
            is_optional => 1,
            is => 'Text',
            doc => 'Override the default model name ({prefix}-{reference sequence model} by default)',
        },
    ],
    has_optional => [
        source_name => {
            is => 'Text',
            doc => 'The name of the source of the imported variants (e.g., dbsnp, 1kg)',
        },
        snv_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            doc => 'The result for snvs to import',
        },
        indel_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            doc => 'The result for indels to import',
        },
        sv_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            doc => 'The result for svs to import',
        },
        cnv_result => {
            is => 'Genome::Model::Tools::DetectVariants2::Result::Base',
            doc => 'The result for cnvs to import',
        },
        job_dispatch => {
            default_value => 'inline',
            doc => 'dispatch specification: an LSF queue or "inline"',
        },
        server_dispatch => {
            default_value => 'inline',
            doc => 'dispatch specification: an LSF queue or "inline"',
        },
        processing_profile => {
            is => 'Genome::ProcessingProfile',
            id_by => 'processing_profile_id',
            doc => 'the processing profile to use (normally selected automatically)',
        },
    ],
    has_transient_optional => [
        _reference => {
            is => 'Genome::Model::Build::ReferenceSequence',
        },
        _subject => {
            is => 'Genome::Subject',
        },
    ],
    has_transient_optional_output => [
        build => {
            is => 'Genome::Model::Build::ImportedVariationList',
            doc => 'the build created by this command',
        },
    ],
};

sub _resolve_and_validate_variant_reference_and_subject {
    my $self = shift;

    my $reference;
    for my $type ('snv', 'indel', 'sv', 'cnv') {
        my $result_accessor = $type . '_result';
        my $result = $self->$result_accessor;
        next unless $result;

        unless($reference) {
            $reference = $result->reference_build;
        } else {
            unless($reference->is_compatible_with($result->reference_build)) {
                die $self->error_message('The provided variants were not all based on the same reference sequence.');
            }
        }
    }

    unless($reference) {
        die $self->error_message('At least one result must be supplied to create an imported variation list build.');
    }

    $self->_reference($reference);
    $self->_subject($reference->model->subject);
    return 1;
}

sub execute {
    my $self = shift;

    unless(defined($self->prefix) || defined($self->model_name)) {
        $self->error_message("Please specify one of 'prefix' or 'model_name'");
        return;
    }

    if (defined($self->prefix) and ($self->prefix eq '' or $self->prefix =~ / /)) {
        $self->error_message("Invalid value for prefix '" . $self->prefix . "'. Please specify a non-empty string containing no spaces.");
        return;
    }

    # make sure we got at least one of --snv-feature-list, --indel-feature-list and 
    # verify that reference sequences are defined and match
    $self->_resolve_and_validate_variant_reference_and_subject;

    my $model = $self->_get_or_create_model();
    unless ($model) {
        $self->error_message("Failed to get or create model.");
        return;
    }
    $self->result_model_id($model->id);

    return $self->_create_build($model);
}

sub _check_existing_builds {
    my ($self, $model) = @_;

    if($model->type_name ne 'imported variation list') {
        $self->error_message("A model with the name '" . $self->model_name . "' already exists, but it is not an imported variation list.");
        return;
    }

    if ($model->reference->id != $self->_reference->id) {
        $self->error_message("Existing model '" . $model->__display_name__ . "' has reference sequence " . $model->reference->__display_name__ .
            " which conflicts with specified value of " . $self->_reference->__display_name);
        return;
    }

    if (defined $self->source_name and $model->source_name ne $self->source_name) {
        $self->error_message("Existing model '" . $model->__display_name__ . "' has source name " . $model->source_name
            . " which conflicts with specified value of " . $self->source_name);
    }

    my @builds = Genome::Model::Build::ImportedVariationList->get(model_id => $model->id, version => $self->version);
    if (scalar(@builds) > 0) {
        my $plural = scalar(@builds) > 1 ? 's' : ''; 
        $self->error_message("Existing build$plural of this model found: " . join(', ', map{$_->__display_name__} @builds));
        return;
    }

    $self->status_message('Using existing model of name "' . $model->name . '" and id ' . $model->genome_model_id . '.');

    return $model;
}

sub _get_or_create_model {
    my $self = shift;

    # * Generate a model name if one was not provided
    unless(defined($self->model_name)) {
        $self->model_name($self->prefix . "-" . $self->_reference->name);
        $self->status_message('Generated model name "' . $self->model_name . '".');
    }

    # * Make a model if one with the appropriate name does not exist.  If one does, check whether making a build for it would duplicate an
    #   existing build.
    my @models = Genome::Model->get('name' => $self->model_name);
    my $model;

    if(scalar(@models) > 1) {
        $self->error_message("More than one model (" . scalar(@models) . ") found with the name \"" . $self->model_name . "\".");
        return;
    } elsif(scalar(@models) == 1) {
        # * We're going to want a new build for an existing model, but first we should see if there are already any builds
        #   of the same version for the existing model.  If so, we ask the user to confirm that they really want to make another.
        $model = $self->_check_existing_builds($models[0]);
    } else {
        # * We need a new model
        
        my $ivl_pp = Genome::ProcessingProfile->get(name=>$pp_name);
        unless($ivl_pp){
            $self->error_message("Could not locate ImportedVariationList ProcessingProfile by name \"$pp_name\"");
            die $self->error_message;
        }

        my %create_params = (
            name => $self->model_name,
            reference => $self->_reference,
            subject_class_name => $self->_subject->class,
            subject_id => $self->_subject->id,
            processing_profile_id => $ivl_pp->id,
        );

        $create_params{source_name} = $self->source_name if $self->source_name;
        $model = Genome::Model::ImportedVariationList->create(%create_params);

        if($model) {
            if(my @problems = $model->__errors__){
                $self->error_message( "Error creating model:\n\t".  join("\n\t", map({$_->desc} @problems)) );
                return;
            } else {
                $self->status_message('Created model of name "' . $model->name . '" and id ' . $model->genome_model_id . '.');
            }
        } else {
            $self->error_message("Failed to create a new model.");
            return;
        }
    }

    return $model;
}

sub _create_build {
    my $self = shift;
    my $model = shift;

    my %build_parameters = (
        model_id => $model->id,
        version => $self->version,
    );

    for my $type ('snv', 'indel', 'sv', 'cnv') {
        my $result_accessor = $type . '_result';
        my $result = $self->$result_accessor;
        next unless $result;

        $build_parameters{$result_accessor} = $result;
    }

    my $build = Genome::Model::Build::ImportedVariationList->create(%build_parameters);
    if($build) {
        my $msg = 'Created build of id ' .$build->build_id;
        $msg .= ' with data directory "' .$build->data_directory .'".' if $build->data_directory;
        $self->status_message($msg);
        $self->build($build);
    } else {
        $self->error_message("Failed to create build for model " . $model->genome_model_id . ".");
        return;
    }

    my @dispatch_parameters;
    if(defined($self->server_dispatch)) {
        push @dispatch_parameters,
            server_dispatch => $self->server_dispatch;
    }

    if(defined($self->job_dispatch)) {
        push @dispatch_parameters,
            job_dispatch => $self->job_dispatch;
    }

    $self->status_message('Starting build.');
    if($build->start(@dispatch_parameters)) {
        $self->status_message('Started build (build is complete if it was run inline).');
    } else {
        $self->error_message("Failed to start build " . $build->build_id . " for model " . $model->genome_model_id . ".");
        return;
    }

    return 1;
}
1;
