package Genome::Model::Command::Define::ImportedReferenceSequence;
use strict;
use warnings;
use Genome;

class Genome::Model::Command::Define::ImportedReferenceSequence {
    is => 'Genome::Model::Command::Define::HelperDeprecated',
    has_input => [
        fasta_file => {
            is => 'Text',
            len => 1000,
            doc => "The full path and filename of the reference sequence fasta file to import."
        },
        processing_profile_id => {
            doc => 'The ID of the procesing profile to use.  There is really only one so it is typically not specified.',
            default_value => '1990904',
        },
        species_name => {
            is => 'Text',
            len => 64,
            doc => 'The species name of the reference.  This value must correspond to a species name found in the gsc.organism_taxon table.'
        },
    ],
    has_optional_input => [
        append_to => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'When specified, the newly created build will contain all sequences from the "append_to" build, followed by those from the fasta file specified.',
        },
        combine => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            is_many => 1,
            doc => 'When specified, this reference combines two or more other references',
        },
        sequence_uri => {
            is => 'Text',
            doc => 'URI to the sequence gzip file to write into BAM headers for alignments against this reference.'
        },
        use_default_sequence_uri => {
            is => 'Boolean',
            doc => 'Use a default generated URI for the BAM header.',
            default_value => 0,
        },
        assembly_name => {
            is => 'Text',
            doc => 'Assembly name to store in the SAM header.  Autoderived if not specified.',
            is_optional => 1, 
        },
        derived_from => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'The reference sequence build from which this one is derived (if any).',
        },
        prefix => {
            is => 'Text',
            doc => 'The source of the sequence, such as "NCBI".  May not have spaces.'
        },
        version => {
            is => 'Text',
            len => 128,
            doc => 'The version number and/or description of the reference.  May not have spaces.  This may be, for example '.
                   '"37" or "37_human_contamination".'
        },
        model_name => {
            is => 'Text',
            len => 255,
            doc => '$PREFIX-$SPECIES_NAME unless otherwise specified.'
        },
        subject_name => {
            is_optional => 1,            
            doc => 'Copied from species_name.'
        },
        skip_bases_files => {
            is => 'Boolean',
            is_optional => 1,
            doc => 'If set, bases files for quick sequence/position look up will not be created.  Useful for importing files with a large number of sequences',
        },
        on_warning => {
            valid_values => ['prompt', 'exit', 'continue'],
            default_value => 'prompt',
            doc => 'The action to take when emitting a warning.'
        },
        job_dispatch => {
            default_value => 'inline',
            doc => 'dispatch specification: an LSF queue or "inline"'
        },
        server_dispatch => {
            default_value => 'inline',
            doc => 'dispatch specification: an LSF queue or "inline"'
        },
        is_rederivable => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Indicates if this reference could be rederived from other internal results or if it is an external import',
        },
    ],
    has_transient => [
        result_build_id => {
            is => 'Text',
            is_optional => 1,
            doc => 'newly created build ID of reference sequence model',
        },
    ],
};

sub help_synopsis {
    return "genome model define imported-reference-sequence --species-name=human --prefix=NCBI --fasta-file=/gscuser/person/fastafile.fasta\n"
}

sub help_detail {
    #return Genome::Model::ImportedReferenceSequence->_help_detail_for_model_define;
}

sub _prompt_to_continue {
    my $self = shift;
    my $str = shift;

    if($self->on_warning eq 'exit') {
        $self->error_message($str);
        return;
    } elsif($self->on_warning eq 'continue') {
        $self->warning_message($str);
        return 1;
    } else {
        my $answer = Genome::Command::Base->_ask_user_question($str . " Continue anyway?");

        if($answer and $answer eq 'yes') {
            $self->status_message('Continuing.');
            return 1;
        } else {
            $self->status_message('Exiting.');
            return;
        }
    }
}

sub execute {
    my $self = shift;

    if ((!defined $self->sequence_uri && !$self->use_default_sequence_uri) || 
        (defined $self->sequence_uri && $self->use_default_sequence_uri)) {
        $self->error_message('Please specify one (and only one) of --sequence-uri or --use-default-sequence-uri.');
        return;
    }

    unless (-s $self->fasta_file) {
        $self->error_message('Input fasta file: '.$self->fasta_file.' is not valid.');
        return;
    }
    
    if(defined($self->prefix) && $self->prefix =~ /\s/) {
        $self->error_message("The prefix argument value must not contain spaces.");
        return;
    }

    unless(defined($self->model_name) || defined($self->species_name)) {
        $self->error_message("Either model name or species name must be supplied.  For a new model, species name is always required.");
        return;
    }

    if (defined $self->append_to and not defined $self->derived_from) {
        $self->derived_from($self->append_to);
    }

    my @combine = $self->combine();
    if (@combine) {
        if (@combine < 2) {
            die $self->error_message("The 'combine' option expects two or more existing references.  Found " . scalar(@combine));
        }
        $self->derived_from($combine[0]);
    }

    # * Verify that species name matches a taxon
    my $taxon;
    if(defined($self->species_name)) {
        my @taxons = Genome::Taxon->get('species_name' => $self->species_name);

        unless(scalar @taxons) {
            $self->error_message("No Genome::Taxon found with species name: " . $self->species_name . ".");
            return;
        }

        if(scalar(@taxons) > 1) {
            $self->error_message("Multiple Genome::Taxon instances found with species name '" . $self->species_name . "'.  This code was written " .
                   "with the assumption that species name uniquely identifies each Genome::Taxon instance.  If strain name or " .
                   "another other field is required in addition to species name to uniquely identify some Genome::Taxon instances, " .
                   "this code should be updated to take strain name or whatever other field as an argument in addition to " .
                   "species name.");
            return;
        }
        $taxon = $taxons[0];
    }

    my $model = $self->_get_or_create_model($taxon);
    return unless $model;

    $self->result_model_id($model->id);

    my $build = $self->_create_build($model);
    return unless $build;

    $self->result_build_id($build->id);

    return 1;
}

sub _check_existing_builds {
    my $self = shift;
    my $model = shift;
    my $taxon = shift;

    if($model->type_name ne 'imported reference sequence') {
        $self->error_message("A model with the name '" . $self->model_name . "' already exists, but it is not an imported reference sequence.");
        return;
    }

    if($model->subject_class_name ne 'Genome::Taxon' || $model->subject_id ne $taxon->taxon_id) {
        $self->error_message("A model with the name '" . $self->model_name . "' already exists but has a different subject class name or a " .
               "subject ID other than that corresponding to the species name supplied.");
        return;
    }

    my @builds = $model->builds;
    if(scalar(@builds) > 0 && !defined($self->version)) {
        my $continue = $self->_prompt_to_continue("At least one build already exists for the model of name '" . $model->name . "' and imported reference version was not specified.");
        return unless $continue;
    }

    my @conflicting_builds;
    foreach my $build (@builds) {
        if(defined($self->version)) {
            if(defined($build->version) && $build->version eq $self->version) {
                push @conflicting_builds, $build;
            }
        } else {
            if(not defined($build->version)) {
                push @conflicting_builds, $build;
            }
        }
    }

    if( scalar(@conflicting_builds) ) {
        my $prompt = 'One or more builds of this model already exist for this version identifier "' . ($self->version || '') . '": ' .
            join(', ', map({$_->build_id()} @conflicting_builds));
        
        my $continue = $self->_prompt_to_continue($prompt);
        return unless $continue;
    }

    $self->status_message('Using existing model of name "' . $model->name . '" and id ' . $model->genome_model_id . '.');

    return $model;
}

sub _get_or_create_model {
    my $self = shift;
    my $taxon = shift;

    # * Generate a model name if one was not provided
    unless(defined($self->model_name)) {
        my $transformed_species_name = $self->species_name;
        $transformed_species_name =~ s/\s/_/g;

        my $model_name = $transformed_species_name;

        if(defined $self->prefix) {
            $self->model_name( join('-', $self->prefix, $transformed_species_name) );
        } else {
            $self->model_name( $transformed_species_name );
        }

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
        
        $model = $self->_check_existing_builds($models[0], $taxon);
    } else {
        # * We need a new model
        
        my $irs_pp = Genome::ProcessingProfile->get(name=>"chromosome-fastas");
        unless($irs_pp){
            $self->error_message("Could not locate ImportedReferenceSequence ProcessingProfile by name \"chromosome-fastas\"");
            die $self->error_message;
        }

        $model = Genome::Model::ImportedReferenceSequence->create(
            'subject_type' => 'species_name',
            'subject_name' => $self->species_name,
            'subject_class_name' => 'Genome::Taxon',
            'subject_id' => $taxon->taxon_id,
            'processing_profile_id' => $irs_pp->id,
            'name' => $self->model_name,
            'is_rederivable' => $self->is_rederivable,
        );

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

    my @build_parameters = (
        model_id => $model->id,
        fasta_file => $self->fasta_file,
    );

    if ($self->use_default_sequence_uri) {
        push(@build_parameters, generate_sequence_uri => 1);
    }

    # single-value
    for my $build_arg (qw/sequence_uri assembly_name derived_from append_to version skip_bases_files prefix/) {
        if (defined $self->$build_arg) {
            push @build_parameters, $build_arg => $self->$build_arg;
        }
    }

    my $build = Genome::Model::Build->create(@build_parameters);
    if($build) {
        $self->status_message('Created build of id ' . $build->build_id);
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
        if($self->server_dispatch eq 'inline') {
            $self->status_message('Reference imported. ID: <' . $build->id . '>, data_directory: <' . $build->data_directory . '>.');
        } else {
            $self->status_message('Started build.');
        }
    } else {
        $self->error_message("Failed to start build " . $build->build_id . " for model " . $model->genome_model_id . ".");
        return;
    }

    return $build;
}

1;
