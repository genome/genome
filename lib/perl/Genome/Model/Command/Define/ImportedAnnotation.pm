package Genome::Model::Command::Define::ImportedAnnotation;
use strict;
use warnings;
use Genome;
use File::Basename;

class Genome::Model::Command::Define::ImportedAnnotation {
    is => 'Genome::Model::Command::Define::Helper',
    has_input => [
        species_name => {
            is => 'Text',
            doc => 'species name for annotation build (mouse, human)',
        },
        version => {
            is => 'Text',
            doc => 'Annotation version.  Generally in the form <ensembl_release_number>_<reference_number><letter_designator> (ex: 58_37c)',
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
            doc => 'Reference sequence build for the annotation version',
        },
        build_name => {
            is => 'Text',
            doc => 'human meaningful name of the build',
        },
        gtf_file => {
            is => 'Path',
            is_optional => 1,
            doc => 'gtf file for rnaSeq annotation.  ONLY use for models that are rna_seq_only',
        },
        annotation_import_version => {
            is => 'Text',
        },
    ],
    has_transient => [
        result_build_id => {
            is => 'Text',
            is_optional => 1,
            doc => 'newly created build ID of ImportedAnnotation model',
        },
    ],
};

sub help_synopsis {
    return "genome model define imported-annotation --species-name=human --version=58_37c --reference-sequence-build=GRCh37-lite-build37 --build_name=NCBI-human.ensembl/58_37c --processing-profile=imported-annotation.ensembl --model_name=NCBI-human.ensembl" 
}

sub help_detail {
    return 'creates a new annotation build';
}

sub execute {
    my $self = shift;

    my $model = $self->_get_or_create_model();
    return unless $model;

    $self->result_model_id($model->id);
    my $build = $self->_create_build($model);
    return unless $build;

    $self->result_build_id($build->id);

    return 1;
}

sub _get_or_create_model {
    my $self = shift;
    my $model_name = $self->model_name;
    my $processing_profile = $self->processing_profile;
    my $species_name = $self->species_name;
    my $model;

    #Try to find a model with the same name
    if($model_name){
        $model = Genome::Model::ImportedAnnotation->get(name => $model_name);
        if ($model){
            #sanity check the model prior to using it
            unless($model->processing_profile->id eq $processing_profile->id){
                $self->error_message("Found model with name: $model_name, but specified processing profile " . $processing_profile->id . " does not match processing profile found on model (" . $model->processing_profile->id . ")");
                return;
            }
            unless($model->species_name eq $species_name){
                $self->error_message("Found model with name: $model_name, but specified species name " . $species_name . " does not match species name found on model (" . $model->species_name . ")");
                return;
            }
            unless($model->annotation_import_version eq $self->annotation_import_version) {
                $self->error_message("Found model with name: $model_name, but specified annotation_import_version does not match annotation_import_version found on model (" .$model->annotation_import_version .") Please update the input on the model");
                return;
            }
            return $model;
        }
    }

    #Try to find a model with the same species and processing_profile (annotation source)
    $model = Genome::Model::ImportedAnnotation->get(species_name => $self->species_name, processing_profile => $self->processing_profile, annotation_import_version => $self->annotation_import_version, reference_sequence_id => $self->reference_sequence_build->model->id);
    return $model if $model;

    #Generate a name if one wasn't specified
    unless($model_name){
        $model_name = join('_', $self->reference_sequence_build->name, $species_name, $self->version);
    }   

    $model = Genome::Model::ImportedAnnotation->create(
                                                    reference_sequence => $self->reference_sequence_build->model,
                                                    name => $model_name,
                                                    processing_profile => $processing_profile,
                                                    subject => $self->reference_sequence_build->model->subject,
                                                );
    $model->species_name($self->species_name);
    $model->annotation_import_version($self->annotation_import_version);
    return $model;
}

sub _create_build {
    my $self = shift;
    my $model = shift;

    my @build_parameters = (
        model_id => $model->id,
        version => $self->version,
        reference_sequence_id => $self->reference_sequence_build->id,
    );

    my $build = Genome::Model::Build::ImportedAnnotation->get(@build_parameters);
    if($build){
        $self->error_message('Matching build already exists with id: ' . $build->id . ', exiting!');
        return;
    }
    
    push(@build_parameters, (name => $self->build_name));

    $build = Genome::Model::Build::ImportedAnnotation->create(@build_parameters);
    if($build) {
        $self->status_message('Created build with id ' . $build->build_id . ' with data directory "' . $build->data_directory . '".');
    } else {
        $self->error_message("Failed to create build for model " . $model->id . ".");
        return;
    }

    $self->status_message('Starting build.');
    my $rv;
    if ($build->processing_profile->rna_seq_only) {
        $rv = $build->start(server_dispatch => "inline");
    }
    else {
        $rv = $build->start;
    }
    if($rv){
        if ($build->processing_profile->rna_seq_only) {
            my $annotation_file_path = $build->_resolve_annotation_file_name('all_sequences','gtf',$build->reference_sequence_id,0,0);
            my $dirname = dirname($annotation_file_path);
            unless (-d $dirname) {
                Genome::Sys->create_directory($dirname);
            }
            Genome::Sys->copy_file($self->gtf_file, $annotation_file_path);
            my $bed_file_path = $build->_resolve_annotation_file_name("all_sequences", "bed", $build->reference_sequence_id, 0, 0);
            my $rv = Genome::Model::Tools::RefCov::GtfToBed->execute(
                bed_file => $bed_file_path,
                gff_file => $annotation_file_path,
            );
            unless($rv) {
                $self->error_message("Failed to create bed file from gff file");
                return;
            }
            my $squashed_bed_file_path = $build->_resolve_annotation_file_name("all_sequences-squashed", "bed", $build->reference_sequence_id, 0, 0);
            $rv = Genome::Model::Tools::BedTools::MergeByGene->execute(
                input_file => $bed_file_path,
                output_file => $squashed_bed_file_path,
            );
            unless($rv) {
                $self->error_message("Failed to create squashed bed file");
                return;
            }
            my $gff_file_path = $build->_resolve_annotation_file_name("all_sequences", "gff", $build->reference_sequence_id, 0, 0);
            Genome::Sys->create_symlink($annotation_file_path, $gff_file_path);
        }
        $self->status_message('Started build (build is complete if it was run inline).');
    } else {
        $self->error_message("Failed to start build " . $build->build_id . " for model " . $model->genome_model_id . ".");
        return;
    }

    return $build;
}

1;
