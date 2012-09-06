package Genome::Model::GenePrediction::Bacterial;

use strict;
use warnings;

use Genome;
use Carp "confess";

class Genome::Model::GenePrediction::Bacterial {
    is => 'Genome::Model::GenePrediction',
    has_param => [
        skip_core_gene_check => {
            is => 'Boolean',
            doc => 'If set, the core gene check is not performed',
            is_optional => 1,
            default => 0,
        },
        minimum_sequence_length => {
            is => 'Number',
            doc => 'Minimum contig sequence length',
            is_optional => 1,
            default => 200,
        },
        runner_count => {
            is => 'Number',
            doc => 'Number of runners for the gene prediction step',
            is_optional => 1,
            default => 50,
        }, 
        skip_acedb_parse => {
            is => 'Boolean',
            doc => 'If set, skip aceDB parsing in bap project finish',
            is_optional => 1,
        },
        keggscan_version => {
            is => 'Number',
            doc => 'Version of KEGGScan to use',
            is_optional => 1,
        },
        interpro_version => {
            is => 'Text',
            doc => 'Version of iprscan and data to use',
            is_optional => 1,
        },
    ],
    has_optional => [
        dev => {
            is => 'Boolean',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'dev' ],
            doc => 'If set, dev databases are used instead of production databases',
        },
        locus_id => {
            via => 'subject',
            to => 'locus_tag',
        },
        run_type => {
            is => 'String', # TODO Does this affect processing? Why do we need to note it?
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'run_type' ],
            doc => 'A three letter identifier appended to locus id, (DFT, FNL, etc)',
        },
        locus_suffix => {
            is => 'Text',
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'locus_suffix', value_class_name => 'UR::Value::Text' ],
            doc => 'appended to the locus_id (and included in the locus tag). (Used for testing.)',
            is_mutable => 1,
            is_many => 0,
        },
        assembly_version => {
            is => 'String', # TODO Can this be removed or derived from the assembly in some way?
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'assembly_version' ],
            doc => 'This notes the assembly version, but doesn\'t really seem to change...',
        },
        project_type => {
            is => 'String', # TODO What is this? Why do we need it?
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'project_type' ],
            doc => 'The type of project this data is being generated for (HGMI, for example)',
        },
        pipeline_version => {
            is => 'String', # TODO Can this be removed? Why do we need it?
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'pipeline_version' ],
            doc => 'Apparently, this notes the pipeline version.', 
        },
        acedb_version => {
            is => 'String', # TODO If we can figure out a way to automate switching to a new db, this can go away
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'acedb_version' ],
            doc => 'Notes the version of aceDB that results should be uploaded to',
        },
        nr_database_location => {
            is => 'Path', # TODO Once using local NR is fully tested and trusted, this param can be removed
            via => 'inputs',
            to => 'value_id',
            where => [ name => 'nr_database_location' ],
            doc => 'The NR database that should be used by default, may be overridden by local copies',
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
    for my $key (keys %params) {
        $self->add_input(
            value_class_name => 'UR::Value',
            value_id => $params{$key},
            name => $key,
        );
    }

    return $self;
}

sub _execute_build {
    my ($self, $build) = @_;

    $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1;

    my $model = $build->model;
    $self->status_message("Executing build logic for " . $self->__display_name__ . ":" . $build->__display_name__);

    # TODO Make explicit build links between this build and the assembly build for tracking

    my $config_file_path = $build->create_config_file;
    unless (-s $config_file_path) {
        $self->error_message("Configuration file not found at expected location: $config_file_path");
        confess;
    }

    $self->status_message("Configuration file created at $config_file_path, creating hap command object");

    Genome::Sys->create_symlink('/gscmnt/278/analysis/HGMI/Acedb', $build->data_directory . '/Acedb'); #FIXME Replace this hardcoded path

    my $hap_object = Genome::Model::Tools::Hgmi::Hap->create(
        config => $config_file_path,
        dev => $model->dev,
        skip_core_check => $self->skip_core_gene_check,
        skip_protein_annotation => (not ($self->keggscan_version || $self->interpro_version)), # TODO Include protein annotation if we have a pp. params suggesting we want it
        keggscan_version => $self->keggscan_version,
        interpro_version => $self->interpro_version,
    );
    unless ($hap_object) {
        $self->error_message("Could not create hap command object!");
        confess;
    }

    $self->status_message("Hap command object created, executing!");

    # THIS IS IMPORTANT! Hap creates forked processes as part of the prediction step, and these child
    # processes get a REFERENCE to open db handles, which get cleaned up and closed during cleanup of
    # the child process. This causes problems in this process, because it expects the handle to still be
    # open. Attempting to use that handle results in frustrating errors like this:
    # DBD::Oracle::db rollback failed: ORA-03113: end-of-file on communication channel (DBD ERROR: OCITransRollback)
    if (Genome::DataSource::GMSchema->has_default_handle) {
        $self->status_message("Disconnecting GMSchema default handle.");
        Genome::DataSource::GMSchema->disconnect_default_dbh();
    }

    my $hap_rv = $hap_object->execute;
    unless ($hap_rv) {
        $self->error_message("Trouble executing hap command!");
        confess;
    }

    $self->status_message("Hap executed and no problems detected!");
    return 1;
}

1;

