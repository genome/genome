package Genome::Model::GenePrediction::Bacterial;

use strict;
use warnings;

use Genome;
use Carp "confess";
use File::Spec;

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
        locus_id => {
            via => 'subject',
            to => 'locus_tag',
        },
    ],
    has_optional_input => [
        dev => {
            is => 'Boolean',
            doc => 'If set, dev databases are used instead of production databases',
        },
        run_type => {
            is => 'String', # TODO Does this affect processing? Why do we need to note it?
            doc => 'A three letter identifier appended to locus id, (DFT, FNL, etc)',
        },
        locus_suffix => {
            is => 'Text',
            doc => 'appended to the locus_id (and included in the locus tag). (Used for testing.)',
        },
        assembly_version => {
            is => 'String', # TODO Can this be removed or derived from the assembly in some way?
            doc => 'This notes the assembly version, but doesn\'t really seem to change...',
        },
        project_type => {
            is => 'String', # TODO What is this? Why do we need it?
            doc => 'The type of project this data is being generated for (HGMI, for example)',
        },
        pipeline_version => {
            is => 'String', # TODO Can this be removed? Why do we need it?
            doc => 'Apparently, this notes the pipeline version.',
        },
        acedb_version => {
            is => 'String', # TODO If we can figure out a way to automate switching to a new db, this can go away
            doc => 'Notes the version of aceDB that results should be uploaded to',
        },
        nr_database_location => {
            is => 'FilePath', # TODO Once using local NR is fully tested and trusted, this param can be removed
            doc => 'The NR database that should be used by default, may be overridden by local copies',
        },
    ],
};

sub create {
    my $class = shift;

    my ($bx, @extra) = $class->define_boolexpr(@_);
    my %params = @extra;

    #Grab the other options specific to creating GenePrediction models
    # Anything left in the params will be made into an input on the model
    my %additional_create_options;
    for my $option ('create_assembly_model', 'assembly_processing_profile_name', 'start_assembly_build', 'assembly_contigs_file') {
        if(exists $params{$option}) {
            $additional_create_options{$option} = delete $params{$option};
        }
    }

    my $self = $class->SUPER::create($bx->params_list, %additional_create_options);
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

    $self->fatal_message('This build type is no longer supported.');
}

sub _resolve_resource_requirements_for_build {
    # This is called during '_resolve_workflow_for_build' to specify the lsf resource requirements of the one-step
    # workflow that is generated from '_execute_build'.
    my ($build) = @_;
    return "-R 'rusage[tmp=10000:mem=2000]' -M 2000000";
}

1;

