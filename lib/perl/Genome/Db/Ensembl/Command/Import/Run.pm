package Genome::Db::Ensembl::Command::Import::Run;

use strict;
use warnings;
use Genome;

class Genome::Db::Ensembl::Command::Import::Run {
    is => 'Command::V2',
    doc => 'Import a version of ensembl annotation',
    has => [
        data_set => {
            is => 'Text',
            doc => 'Ensembl data set to import (ex )',
            is_optional => 1,
        },
        imported_annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'Imported anntation build',
        },
        software_version => {
            is => 'Text',
        },
    ],
};

sub help_brief {
}

sub help_detail {
    return <<EOS
EOS
}

sub execute {
    my $self = shift;
    my $data_set = $self->data_set;
    my $build = $self->imported_annotation_build;
    my $data_directory = $build->data_directory;
    my $version= $build->ensembl_version;
    my $species_name = $build->species_name;
    my $reference_build_id = $build->reference_sequence_id;
    
    #download the Ensembl API to $build->data_directory
    my $api_version = $self->ensembl_version_string($version);
    my $api_result = Genome::Db::Ensembl::Api->get_or_create(
        version => $api_version,
        test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME},
    );
    #link result to build
    $api_result->add_user(label => 'api', user => $build);
    my @packages = $api_result->package_names;
    foreach my $package (@packages) {
        Genome::Sys->create_symlink(
            $api_result->output_dir."/$package",
            $data_directory."/$package",
        );
    }

    $self->status_message("Using Ensembl API result: ".$api_result->id);
    
    my $annotation_data_directory = join('/', $data_directory, 'annotation_data');
    unless(-d $annotation_data_directory){
        Genome::Sys->create_directory($annotation_data_directory);
        unless (-d $annotation_data_directory) {
            $self->error_message("Failed to create new annotation data dir: " . $annotation_data_directory);
            return;
        }
    }

    my $log_file = $data_directory . "/" . 'ensembl_import.log';
    my $dump_file = $data_directory . "/" . 'ensembl_import.dump';
    
    my $command = join(" " , 
        "genome db ensembl import create-annotation-structures", 
        "--data-directory $annotation_data_directory",
        "--version $version",
        "--species $species_name",
        "--data-set $data_set", 
        "--reference-build-id $reference_build_id",
        "--log-file $log_file",
        "--dump-file $dump_file",
        "--build-id ".$build->id,
        "--software-version ".$self->software_version);
    $build->prepend_api_path_and_execute(cmd => $command);

}

sub ensembl_version_string {
    my $self = shift;
    my $ensembl = shift;

    # <ens version>_<ncbi build vers><letter>
    # 52_36n

    my ( $e_version_number, $ncbi_build ) = split( /_/x, $ensembl );
    return $e_version_number;
}

1;
