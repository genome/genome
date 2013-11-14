package Genome::Db::Ensembl::Command::Import::CreateAnnotationStructures;

use strict;
use warnings;

use Genome;

class Genome::Db::Ensembl::Command::Import::CreateAnnotationStructures {
    is  => 'Command',
    has => [
        data_directory => {
            is => 'Path',
            doc => "ImportedAnnotation destination location",
        },
        species => {
            is => 'Text',
            doc => 'Species of annotation to import (mouse, human currently suported)',
            valid_values => [qw(mouse human)],
        },
        data_set => {
            is => 'Text',
            doc => 'Ensembl data set to import',
            default => 'Core',
        },
        log_file => {
            is => 'Path',
            doc => 'optional file to record very detailed information about transcripts, substructures, proteins and    
            genes imported',
            is_optional => 1,
        },
        dump_file => {
            is => 'Path',
            doc => 'optional file to store object cache dumps from UR.  Used to diagnose import problems',
            is_optional => 1,
        },
        build => {
            is => "Genome::Model::Build",
            id_by => "build_id",
        },
        version => {
            is => "Text",
            doc => "Version of ensembl db to import",
        },
        software_version => {
            is => 'Text',
        },
    ],
    has_input =>  [
        reference_build_id => {
            is => 'Text',
        },
        build_id => {
            is => 'Text',
            doc => 'User of the software results',
        },
    ],
};

sub help_brief
{
    "Import ensembl annotation to the file based data sources";
}

sub help_synopsis
{
    return <<EOS

EOS
}

sub help_detail
{
    return <<EOS
EOS
}

sub execute
{
    my $self = shift;
    my $result = Genome::Db::Ensembl::AnnotationStructures->get_or_create(
        reference_build_id => $self->reference_build_id,
        version => $self->version,
        species => $self->species,
        data_set => $self->data_set,
        dump_file => $self->dump_file,
        log_file => $self->log_file,
        software_version => $self->software_version,
        test_name => ($ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef),
    );

    $result->add_user(label => 'annotation_structures', user => $self->build);

    $self->status_message("Using AnnotationStructures result ".$result->id);

    foreach my $result_path ($result->result_paths) {
        Genome::Sys->create_symlink(
            $result->output_dir."/".$result_path,
            $self->data_directory."/".$result_path
        );
    }
    return 1;
}

1;

