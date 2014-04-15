package Genome::Db::Ensembl::Command::Vep;

use strict;
use warnings;
use Genome;

class Genome::Db::Ensembl::Command::Vep {
    is => 'Genome::Db::Ensembl::Command::Vep::Base',
    has => [
        ensembl_annotation_build_id => {
            is => 'String',
            doc => "ID of ImportedAnnotation build with the desired ensembl version. \n  Current default build is: $ENV{GENOME_DB_ENSEMBL_DEFAULT_IMPORTED_ANNOTATION_BUILD}, \n  Ensembl 67_37l_v2 is: 124434505)",
        },

    ],
};

sub run_with_api {
    my $self = shift;
    my %params = @_;

    $self->annotation_build->prepend_api_path_and_execute(
        %params
    );
}

sub ensembl_version {
    my $self = shift;
    return Genome::Db::Ensembl::Command::Import::Run->ensembl_version_string($self->annotation_build->ensembl_version);
}

sub annotation_build {
    my $self = shift;
    unless($self->ensembl_annotation_build_id) {
        die $self->error_message("No ensembl annotation build specified");
    }
    my $build = Genome::Model::Build::ImportedAnnotation->get($self->ensembl_annotation_build_id);

    unless ($build) {
        die $self->error_message("Could not find ImportedAnnotation build with id ".$self->ensembl_annotation_build_id);
    }
    return $build;
}

1;

