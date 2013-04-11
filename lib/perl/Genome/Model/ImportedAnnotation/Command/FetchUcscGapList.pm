package Genome::Model::ImportedAnnotation::Command::FetchUcscGapList;

use strict;
use warnings;

use Genome;

class Genome::Model::ImportedAnnotation::Command::FetchUcscGapList {
    is => 'Command::V2',
    has_input => [
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'Annotation build under construction',
        },
    ],
    has_output => [
        gap_feature_list => {
            is => 'Genome::FeatureList',
        },
    ],
    doc => "Fetches gap data for a reference annotation build",
};

sub execute {
    my $self = shift;

    my ($reference_sequence, $reference_name) = $self->_resolve_reference();

    my $gap_filename = Genome::Sys->create_temp_file_path();
    Genome::Db::Ucsc::Command::GapList->execute(
        filename => $gap_filename,
        reference_name => $reference_name,
    );
    my $file_md5 = Genome::Sys->md5sum($gap_filename);

    $self->gap_feature_list(Genome::FeatureList->create(
        name => $self->_resolve_feature_list_name(),
        format => 'true-BED',

        file_content_hash => $file_md5,
        reference => $reference_sequence,
        file_path => $gap_filename,

        description => 'Downloaded from UCSC database',
        source => 'UCSC',
    ));

    return $self->gap_feature_list;
}

sub _resolve_reference {
    my $self = shift;

    my $reference_sequence = $self->annotation_build->reference_sequence;
    unless ($reference_sequence) {
        $self->error_message(sprintf(
                "Failed to get reference sequence for annotation build %s",
                $self->annotation_build->id));
        die $self->error_message;
    }

    my $reference_name = $reference_sequence->name;
    unless ($reference_name) {
        $self->error_message(sprintf(
                "No reference name associated with reference %s for annoation build %s",
                $reference_sequence->id, $self->annoation_build->id));
        die $self->error_message;
    }
    if ($reference_name eq "GRCh37-lite-build37" or $reference_name eq "g1k-human-build37") {
        $reference_name = "hg19";
    }
    elsif ($reference_name eq "NCBI-human-build36") {
        $reference_name = "hg18";
    }

    return ($reference_sequence, $reference_name)
}

sub _resolve_feature_list_name {
    my $self = shift;
    return $self->annotation_build->name . " - gap list";
}

1;
