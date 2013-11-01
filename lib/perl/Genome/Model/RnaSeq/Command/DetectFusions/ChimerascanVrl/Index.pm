package Genome::Model::RnaSeq::Command::DetectFusions::ChimerascanVrl::Index;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::ChimerascanVrl::Index {
    is => 'Command::V2',
    doc => 'create the annotation index used inside chimerascan-vrl',
    has_param => [
        version => {
            is => 'Text',
            doc => 'the version of chimerascan-vrl to use',
        },
        bowtie_version => {
            is => 'Text',
            doc => 'the version of bowtie chimerscan will use',
        },
        picard_version => {
            is => 'Text',
            doc => 'the version of picard used to manipulate BAM files',
        },
    ],
    has_input => [
        reference_build =>  {
            is => 'Genome::Model::Build::ReferenceSequence',
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
        },
        build => {
            is => 'Genome::Model::Build',
            doc => 'optional build to mark as a user of the created result',
            is_optional => 1,
        },
    ],
};

sub help_synopsis {
    return <<EOS
 genome model rna-seq detect-fusions chimerascan-vrl-index --version=0.4.6 --bowtie-version=0.12.7 --picard-version 1.2.3 --reference-build=1234 --annotation-build=1234

EOS
}

sub help_detail {
    return <<EOS
Create the annotation index used inside chimerascan-vrl
EOS
}

sub execute {
    my $self = shift;

    $self->status_message( 'Starting to create chimerascan-vrl index!' );
    my $result = Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Index->get_or_create(
            test_name => $ENV{GENOME_ALIGNER_INDEX_TEST_NAME} || undef,
            version => $self->version,
            bowtie_version => $self->bowtie_version,
            reference_build => $self->reference_build,
            annotation_build => $self->annotation_build,
            picard_version => $self->picard_version,
    );

    die( 'Failed to generate result!' ) unless $result;

    if( $self->build ){
        $self->status_message( 'Registering ' . $self->build->id . ' as a user of the generated result' );
        $result->add_user( user => $self->build, label => 'uses' );
    }

    return 1;
}

1;
