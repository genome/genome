package Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::IndexBase;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::IndexBase {
    is => 'Command::V2',
    is_abstract => 1,
    doc => 'create the annotation index used inside chimerascan',
    has_input => [
        detector_version => {
            is => 'Text',
            doc => 'the version of chimerascan to use',
        },
        detector_params => {
            is => 'Text',
            doc => 'parameters for the chosen fusion detector',
        },
        bowtie_version => {
            is => 'Text',
            doc => 'the version of bowtie chimerascan will use',
        },
        build => {
            is_output => 1,
            is => "Genome::Model::Build::RnaSeq",
        },
    ],
    has => [
        lsf_resource => {
            default_value => "-R 'select[type==LINUX64 && mem>32000] span[hosts=1] rusage[mem=32000]' -M 32000000 -n 2",
            is_param => 1,
            is_optional => 1,
            doc => 'default LSF resource expectations',
        },
    ],
};

sub result_class_name {
    die "Abstract";
}

sub execute {
    my $self = shift;

    $self->status_message( 'Starting to create index!' );
    my $result_class_name = $self->result_class_name;
    my $result = $result_class_name->get_or_create(
            test_name => $ENV{GENOME_ALIGNER_INDEX_TEST_NAME} || undef,
            version => $self->detector_version,
            bowtie_version => $self->bowtie_version,
            reference_build => $self->build->reference_sequence_build,
            annotation_build => $self->build->annotation_build,
            picard_version => $self->build->processing_profile->picard_version,
    );

    die( 'Failed to generate result!' ) unless $result;

    if( $self->build ){
        $self->status_message( 'Registering ' . $self->build->id . ' as a user of the generated result' );
        $result->add_user( user => $self->build, label => 'uses' );
    }

    return 1;
}

1;
