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
    has_optional_output => [
        software_result => {
            is => 'Genome::SoftwareResult',
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

    unless($self->_fetch_result('get_or_create')){
        die("Unable to create a software result for " . $self->result_class_name);
    }

    return 1;
}

sub shortcut {
    my $self = shift;
    return $self->_fetch_result('get_with_lock');
}

sub _fetch_result {
    my $self = shift;
    my $method = shift;

    my $result_class_name = $self->result_class_name;
    my $result = $result_class_name->$method(
            test_name => $ENV{GENOME_ALIGNER_INDEX_TEST_NAME} || undef,
            version => $self->detector_version,
            bowtie_version => $self->bowtie_version,
            reference_build => $self->build->reference_sequence_build,
            annotation_build => $self->build->annotation_build,
            picard_version => $self->build->processing_profile->picard_version,
    );

    if ($result){
        $self->software_result($result);
        return 1;
    }

    return 0;
}

1;
