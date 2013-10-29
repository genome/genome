package Genome::Model::RnaSeq::Command::DetectFusions::ChimerascanVrl;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::ChimerascanVrl {
    is => 'Genome::Model::RnaSeq::Command::DetectFusions::Base',
    doc => 'run the chimerascan-vrl transcript fusion detector',
    has => [
        lsf_resource => {
            default_value => "-R 'select[type==LINUX64 && mem>32000] span[hosts=1] rusage[mem=32000]' -M 32000000 -n 2",
            is_param => 1,
            is_optional => 1,
            doc => 'default LSF resource expectations',
        },
    ],
};

sub help_synopsis {
    return <<EOS
 genome model rna-seq detect-fusions chimerascan-vrl build-id=1234
EOS
}

sub help_detail {
    return <<EOS
Run the chimerascan-vrl gene fusion detector.

It is used by the RNASeq pipeline to perform fusion detection when the fusion detection strategy is set to something like:
 'chimerascan-vrl 0.4.6 [-p 2 --reuse-bam=1 --bowtie-version=0.12.7]'

EOS
}

sub execute {
    my $self = shift;

    unless($self->_fetch_result('get_or_create')){
        die("Unable to create a software result for ChimerascanVrl");
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

    my $result = Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::VariableReadLength::Result->$method(
            test_name => $ENV{GENOME_SOFTWARE_RESULT_TEST_NAME} || undef,
            version => $self->version,
            alignment_result => $self->build->alignment_result,
            detector_params => $self->detector_params,
            annotation_build => $self->build->annotation_build,
            picard_version => $self->build->processing_profile->picard_version,
            original_bam_paths => [map {$_->bam_path} $self->build->instrument_data],
    );

    if ($result){
        $self->_link_build_to_result($result);
        return 1;
    }

    return 0;
}

1;
