package Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::DetectorBase;

use strict;
use warnings;

use Genome;

class Genome::Model::RnaSeq::Command::DetectFusions::Chimerascan::DetectorBase {
    is => 'Command::V2',
    is_abstract => 1,
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
            is => "Genome::Model::Build::RnaSeq",
            is_output => 1,
        },
        reuse_bam => {
            is => 'Boolean',
            doc => 'Should we reuse the bams from alignment (experimental)',
        },
    ],
    has_optional_output => [
        software_result => {
            is => 'Genome::SoftwareResult',
        },
        bedpe_file => {
            is => 'Path',
        },
    ],
    has_param => [
        lsf_queue => {
            default_value => Genome::Config::get('lsf_queue_build_worker_alt'),
            is_optional => 1,
            doc => 'queue to use when running in a workflow',
        },
        lsf_resource => {
            default_value => "-R 'select[mem>32000 && tmp>50000] span[hosts=1] rusage[mem=32000,tmp=50000]' -M 32000000 -n 2",
            is_optional => 1,
            doc => 'default LSF resource expectations',
        },
    ],
};

sub result_class {
    die "Abstract, must be defined in subclasses";
}

sub execute {
    my $self = shift;

    unless($self->_fetch_result('get_or_create')){
        die("Unable to create a software result for " . $self->result_class);
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

    my $result_class = $self->result_class;
    my $result = $result_class->$method(
            test_name => Genome::Config::get('software_result_test_name') || undef,
            version => $self->detector_version,
            reuse_bam => $self->reuse_bam,
            bowtie_version => $self->bowtie_version,
            alignment_result => $self->build->alignment_result,
            detector_params => $self->detector_params,
            annotation_build => $self->build->annotation_build,
            picard_version => $self->build->processing_profile->picard_version,
            original_bam_paths => [map {$_->bam_path} $self->build->instrument_data],
            users => Genome::SoftwareResult::User->user_hash_for_build($self->build),
    );

    if ($result){
        $self->software_result($result);
        $self->bedpe_file($result->bedpe_file);
        return 1;
    }

    return 0;
}

1;
