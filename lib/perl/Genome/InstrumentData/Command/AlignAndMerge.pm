package Genome::InstrumentData::Command::AlignAndMerge;

use strict;
use warnings;
use Genome;

use List::Util qw(sum);
use POSIX qw(ceil);

class Genome::InstrumentData::Command::AlignAndMerge {
    is => ['Command::V2'],
    has => [
        name => {
            is => 'Text',
            doc => 'The aligner to use',
        },
        version => {
            is => 'Text',
            doc => 'Version of the aligner to use',
        },
        params => {
            is => 'Text',
            doc => 'Parameters to pass to the aligner',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'Id of the data to align',
        },
        reference_sequence_build => {
            is => 'Genome::Model::Build::ReferenceSequence',
            doc => 'Id of the reference to which to align',
        },
        picard_version => {
            is => 'Text',
            doc => 'The version of Picard to use when needed by aligners/filters',
        },
        samtools_version => {
            is => 'Text',
            doc => 'The version of Samtools to use when needed by aligners/filters',
        },
        result_users => {
            is => 'HASH',
            doc => 'mapping of labels to user objects. Will be added to any generated results',
        },
    ],
    has_optional_input => [
        bedtools_version => {
            is => 'Text',
            doc => 'The version of Bedtools to use when needed by aligners/filters',
        },
        force_fragment => {
            is => 'Boolean',
            default_value => 0,
            doc => 'Treat all reads as fragment reads',
        },
        trimmer_name => {
          is => 'Text',
          doc => 'name of the read trimmer to use before alignment is performed'
        },
        trimmer_params => {
          is => 'Text',
          doc => 'params to use with the read trimmmer'
        },
        trimmer_version => {
          is => 'Text',
          doc => 'version of the read trimmer to use'
        },
    ],
    has_optional_output => [
        result_id => {
            is => 'Text',
            doc => 'The ID of the result generated/found when running the command',
        },
        alignment_result => {
            is => 'Genome::InstrumentData::AlignedBamResult::Merged',
            id_by => 'result_id',
            doc => 'The result generated/found when running the command',
        },
        per_lane_alignment_result_ids => {
            is => 'Text',
            doc => 'The IDs of the per-lane results generated/found when running the command',
            is_many => 1,
        },
        per_lane_alignment_results => {
            is => 'Genome::InstrumentData::AlignmentResult',
            is_many => 1,
            calculate => q{
                return Genome::InstrumentData::AlignmentResult->get([$self->per_lane_alignment_result_ids]);
            },
            doc => 'The per-lane results generated/found when running the command',
        },
    ],
};

sub execute {
    my $self = shift;

    for my $param (qw(trimmer_name trimmer_params trimmer_version)) {
        if (defined($self->$param)) {
            $self->error_message('Optional parameter (%s) not supported at this time.', $param);
            return 0;
        }
    }

    my $result = $self->_process_alignments('get_or_create');
    unless ($result) {
        $self->error_message("Error finding or generating alignments!");
        return 0;
    }
    $self->result_id($result->id);

    my @per_lane_results = $self->_process_per_lane_alignments('get_or_create');
    unless (scalar(@per_lane_results) == scalar($self->instrument_data)) {
        $self->error_message("Error finding or generating per-lane alignments!");
        return 0;
    }
    $self->per_lane_alignment_result_ids([map { $_->id } @per_lane_results]);

    return $result;
}

sub shortcut {
    my $self = shift;

    my $result = $self->_process_alignments('get_with_lock');
    unless ($result) {
        return undef;
    }
    $self->result_id($result->id);

    my @per_lane_results = $self->_process_per_lane_alignments('get_with_lock');
    unless (scalar(@per_lane_results) == scalar($self->instrument_data)) {
        return undef;
    }
    $self->per_lane_alignment_result_ids([map { $_->id } @per_lane_results]);

    return $result;
}

sub _process_alignments {
    my $self = shift;
    my $mode = shift;

    my $merged_class = $self->class->merged_result_class($self->name);
    my %params = $self->_alignment_params;
    my $result = $merged_class->$mode(
        instrument_data => [$self->instrument_data],
        %params,
    );

    return $result;
}

sub _process_per_lane_alignments {
    my $self = shift;
    my $mode = shift;

    my $per_lane_class = $self->class->per_lane_result_class($self->name);
    my %params = $self->_alignment_params;
    my @per_lane_results;
    for my $instrument_data ($self->instrument_data) {
        push @per_lane_results, $per_lane_class->$mode(
            instrument_data => $instrument_data,
            samtools_version => $self->samtools_version,
            picard_version => $self->picard_version,
            merged_alignment_result_id => $self->result_id,
            %params,
        );
    }

    return @per_lane_results;
}

sub _alignment_params {
    my $self = shift;

    return (
        reference_build => $self->reference_sequence_build,
        users => $self->result_users,
        aligner_name => $self->name,
        aligner_version => $self->version,
        aligner_params => $self->params,
        test_name => Genome::Config::get('software_result_test_name') || undef,
    );
}

sub merged_result_class {
    my $class = shift;
    my $name = shift;
    return 'Genome::InstrumentData::AlignmentResult::Merged::' . ucfirst($name);
}

sub per_lane_result_class {
    my $class = shift;
    my $name = shift;
    return 'Genome::InstrumentData::AlignmentResult::' . ucfirst($name);
}

sub lsf_resource_string_for_aligner_and_instrument_data {
    my $class = shift;
    my $aligner_name = shift;
    my @instrument_data = @_;

    my $merged_result_class = Genome::InstrumentData::Command::AlignAndMerge->merged_result_class($aligner_name);
    my $estimated_gtmp_bytes = sum(map { $merged_result_class->estimated_gtmp_for_instrument_data($_) } @instrument_data);
    return $class->_format_lsf_resource_string($estimated_gtmp_bytes);
}

sub _format_lsf_resource_string {
    my $class = shift;
    my $gtmp_bytes = shift;

    my $cpus = 8;
    my $mem_gb = 60;
    my $queue = Genome::Config::get('lsf_queue_alignment_default');

    my $gtmp_kb = ceil($gtmp_bytes / 1024);
    my $gtmp_mb = ceil($gtmp_kb / 1024);
    my $gtmp_gb = ceil($gtmp_mb / 1024);

    my $mem_mb = $mem_gb * 1024;
    my $mem_kb = $mem_mb * 1024;

    my $select  = "select[ncpus >= $cpus && mem >= $mem_mb && gtmp >= $gtmp_gb] span[hosts=1]";
    my $rusage  = "rusage[mem=$mem_mb, gtmp=$gtmp_gb]";
    my $options = "-M $mem_kb -n $cpus -q $queue";

    my $required_usage = "-R \'$select $rusage\' $options";

    return $required_usage;
}

1;
