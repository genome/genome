package Genome::InstrumentData::AlignmentResult::Merged::Speedseq;

use strict;
use warnings;
use Genome;
use Genome::InstrumentData::AlignmentResult::Merged::Helpers qw(
    create_bam_md5
    resolve_allocation_subdirectory
    resolve_alignment_subdirectory
    resolve_allocation_disk_group_name
);
use File::stat;
use Genome::Utility::Text;

class Genome::InstrumentData::AlignmentResult::Merged::Speedseq {
    is => ['Genome::InstrumentData::AlignedBamResult::Merged', 'Genome::SoftwareResult::WithNestedResults'],
    has => [
        aligner => {
            calculate_from => [qw/aligner_name aligner_version aligner_params/],
            calculate => q|no warnings; "$aligner_name $aligner_version $aligner_params"|
        },
    ],
    has_input => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
        },
        reference_build => {
            is => 'Genome::Model::Build::ImportedReferenceSequence',
        },
        aligner_name => {
            is => 'Text',
            doc => 'the name of the aligner to use, maq, blat, newbler etc.',
        },
        aligner_version => {
            is => 'Text',
            doc => 'the version of the aligner to use, i.e. 0.6.8, 0.7.1, etc.',
            is_optional=>1,
        },
        aligner_params => {
            is => 'Text',
            is_optional=>1,
            doc => 'any additional params for the aligner in a single string',
        },
    ],
    has_param => [
        samtools_version => {
            is=>'Text',
            is_optional=>1,
            doc=>'Version of samtools to use when creating BAM files',
        },
        picard_version => {
            is=>'Text',
            is_optional=>1,
            doc=>'Version of picard to use when creating bam files',
        },
    ],
    has_calculated => [
        _final_bam_file => {
            is => 'Text', calculate_from => ['temp_staging_directory', 'id',],
            calculate => q{ return join('/', $temp_staging_directory, $id . '.bam'); },
        },
        merged_alignment_bam_path => {
            is => 'Text', calculate_from => ['output_dir', 'id'],
            calculate => q{ return join('/', $output_dir, $id . '.bam'); }
        },
    ]
};

sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_);
    return unless ($self);

    $self->_prepare_staging_directory;
    my $temp_directory = Genome::Sys->create_temp_directory();

    my $aligner_index = Genome::InstrumentData::AlignmentResult::get_reference_sequence_index($self);

    my %params = (
        bams => [ map {$_->bam_path} $self->instrument_data ],
        reference_fasta => $aligner_index->full_consensus_path('fa'),
        output_prefix => File::Spec->join($self->temp_staging_directory, $self->id),
        version => $self->aligner_version,
        temp_directory => $temp_directory,
    );
    my %aligner_params = $class->_parse_aligner_params($self->aligner_params);
    $aligner_params{threads} = $self->_available_cpu_count;
    my $command = Genome::Model::Tools::Speedseq::Realign->create(%params, %aligner_params);
    unless ($command->execute) {
        die $self->error_message('Failed to execute Speedseq realign for instrument data: ' . join(', ', map {$_->id} $self->instrument_data));
    }

    my $final_bam = $self->_final_bam_file,
    $self->debug_message("Indexing the final BAM file...");
    my $index_cmd = Genome::Model::Tools::Sam::IndexBam->create(
        bam_file    => $final_bam,
        use_version => $self->samtools_version,
    );
    my $index_cmd_rv = $index_cmd->execute;

    if($index_cmd_rv ne 1) {
        #not failing here because this is not a critical error.  this can be regenerated manually if needed.
        $self->warning_message('Failed to create bam index for ' . $final_bam);
    }

    $self->create_bam_md5($final_bam);

    $self->_prepare_output_directory;
    $self->_promote_data;
    $self->_reallocate_disk_allocation;

    return $self;
}

sub aligner_name_for_aligner_index {
    my $self = shift;
    return $self->aligner_name;
}

sub estimated_gtmp_for_instrument_data {
    my $class = shift;
    my $instrument_data = shift;
    my $bam_path = $instrument_data->bam_path();

    my $st = stat($bam_path);
    unless ($st) {
        $class->warning_message('Unable to find bam file at %s', $bam_path);
        return 1; #This job will fail when scheduled if the BAM is missing, so won't need much gtmp!
    }

    return 3 * $st->size();
}

sub _modify_params_for_lookup_hash {
    my $class = shift;
    my $params_ref = shift;

    my $aligner_param_str = delete $params_ref->{aligner_params};
    return unless $aligner_param_str;

    my %aligner_params = $class->_parse_aligner_params($aligner_param_str);

    delete $aligner_params{sort_memory};
    delete $aligner_params{verbose};
    delete $aligner_params{threads};

    $aligner_param_str = join(',', map(
        join(' => ', $_, Genome::Utility::Text::wrap_as_string($aligner_params{$_})),
        sort keys %aligner_params)
    );
    $params_ref->{aligner_params} = $aligner_param_str;

    return $params_ref;
}

sub _parse_aligner_params {
    my $class = shift;
    my $param_str = shift;

    local $@;
    my %params_parsed = eval($param_str);
    my $error = $@;
    if($error) {
        die $class->error_message('Failed to parse aligner params: %s', $error);
    }

    return %params_parsed;
}

1;
