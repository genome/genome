package Genome::InstrumentData::AlignmentResult::Merged::External;

use strict;
use warnings;

use Genome;

use File::Basename qw();
use File::Spec;


class Genome::InstrumentData::AlignmentResult::Merged::External {
    is => ['Genome::SoftwareResult::StageableSimple', 'Genome::InstrumentData::AlignedBamResult::Merged'],
    has_input => [
        sample => {
            is => 'Genome::Sample',
            doc => 'The sample from which the data for these alignments were taken',
        },
        file_content_hash => {
            is => 'Text',
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'The instrument data used for these alignments (if any)',
            is_optional => 1,
            is_many => 1,
        },
        annotation_build => {
            is => 'Genome::Model::Build::ImportedAnnotation',
            doc => 'The annotation build used for these alignments (if any)',
            is_optional => 1,
        },
    ],
    has_optional_metric => [
        original_file_path => {
            is => 'Text',
        },
        description => {
            is => 'Text',
            doc => 'brief information about the alignment data',
        },
        filetype => {
            is => 'Text',
            doc => 'the type of alignment (expressed as file extension)',
        },
    ],
    has_constant => [
        bam_path => {
            is_output => 1,
            calculate_from => [qw(id filetype output_dir)],
            calculate => q{
                return File::Spec->join($output_dir, ($id . $filetype))
            },
        },
    ],
    has => [
        sample_name => {
            via => 'sample',
            to => 'name',
        },
    ],
};

sub _run {
    my $self = shift;

    my $original_file = $self->original_file_path;
    my (undef, undef, $suffix) = File::Basename::fileparse($original_file, '.bam', '.cram');
    $self->filetype($suffix);

    my $staging_dir = $self->temp_staging_directory;
    my $staging_file = File::Spec->join($staging_dir, $self->id . $suffix);

    Genome::Sys->create_symlink($original_file, $staging_file);

    return 1;
}

sub _needs_symlinks_followed_when_syncing {
    return 1;
}

sub resolve_allocation_subdirectory {
    my $self = shift;

    return File::Spec->join('model_data', 'external-merged-alignment', $self->id);
}


sub _modify_params_for_lookup_hash {
    my ($class, $params_ref) = @_;

    my $original_file_path = delete $params_ref->{'original_file_path'};
    my $specified_checksum = $params_ref->{'file_content_hash'};
    $params_ref->{'file_content_hash'} = $class->_calculate_and_compare_md5_hashes(
        $original_file_path, $specified_checksum);
}

sub _calculate_and_compare_md5_hashes {
    my ($class, $original_file_path, $specified_checksum) = @_;

    my $checksum = $specified_checksum;
    if (defined($original_file_path) and -e $original_file_path) {
        $checksum = Genome::Sys->md5sum($original_file_path);
        if (defined($specified_checksum) and $specified_checksum ne $checksum) {
            $class->fatal_message(
                'file_content_hash does not match md5sum output for original_file_path');
        }
    }

    return $checksum;
}

sub _gather_params_for_get_or_create {
    my $class = shift;

    my $bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);

    if($bx->specifies_value_for('original_file_path')) {
        my $original_file_path = $bx->value_for('original_file_path');
        my $specified_checksum;
        if ($bx->specifies_value_for('file_content_hash')) {
            $specified_checksum = $bx->value_for('file_content_hash');
        }
        $bx = $bx->add_filter('file_content_hash',
            $class->_calculate_and_compare_md5_hashes(
                $original_file_path, $specified_checksum));
    }

    return $class->SUPER::_gather_params_for_get_or_create($bx->params_list);
}

1;
