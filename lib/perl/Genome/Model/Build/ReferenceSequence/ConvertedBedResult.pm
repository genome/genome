package Genome::Model::Build::ReferenceSequence::ConvertedBedResult;

use strict;
use warnings;
use Genome;

class Genome::Model::Build::ReferenceSequence::ConvertedBedResult {
    is => 'Genome::SoftwareResult::DiskAllocationStaged',
    has_input => [
        source_reference => {
            is => 'Genome::Model::Build::ReferenceSequence',
        },
        target_reference => {
            is => 'Genome::Model::Build::ReferenceSequence',
        },
        source_bed => {
            is => 'Path',
        },
    ],
    has => [
        target_bed => {
            is_calculated => 1,
            calculate_from => ['disk_allocations', '_target_bed'],
            calculate => q|
                return File::Spec->join($disk_allocations->absolute_path, $_target_bed);
            |,
            doc => 'Final converted bed path',
        },
    ],
    has_transient => [
        _target_bed => {
            is_calculated => 1,
            calculate => q| return 'converted.bed' |,
        }
    ],
};

sub _generate_result {
    my ($self, $staging_directory) = @_;
    my $converted_file_path = File::Spec->join($staging_directory, $self->_target_bed);
    $DB::single=1;
    my $converted_bed_file = Genome::Model::Build::ReferenceSequence::Converter->convert_bed(
        $self->source_bed, $self->source_reference, $converted_file_path, $self->target_reference);

    unless (-s $converted_bed_file) {
        die $self->error_message("Failed to convert from bed (%s) with reference (%s) to bed (%s) with reference (%s)",
            $self->source_bed, $self->source_reference, $converted_file_path, $self->target_reference);
    }
    return 1;
}

# we are going to need about as much space for the target_bed as the source_bed used
sub _staging_kilobytes_requested {
    my $self = shift;
    my @stat = stat($self->source_bed);
    my $size = $stat[7]; # in bytes
    return ($size / 1024);
}

1;

