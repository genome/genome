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
        source_md5 => {
            is => 'Text',
        },
    ],
    has_transient => [
        source_bed => {
            is => 'Path',
            is_optional => 1,
        },
        _target_bed => {
            is_calculated => 1,
            calculate => q| return 'converted.bed' |,
        }
    ],
};

sub target_bed {
    my $self = shift;
    my $bed = File::Spec->join($self->disk_allocations->absolute_path, $self->_target_bed);
    unless (-s $bed) {
        die $self->error_message("When trying to access target_bed ($bed) , we found that it does not exist or has no size");
    }
    return $bed;
};

sub _generate_result {
    my ($self, $staging_directory) = @_;

    $self->_validate_source_bed;

    my $converted_file_path = File::Spec->join($staging_directory, $self->_target_bed);
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

sub _validate_source_bed {
    my $self = shift;

    my $bed_md5 = Genome::Sys->md5sum($self->source_bed);
    unless($bed_md5 eq $self->source_md5) {
        die $self->error_message(
            'Source BED <%s> failed MD5 check. Expected %s, got %s.',
            $self->source_bed, $self->source_md5sum, $bed_md5,
        );
    }

    return 1;
}

1;

