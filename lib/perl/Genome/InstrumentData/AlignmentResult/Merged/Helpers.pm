package Genome::InstrumentData::AlignmentResult::Merged::Helpers;

use strict;
use warnings;
use Genome;
use Sys::Hostname;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    create_bam_md5
    resolve_allocation_subdirectory
    resolve_alignment_subdirectory
    resolve_allocation_disk_group_name
);

sub create_bam_md5 {
    my $self = shift;

    my $bam_file = shift;
    my $md5_file = $bam_file.'.md5';
    my $cmd = "md5sum $bam_file > $md5_file";

    $self->debug_message("Creating md5 file for the BAM file...");

    Genome::Sys->shellcmd(
        cmd                        => $cmd,
        input_files                => [$bam_file],
        output_files               => [$md5_file],
        skip_if_output_is_present  => 0,
    );

    return 1;
}

sub resolve_allocation_subdirectory {
    return $_[0]->resolve_alignment_subdirectory;
}

sub resolve_alignment_subdirectory {
    my $self = shift;

    my $hostname = hostname;
    my $user = $ENV{'USER'};
    my $base_dir = sprintf("merged-alignment-%s-%s-%s-%s", $hostname, $user, $$, $self->id);
    # TODO: the first subdir is actually specified by the disk management system.
    my $directory = join('/', 'build_merged_alignments', $base_dir);
    return $directory;
}

sub resolve_allocation_disk_group_name {
    $ENV{GENOME_DISK_GROUP_MODELS};
}

1;
