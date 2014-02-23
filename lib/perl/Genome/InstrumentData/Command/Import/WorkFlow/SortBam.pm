package Genome::InstrumentData::Command::Import::WorkFlow::SortBam;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::WorkFlow::SortBam { 
    is => 'Command::V2',
    has_input => [
        unsorted_bam_path => {
            is => 'Text',
            doc => 'The path of the unsorted bam to sort.',
        }
    ],
    has_output => [ 
        sorted_bam_path => {
            is => 'Text',
            calculate_from => [qw/ sorted_bam_prefix /],
            calculate => q( return $sorted_bam_prefix.'.bam'; ),
            doc => 'The path of the sorted bam.',
        },
    ],
    has_optional_calculated => [
        sorted_bam_prefix => {
            is => 'Text',
            calculate_from => [qw/ unsorted_bam_path /],
            calculate => q(
                my $sorted_bam_prefix = $unsorted_bam_path;
                $sorted_bam_prefix =~ s/\.bam$/.sorted/;
                return $sorted_bam_prefix;
            ),
            doc => 'The prefix of the sorted bam.',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->debug_message('Sort bams...');

    my $sort_ok = $self->_sort_bam;
    return if not $sort_ok;

    my $verify_read_count_ok = $self->_verify_read_count;
    return if not $verify_read_count_ok;

    my $cleanup_ok = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->remove_paths_and_auxiliary_files($self->unsorted_bam_path);
    return if not $cleanup_ok;

    $self->debug_message('Sort bams...done');
    return 1;
}

sub _sort_bam {
    my $self = shift;

    my $unsorted_bam_path = $self->unsorted_bam_path;
    $self->debug_message("Unsorted bam path: $unsorted_bam_path");

    my $sorted_bam_prefix = $self->sorted_bam_prefix;
    $self->debug_message("Sorted bam prefix: $sorted_bam_prefix");

    my $sorted_bam_path = $self->sorted_bam_path;
    $self->debug_message("Sorted bam path: $sorted_bam_path");

    my $cmd = "samtools sort -m 3000000000 -n $unsorted_bam_path $sorted_bam_prefix";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $sorted_bam_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run samtools sort!');
        return;
    }

    return 1;
}

sub _verify_read_count {
    my $self = shift;
    $self->debug_message('Verify read count...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $unsorted_flagstat = $helpers->load_or_run_flagstat($self->unsorted_bam_path);
    return if not $unsorted_flagstat;

    my $sorted_flagstat = $helpers->load_or_run_flagstat($self->sorted_bam_path);
    return if not $sorted_flagstat;

    $self->debug_message('Sorted bam read count: '.$sorted_flagstat->{total_reads});
    $self->debug_message('Unsorted bam read count: '.$unsorted_flagstat->{total_reads});

    if ( $sorted_flagstat->{total_reads} != $unsorted_flagstat->{total_reads} ) {
        $self->error_message('Sorted and unsorted bam read counts do not match!');
        return;
    }

    $self->debug_message('Verify read count...done');
    return 1;
}

1;

