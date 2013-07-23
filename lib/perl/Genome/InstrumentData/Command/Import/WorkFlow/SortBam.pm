package Genome::InstrumentData::Command::Import::WorkFlow::SortBam;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::WorkFlow::SortBam { 
    is => 'Command::V2',
    has_input => [
        unsorted_bam_paths => {
            is => 'Text',
            is_many => 1,
            doc => 'The path of the unsorted bam to sort.',
        }
    ],
    has_output => [ 
        sorted_bam_paths => {
            is => 'Text',
            is_many => 1,
            doc => 'The path of the sorted bam.',
        },
    ],
};

sub sorted_bam_prefix_for_unsorted_bam_path {
    my ($self, $unsorted_bam_path) = @_;
    Carp::confess('No unsorted bam path to derive sorted bam prefix!') if not $unsorted_bam_path;
    my $sorted_bam_prefix = $unsorted_bam_path;
    $sorted_bam_prefix =~ s/\.bam//;
    $sorted_bam_prefix .= '.sorted';
    return $sorted_bam_prefix;
}

sub execute {
    my $self = shift;
    $self->status_message('Sort bams...');

    my @unsorted_bam_paths = $self->unsorted_bam_paths;
    $self->status_message('Unsorted bam count: '.@unsorted_bam_paths);
    my @sorted_bam_paths;
    my $cnt = 0;
    for my $unsorted_bam_path ( @unsorted_bam_paths ) {
        $self->status_message('Bam #'.++$cnt);

        my $sorted_bam_path = $self->_sort_bam($unsorted_bam_path);
        return if not $sorted_bam_path;

        my $verify_read_count_ok = $self->_verify_read_count($unsorted_bam_path, $sorted_bam_path);
        return if not $verify_read_count_ok;

        push @sorted_bam_paths, $sorted_bam_path;
    }
    $self->sorted_bam_paths(\@sorted_bam_paths);

    $self->status_message('Sort bams...done');
    return 1;
}

sub _sort_bam {
    my ($self, $unsorted_bam_path) = @_;

    $self->status_message("Unsorted bam path: $unsorted_bam_path");
    my $sorted_bam_prefix = $self->sorted_bam_prefix_for_unsorted_bam_path($unsorted_bam_path);
    $self->status_message("Sorted bam prefix: $sorted_bam_prefix");
    my $sorted_bam_path = $sorted_bam_prefix.'.bam';
    $self->status_message("Sorted bam path: $sorted_bam_path");
    my $cmd = "samtools sort -m 3000000000 -n $unsorted_bam_path $sorted_bam_prefix";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $sorted_bam_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run samtools sort!');
        return;
    }

    return $sorted_bam_path;
}

sub _verify_read_count {
    my ($self, $unsorted_bam_path, $sorted_bam_path) = @_;
    $self->status_message('Verify read count...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $unsorted_flagstat = $helpers->load_flagstat($unsorted_bam_path.'.flagstat');
    return if not $unsorted_flagstat;

    my $sorted_flagstat = $helpers->validate_bam($sorted_bam_path);
    return if not $sorted_flagstat;

    $self->status_message('Sorted bam read count: '.$sorted_flagstat->{total_reads});
    $self->status_message('Unsorted bam read count: '.$unsorted_flagstat->{total_reads});

    if ( $sorted_flagstat->{total_reads} != $unsorted_flagstat->{total_reads} ) {
        $self->error_message('Sorted and unsorted bam read counts do not match!');
        return;
    }

    $self->status_message('Verify read count...done');
    return 1;
}

1;

