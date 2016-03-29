package Genome::InstrumentData::Command::Import::WorkFlow::SortBam;

use strict;
use warnings;

require File::Basename;

use Genome;

class Genome::InstrumentData::Command::Import::WorkFlow::SortBam { 
    is => 'Command::V2',
    roles => [qw/
        Genome::InstrumentData::Command::Import::WorkFlow::Role::WithWorkingDirectory
        Genome::InstrumentData::Command::Import::WorkFlow::Role::RemovesInputFiles
    /],
    has_input => [
        bam_path => {
            is => 'FilePath',
            doc => 'The path of the unsorted bam to sort.',
        }
    ],
    has_output => [ 
        output_bam_path => {
            is => 'FilePath',
            calculate_from => [qw/ bam_path /],
            calculate => q| return $self->get_working_bam_path_with_new_extension($bam_path, 'sorted'); |,
            doc => 'The path of the sorted bam.',
        },
    ],
    has_optional_calculated => {
        sorted_bam_prefix => {
            is => 'Text',
            calculate_from => [qw/ output_bam_path /],
            calculate => q|
                $output_bam_path =~ s/\.bam$//;
                return $output_bam_path;
            |,
            doc => 'The prefix of the sorted bam.',
        },
    },
};

sub execute {
    my $self = shift;
    $self->debug_message('Sort bams...');

    my $sort_ok = $self->_sort_bam;
    return if not $sort_ok;

    my $verify_read_count_ok = $self->_verify_read_count;
    return if not $verify_read_count_ok;

    $self->debug_message('Sort bams...done');
    return 1;
}

sub _sort_bam {
    my $self = shift;

    my $bam_path = $self->bam_path;
    $self->debug_message("Unsorted bam path: $bam_path");

    my $sorted_bam_prefix = $self->sorted_bam_prefix;
    $self->debug_message("Sorted bam prefix: $sorted_bam_prefix");

    my $output_bam_path = $self->output_bam_path;
    $self->debug_message("Sorted bam path: $output_bam_path");

    my $cmd = Genome::Model::Tools::Picard::SortSam->execute(
        input_file => $bam_path,
        output_file => $output_bam_path,
        sort_order => 'queryname',
        use_version => '1.82',
    );
    if ( not $cmd->result or not -s $output_bam_path ) {
        $self->error_message('Failed to run picard sort sam!');
        return;
    }

    return 1;
}

sub _verify_read_count {
    my $self = shift;
    $self->debug_message('Verify read count...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;
    my $unsorted_flagstat = $helpers->load_or_run_flagstat($self->bam_path);
    return if not $unsorted_flagstat;

    my $sorted_flagstat = $helpers->load_or_run_flagstat($self->output_bam_path);
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

