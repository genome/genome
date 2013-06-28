package Genome::InstrumentData::Command::Import::WorkFlow::SortBam;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::WorkFlow::SortBam { 
    is => 'Command::V2',
    has_input => [
        instrument_data => { 
            is => 'Genome::InstrumentData',
            is_output => 1,
            doc => 'Instrument data.',
        },
        unsorted_bam_path => {
            is => 'Text',
            doc => 'The path of the unsorted bam to sort.',
        }
    ],
    has_output => [ 
        sorted_bam_path => {
            calculate_from => [qw/ instrument_data /],
            calculate => q( return $instrument_data->data_directory.'/tmp/sorted.bam'; ),
            doc => 'The path of the sortred bam.',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Sort bam...');

    my $sort_ok = $self->_sort_bam;
    return if not $sort_ok;

    #my $flagstat = Genome::InstrumentData::Command::Helpers->run_flagstat();

    $self->status_message('Sort bam...done');
    return 1;
}

sub _sort_bam {
    my $self = shift;

    my $unsorted_bam_path = $self->unsorted_bam_path;
    $self->status_message("Unsorted bam path: $unsorted_bam_path");
    my $sorted_bam_path = $self->sorted_bam_path;
    my $sorted_bam_prefix = $sorted_bam_path;
    $sorted_bam_prefix =~ s/\.bam$//;
    $self->status_message("Sorted bam prefix: $sorted_bam_prefix");
    $self->status_message("Sorted bam path: $sorted_bam_path");
    my $cmd = "samtools sort -m 3000000000 -n $unsorted_bam_path $sorted_bam_prefix";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $sorted_bam_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run samtools sort!');
        return;
    }

    return 1;
}

1;

