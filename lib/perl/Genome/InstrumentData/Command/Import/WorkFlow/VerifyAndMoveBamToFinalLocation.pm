package Genome::InstrumentData::Command::Import::WorkFlow::VerifyAndMoveBamToFinalLocation;

use strict;
use warnings;

use Genome;

require File::Path;

class Genome::InstrumentData::Command::Import::WorkFlow::VerifyAndMoveBamToFinalLocation { 
    is => 'Command::V2',
    has_input => [
        instrument_data => { 
            is => 'Genome::InstrumentData',
            is_output => 1,
            doc => 'Instrument data.',
        },
        bam_path => {
            is => 'Text',
            doc => 'The path of the bam to verify and move.',
        }
    ],
    has_output => [ 
        final_bam_path => {
            calculate_from => [qw/ instrument_data /],
            calculate => q( return $instrument_data->data_directory.'/all_sequences.bam'; ),
            doc => 'The path of the final bam.',
        },
    ],
    has_calculated => [
        flagstat_path => {
            calculate_from => [qw/ final_bam_path /],
            calculate => q( return $final_bam_path.'.flagstat'; ),
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Verify and move bam to permanent location...');

    my $helpers = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->get;

    my $bam_path = $self->bam_path;
    my $flagstat_path = $self->flagstat_path;
    my $flagstat = $helpers->validate_bam($bam_path, $flagstat_path);
    return if not $flagstat;

    my $final_bam_path = $self->final_bam_path;
    $self->status_message("Final bam path: $final_bam_path");
    my $move_ok = $helpers->move_file($bam_path, $final_bam_path);
    return if not $move_ok;

    my $instrument_data = $self->instrument_data;
    $instrument_data->add_attribute(attribute_label => 'read_count', attribute_value => $flagstat->{total_reads});
    $instrument_data->add_attribute(attribute_label => 'bam_path', attribute_value => $final_bam_path);

    my $tmp_dir = $instrument_data->data_directory.'/tmp';
    if ( -d $tmp_dir ) {
        $self->status_message('Remove tmp dir...');
        File::Path::rmtree($tmp_dir, 1);
        $self->status_message('Remove tmp dir...done');
    }

    $self->status_message('Reallocate...');
    $instrument_data->allocations->reallocate;# with move??
    $self->status_message('Reallocate...done');

    $self->status_message('Verify and move bam to permanent location...done');
    return 1;
}

1;

