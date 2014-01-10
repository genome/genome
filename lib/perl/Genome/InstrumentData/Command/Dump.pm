package Genome::InstrumentData::Command::Dump;

use strict;
use warnings;

use Genome;
use Cwd;

class Genome::InstrumentData::Command::Dump {
    is => 'Command::V2',
    has => [
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => 'Instrument data to dump.'
        },
        directory => {
            is => 'Text',
            doc => 'Directory to dump to. Defaults to current working directory.',
            default_value => Cwd::getcwd(),
        },
    ],
    doc => 'Dump files from instrument data',
};

sub help_detail {
    return 'Dump fastqs from instrument data';
}

sub execute {
    my $self = shift;

    my $instrument_data = $self->instrument_data;
    if ( not $instrument_data ) {
        $self->error_message('No instrument data given!');
        return;
    }
    $self->debug_message('Instrument data: '.join(' ', map { $instrument_data->$_ } (qw/ id sequencing_platform /)));

    if ( not $instrument_data->can('dump_sanger_fastq_files') ) {
        $self->error_message('Instrument data can not dump fastq files!');
        return;
    }

    unless (Genome::Sys->validate_directory_for_write_access($self->directory)) {
        $self->error_message('Failed to validate directory '. $self->directory .' for write access!');
        return;
    }
    $self->debug_message('Directory: '.$self->directory);

    my @files = $instrument_data->dump_sanger_fastq_files(directory=>$self->directory);
    if ( not @files ) {
        $self->error_message('Failed to dump files!');
        return;
    }

    $self->debug_message('Finished!');
    return 1;
}

1;

