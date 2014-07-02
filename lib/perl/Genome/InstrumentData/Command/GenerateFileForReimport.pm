package Genome::InstrumentData::Command::GenerateFileForReimport;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Reimport;

class Genome::InstrumentData::Command::GenerateFileForReimport { 
    is => 'Command::V2',
    has_input => {
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'Instrument data to use to create file to reimport.',
        },
    },
    has_output => {
        file => {
            is => 'File',
            doc => 'File name [source files tsv] to generate. Use this in the import manager.',
        },
    },
};

sub help_detail {
    return <<HELP;
HELP
}

sub execute {
    my $self = shift;
    $self->status_message('Generate file to reimport instrument data...');

    my @reimports;
    for my $instrument_data ( $self->instrument_data ) {
        my $reimport = Genome::InstrumentData::Reimport->attributes_for_reimport($instrument_data);
        return if not $reimport;

        if ( not $reimport->{source_files} or not -s $reimport->{source_files} ) {
            $self->error_message("No source file for instrument data! %s", $instrument_data->id);
            return;
        }

        push @reimports, $reimport;
    }
    $self->status_message('Found '.@reimports.' instrument data...');

    my @headers = Genome::InstrumentData::Reimport->headers_for_attributes_for_reimports(@reimports);
    return if not @headers;

    my $file = $self->file;
    my $fh = eval{ Genome::Sys->open_file_for_writing($file); };
    if ( not $fh ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to open file! '.$file);
        return 1;
    }

    $fh->print( join("\t", @headers)."\n" );
    for my $reimport ( @reimports ) {
        $fh->print( join("\t", map { $reimport->{$_} // '' } @headers)."\n" );
    }
    $fh->close;

    $self->status_message('Wrote file: '.$file);
    $self->status_message('Success!');
    return 1;
}

1;

