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
    has_optional_input => {
        instrument_data_and_new_source_files => {
            is => 'Text',
            is_many => 1,
            doc => 'Mapping of instrument data ids and new source files. ',
        },
    },
    has_output => {
        file => {
            is => 'File',
            default_value => '-',
            doc => 'File name [source files tsv] to generate. Use this in the import manager. Defaults to STDOUT.',
        },
    },
    has_optional_transient => {
        _instrument_data_and_new_source_files => { is => 'Hash', },
    },
};

sub help_detail {
    return <<HELP;
HELP
}

sub __errors__ { 
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    my %instrument_data_and_new_source_files;
    my ($instdata_id, $source_file, $previous_insdata_id);
    for my $instrument_data_and_new_source_file ( $self->instrument_data_and_new_source_files ) {
        if ( $instrument_data_and_new_source_file =~ /=/ ) {
            ($previous_insdata_id, $source_file) = split(/=/, $instrument_data_and_new_source_file);
        }
        elsif ( defined $previous_insdata_id ) {
            $source_file = $instrument_data_and_new_source_file;
        }
        else {
            return (
                UR::Object::Tag->create(
                    type => 'invalid',
                    properties => [qw/ instrument_data_and_new_source_files /],
                    desc => 'Mal-formed instrument data and new source files! '.$instrument_data_and_new_source_file,
                )
            );
        }

        if ( not -s $source_file ) {
            return (
                UR::Object::Tag->create(
                    type => 'invalid',
                    properties => [qw/ instrument_data_and_new_source_files /],
                    desc => 'Source file does not exist! '.$source_file,
                )
            );

        }

        push @{$instrument_data_and_new_source_files{$previous_insdata_id}}, $source_file;
    }
    $self->_instrument_data_and_new_source_files(\%instrument_data_and_new_source_files);

    return @errors;
}

sub execute {
    my $self = shift;
    $self->status_message('Generate file to reimport instrument data...');

    my $instrument_data_and_new_source_files = $self->_instrument_data_and_new_source_files;
    my @reimports;
    for my $instrument_data ( $self->instrument_data ) {
        my $reimport = Genome::InstrumentData::Reimport->attributes_for_reimport_from_instrument_data($instrument_data);
        return if not $reimport;

        if ( $instrument_data_and_new_source_files->{ $instrument_data->id } ) {
            $reimport->{source_files} = join(',', @{$instrument_data_and_new_source_files->{ $instrument_data->id }});
        }
        elsif ( not $reimport->{source_files} or not -s $reimport->{source_files} ) {
            $self->error_message("No existing source file for instrument data! %s", $instrument_data->id);
            return;
        }

        push @reimports, $reimport;
    }
    $self->status_message('Found '.@reimports.' instrument data...');

    my @headers = Genome::InstrumentData::Reimport->headers_for_reimport_attributes(@reimports);
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

