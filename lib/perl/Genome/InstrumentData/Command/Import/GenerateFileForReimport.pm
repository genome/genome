package Genome::InstrumentData::Command::Import::GenerateFileForReimport;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Reimport;

class Genome::InstrumentData::Command::Import::GenerateFileForReimport { 
    is => 'Command::V2',
    has_input => {
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'Instrument data to use to create file to reimport.',
        },
    },
    has_many_optional_input => {
        downsample_ratios => {
            is => 'Float',
            doc => 'Downsample ratios to add to each instruemnt data reimport.',
        },
        instrument_data_and_new_source_files => {
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

        push @{$instrument_data_and_new_source_files{$previous_insdata_id}}, $source_file;
    }
    $self->_instrument_data_and_new_source_files(\%instrument_data_and_new_source_files);

    return @errors;
}

sub execute {
    my $self = shift;
    $self->status_message('Generate file to reimport instrument data...');

    my @downsample_ratios = $self->downsample_ratios;
    my $instrument_data_and_new_source_files = $self->_instrument_data_and_new_source_files;
    my @reimports;
    for my $instrument_data ( $self->instrument_data ) {
        my $instdata_attrs = Genome::InstrumentData::Reimport->attributes_for_reimport_from_instrument_data($instrument_data);
        return if not $instdata_attrs;

        my @source_files = ( $instrument_data_and_new_source_files->{ $instrument_data->id } )
        ? @{$instrument_data_and_new_source_files->{ $instrument_data->id }}
        : $instdata_attrs->{source_files};

        for my $source_file ( @source_files ) {
            if ( not $source_file or not -s $source_file ) {
                die $self->error_message("Source file does not exist! $source_file");
            }

            if ( @downsample_ratios ) {
                for my $downsample_ratio ( @downsample_ratios ) {
                    push @reimports, {
                        %$instdata_attrs,
                        source_files => $source_file,
                        downsample_ratio => $downsample_ratio,
                    };
                }
            }
            else {
                push @reimports, {
                    %$instdata_attrs,
                    source_files => $source_file,
                };
            }
        }
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

