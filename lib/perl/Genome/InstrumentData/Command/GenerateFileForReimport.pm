package Genome::InstrumentData::Command::GenerateFileForReimport;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Imported::Reimport;

class Genome::InstrumentData::Command::GenerateFileForReimport { 
    is => 'Command::V2',
    has_input => {
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'Instrument data to use to create file to reimport.',
        },
        #skip_errors => {},
    },
    has_output => {
        file => {
            is => 'File',
            doc => 'File name [source files tsv] to generate. Use this in the import manager.',
        },
    },
    has_optional_transient => {
    },
};

sub help_detail {
    return <<HELP;
HELP
}

sub execute {
    my $self = shift;
    $self->status_message('Generate file to reimport instrument data...');

    # Gather the params for each inst data
    my @reimports;
    my $reimported_attribute_label = Genome::InstrumentData::Reimport->attribute_label_for_reimported_from;
    my %headers = ( $reimported_attribute_label => 1 );
    for my $instrument_data ( $self->instrument_data ) {

        my $source_file = $instrument_data->archive_path;
        if ( not $source_file or not -s $source_file ) {
            $self->error_message("No source file for instrument data! %s", $instrument_data->id);
            return;
        }

        my %reimport = ( 
            $reimported_attribute_label => $instrument_data->id,
            library_name => $instrument_data->library->name,
            source_files => $source_file,
        );
        for my $optional_property_name (qw/ run_name subset_name /) {
            next if not $instrument_data->$optional_property_name;
            $reimport{$optional_property_name} = $instrument_data->$optional_property_name;
            $inst_data_attrs{$optional_property_name}++;
        }
        push @reimports, \%reimport;

        for my $attribute ( $instrument_data->attributes ) {
            my $attribute_value = 
            $reimport{ $attribute->attribute_label } = $attribute->attribute_value;
            $inst_data_attrs{ $attribute->attribute_label }++;
        }
    }
    $self->status_message('Found '.@reimports.' instrument data...');

    # Determine the headers
    for my $ignored_attr_label ( Genome::InstrumentData::Reimport->attribute_labels_to_ignore_when_reimporting ) {
        delete $inst_data_attrs{$ignored_attr_label};
    }

    my @headers = sort keys %inst_data_attrs;
    unshift @headers, (qw/ library_name source_files /);

    # Write the file
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

