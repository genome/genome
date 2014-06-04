package Genome::InstrumentData::Command::Import::GenerateSourceFilesTsvForReimport;

use strict;
use warnings;

use Genome;

class Genome::InstrumentData::Command::Import::GenerateSourceFilesTsvForReimport { 
    is => 'Command::V2',
    has_input => {
        instrument_data => {
            is => 'Genome::InstrumentData',
            is_many => 1,
            doc => 'Instrument data to use to create source files TSV to reimport.',
        },
        #skip_errors => {},
    },
    has_output => {
        source_files_tsv => {
            is => 'File',
            doc => 'Source files tsv to generate.',
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
    $self->status_message('Generate source files tsv to reimport instrument data...');

    # Gather the params for each inst data
    my @reimports;
    my %inst_data_attrs = ( derived_from => 1, ); # 1 is placeholder value
    for my $instrument_data ( $self->instrument_data ) {

        my $source_file = $instrument_data->archive_path;
        if ( not $source_file or not -s $source_file ) {
            $self->error_message("No source file for instrument data! %s", $instrument_data->id);
            return;
        }

        my %reimport = ( 
            derived_from => $instrument_data->id,
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
    # TODO move ignored attribute labels to a common class
    for my $ignored_attr_label (qw/ bam_path genotype_file genotype_file_name ignored import_date import_format user_name /) {
        delete $inst_data_attrs{$ignored_attr_label};
    }

    my @headers = sort keys %inst_data_attrs;
    unshift @headers, (qw/ library_name source_files /);


    # Write the source files tsv
    my $source_files_tsv = $self->source_files_tsv;
    my $source_files_tsv_fh = eval{ Genome::Sys->open_file_for_writing($source_files_tsv); };
    if ( not $source_files_tsv_fh ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to open source files tsv! '.$source_files_tsv);
        return 1;
    }

    $source_files_tsv_fh->print( join("\t", @headers)."\n" );
    for my $reimport ( @reimports ) {
        $source_files_tsv_fh->print( join("\t", map { $reimport->{$_} // '' } @headers)."\n" );
    }
    $source_files_tsv_fh->close;

    $self->status_message('Wrote source files tsv: '.$source_files_tsv);
    $self->status_message('Success!');
    return 1;
}

1;

