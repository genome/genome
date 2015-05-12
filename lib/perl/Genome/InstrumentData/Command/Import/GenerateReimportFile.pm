package Genome::InstrumentData::Command::Import::GenerateReimportFile;

use strict;
use warnings;

use Genome;

use Text::CSV;

class Genome::InstrumentData::Command::Import::GenerateReimportFile { 
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
            doc => 'Downsample ratios to add to each instrument data reimport.',
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
        sep_char => { is => 'Text', },
    },
};

sub help_brief {
    return 'generate a file to reimport/import downsampled instrument data';
}

sub help_detail {
    return <<HELP;
Given existing instrument data, this command will generate a file that can be used in the import manager to reimport these instrument data. Reimports are done because the original file(s) imported are unable to be used in pipelines or to import downsampled versions (via given downsample ratios). Additionally, a mapping of instrument data ids and new source files can be given if the imported ones are invalid.
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
        my $instdata_attrs = $self->attributes_for_reimport_from_instrument_data($instrument_data);
        return if not $instdata_attrs;

        my @source_files = ( $instrument_data_and_new_source_files->{ $instrument_data->id } )
        ? @{$instrument_data_and_new_source_files->{ $instrument_data->id }}
        : $instdata_attrs->{'instdata.source_files'};

        for my $source_file ( @source_files ) {
            if ( not $source_file or not -s $source_file ) {
                die $self->error_message("Source file does not exist! $source_file");
            }

            if ( @downsample_ratios ) {
                for my $downsample_ratio ( @downsample_ratios ) {
                    push @reimports, {
                        %$instdata_attrs,
                        'instdata.source_files' => $source_file,
                        'instdata.downsample_ratio' => $downsample_ratio,
                    };
                }
            }
            else {
                push @reimports, {
                    %$instdata_attrs,
                    'instdata.source_files' => $source_file,
                };
            }
        }
    }
    $self->status_message('Found '.@reimports.' instrument data...');
    
    my $file = $self->file;
    my $sep_char = Genome::InstrumentData::Command::Import::CsvParser->resovle_sep_vhar_from_file_extension($file);
    my $writer = Text::CSV->new({
            sep_char => $sep_char,
            empty_is_undef => 1,
            always_quote => ( $sep_char eq ',' ? 1 : 0 ),
        });
    die $self->error_message('Failed to create Text::CSV parser!') if not $writer;

    my $fh = eval{ Genome::Sys->open_file_for_writing($file); };
    if ( not $fh ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to open file! '.$file);
        return 1;
    }

    my %headers = map { $_ => 1 } map { keys %$_ } @reimports;
    for (qw/ library.name instdata.source_files /) { delete $headers{$_}; }
    my @headers = sort keys %headers;
    unshift @headers, (qw/ library.name instdata.source_files /);

    $writer->print($fh, \@headers);
    $fh->print("\n");
    for my $reimport ( @reimports ) {
        $writer->print($fh, [ map { $reimport->{$_} // '' } @headers ]);
        $fh->print("\n");
    }
    $fh->close;
    $self->status_message('Wrote file: '.$file);

    $self->status_message('Success!');
    return 1;
}

sub attribute_labels_to_ignore_when_reimporting {
    (qw/ bam_path base_count genotype_file genotype_file_name
        fragment_count
        ignored import_date import_format is_paired_end
        original_data_path original_data_path_md5
        read_length read_count reference_sequence_build_id
        user_name 
        /);
}

sub attributes_for_reimport_from_instrument_data {
    my ($self, $instrument_data) = @_;

    die 'No instrument data given!' if not $instrument_data;

    my %reimport = ( 
        'library.name' => $instrument_data->library->name,
        'instdata.reimported_from' => $instrument_data->id,
    );

    my $source_file = eval{ $instrument_data->bam_path; };
    if ( not $source_file ) {
        $source_file = eval{ $instrument_data->archive_path; };
    }
    $reimport{'instdata.source_files'} = $source_file if $source_file;

    for my $optional_property_name (qw/ run_name subset_name /) {
        next if not $instrument_data->$optional_property_name;
        $reimport{'instdata.'.$optional_property_name} = $instrument_data->$optional_property_name;
    }

    ATTRIBUTE: for my $attribute ( $instrument_data->attributes ) {
        for my $attribute_label_to_ignore ( $self->attribute_labels_to_ignore_when_reimporting ) {
            next ATTRIBUTE if $attribute->attribute_label eq $attribute_label_to_ignore;
        }
        $reimport{'instdata.'.$attribute->attribute_label} = $attribute->attribute_value;
    }

    return \%reimport;
}

1;

