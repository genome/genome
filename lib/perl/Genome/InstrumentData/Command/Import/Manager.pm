package Genome::InstrumentData::Command::Import::Manager;

use strict;
use warnings;

use Genome;

use IO::File;
use YAML;

class Genome::InstrumentData::Command::Import::Manager {
    is => 'Command::V2',
    doc => 'Manage importing sequence files into Genome',
    has => [
        source_files_tsv => {
            is => 'Text',
            doc => 'Tab separated file of samples and attributes. See detailed help for more info.',
        },
    ],
    has_optional => [
        import_launch_config => {
            is => 'Text',
            doc => 'Command to use to launch imports. If given, imports will be launched for source files needing import. See detailed help for more info.',
        },
        import_list_config => {
            is => 'Text',
            doc => 'Command and column positions for sample anme and status. Use format: cmd;sample_name_col_#;sub tatus_col_#',
        },
        model_params => {
            is => 'Text',
            is_many => 1,
            doc => 'Parameters to use to create a model for the created instrument data. Give as comma separated key value pairs: param1=value1,param2=value2,... The model class will be derived fromm the processing profile.',
        },
    ],
    has_optional_transient => [
        _import_list_command => { is => 'Text', },
        _import_list_sample_name_column => { is => 'Text', },
        _import_list_status_column => { is => 'Text', },
        _import_command_format => { is => 'Text', },
        _model_class => { is => 'Text' },
        _model_params => { is => 'Hash' },
        _imports => { is => 'Array', },
        instrument_data_import_command_substitutions => { 
            is => 'Hash', 
            default_value => { map { $_ => qr/%{$_}/ } (qw/ sample_name /), },
        },
    ],
};

sub help_detail {
    return <<HELP;
Sample Info File
 This is a tab sparated file in the working directory named "info.tsv" with headers. 
 
 Required Columns
  "name"         The name should be dash (-) separated values of the nomenclature, indivdual id/name and sample id/name.
  "source_files" Source files [bam, fastq, sra, etc] to import for the sample. Separate files by comma (,). 
 
 Additional Columns
  The columns are to specify attributes for sample, patient and instrument data. Prefix the
  attribute name with the entity it is to be assigned to. Some attributes do not need to be
  prefixed, like "gender" and "race", because they are known to go to a certain entity.
 
  sample             "s."
  patient/individual "p."
  instrument_data    "i."
HELP
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    my $import_cmd_error = $self->_resolve_import_command;
    if ( $import_cmd_error ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ import_launch_config /],
            desc => $import_cmd_error,
        );
        return @errors;
    }

    my $import_list_config = $self->import_list_config;
    if ( $import_list_config ) {
        my %import_list_config;
        @import_list_config{qw/ command sample_name_column status_column /} = split(';', $import_list_config);
        for my $attr ( keys %import_list_config ) {
            if ( not defined $import_list_config{$attr} ) {
                push @errors, UR::Object::Tag->create(
                    type => 'invalid',
                    properties => [qw/ import_list_config /],
                    desc => "Missing $attr in $import_list_config",
                );
                return @errors;
            }
            $import_list_config{$attr}-- if $attr =~ /col/;
            my $method = '_import_list_'.$attr;
            $self->$method( $import_list_config{$attr} );
        }
    }

    my $load_info_error = $self->_load_source_files_tsv;
    if ( $load_info_error ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ source_files_tsv /],
            desc => $load_info_error,
        );
        return @errors;
    }

    my $model_params_error = $self->_resolve_model_params;
    if ( $model_params_error ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ model_params /],
            desc => $model_params_error,
        );
        return @errors;
    }

    return;
}

sub _load_source_files_tsv {
    my $self = shift;

    my $source_files_tsv = $self->source_files_tsv;
    if ( not -s $source_files_tsv ) {
        return 'Sample info file does not exist! '.$source_files_tsv;
    }

    my $info_reader = Genome::Utility::IO::SeparatedValueReader->create(
        input => $source_files_tsv,
        separator => "\t",
    );
    if ( not $info_reader ) {
        return 'Failed to open sample info file! '.$source_files_tsv;
    }

    my %headers_not_found = ( sample_name => 1, source_files => 1, );
    for my $header ( @{$info_reader->headers} ) {
        delete $headers_not_found{$header};
    }

    if ( %headers_not_found ) {
        return 'No '.join(' ', map { '"'.$_.'"' } keys %headers_not_found).' column in sample info file! '.$self->source_files_tsv;
    }

    my @imports;
    my %seen_source_files;
    while ( my $hash = $info_reader->next ) {
        my $source_files = delete $hash->{source_files};
        if ( $seen_source_files{$source_files} ) {
            $self->error_message('Duplicate source files! '.$source_files);
            return;
        }
        $seen_source_files{$source_files}++;
        my $sample_name = delete $hash->{sample_name};
        if ( not $sample_name ) {
            $self->error_message('No sample name in info file on line '.$info_reader->line_number.'!');
            return;
        }
        my $import = {
            sample_name => $sample_name,
            source_files => $source_files,
            instrument_data_attributes => [],
        };
        push @imports, $import;
        for my $attr ( sort keys %$hash ) {
            my $value = $hash->{$attr};
            next if not defined $value or $value eq '';
            push @{$import->{'instrument_data_attributes'}}, $attr."='".$value."'";
        }
    }
    $self->_imports(\@imports);

    return;
}

sub _resolve_model_params {
    my $self = shift;

    my @model_params = $self->model_params;
    return if not @model_params;

    my %model_params;
    for my $param ( @model_params ) {
        my ($name, $value) = split('=', $param);
        if ( not defined $value ) {
            return "No value for model param $name!";
        }
        $model_params{$name} = $value;
    }

    if ( not $model_params{processing_profile_id} ) {
        return "No processing profile id for model in params! @model_params";
    }
    my $processing_profile_id = delete $model_params{processing_profile_id};
    $model_params{processing_profile} = Genome::ProcessingProfile->get($processing_profile_id);
    if ( not $model_params{processing_profile} ) { 
        return 'Failed to get processing profile for id! '.$processing_profile_id;
    }

    my @processing_profile_class_parts = split('::', $model_params{processing_profile}->class);
    my $model_class = join('::', $processing_profile_class_parts[0], 'Model', $processing_profile_class_parts[2]);
    my $model_meta = $model_class->__meta__;
    if ( not $model_meta ) {
        return 'Failed to get model class meta! '.$model_class;
    }

    for my $param_name ( keys %model_params ) {
        next if $param_name !~ s/_id$//;
        my $param_value_id = delete $model_params{$param_name.'_id'};
        my $param_property = $model_meta->property_meta_for_name($param_name);
        return "Failed to get '$param_name' property from model class! ".$model_class if not $param_property;
        my $param_value_class = $param_property->data_type;
        my $param_value = $param_value_class->get($param_value_id);
        return "Failed to get $param_value_class for id! $param_value_id" if not $param_value;
        $model_params{$param_name} = $param_value;
    }

    $self->_model_class($model_class);
    $self->_model_params(\%model_params);

    return;
}

sub _resolve_import_command {
    my $self = shift;

    my $cmd_format = $self->import_launch_config;
    if ( $cmd_format ) {
        my $substitutions = $self->instrument_data_import_command_substitutions;
        my $sample_name_substitution = $substitutions->{sample_name};
        if ( $cmd_format !~ /$sample_name_substitution/ ) {
            $self->error_message('No sample name substitutions (%{sample_name} in launch command! '.$cmd_format);
            return;
        }
        $cmd_format .= ' ';
        # TODO rm unecessary subs
    }

    $cmd_format .= 'genome instrument-data import basic --sample name=%{sample_name} --source-files %s --import-source-name %s%s',
    $self->_import_command_format($cmd_format);

    return;
}

sub execute {
    my $self = shift;

    my $load_samples = $self->_load_samples;
    return if not $load_samples;

    my $load_instrument_data = $self->_load_instrument_data;
    return if not $load_instrument_data;

    my $load_models = $self->_load_models;
    return if not $load_models;

    my $load_sample_statuses = $self->_load_sample_statuses;
    return if not $load_sample_statuses;

    my $launch_imports = $self->_launch_imports;
    return if not $launch_imports;

    return $self->_output_status;
}

sub _load_samples {
    my $self = shift;

    my $imports = $self->_imports;
    for my $sample ( @$imports ) {
        $sample->{sample} = Genome::Sample->get(name => $sample->{sample_name});
        next if not $sample->{sample};
        if ( not $sample->{sample}->nomenclature ) {
            $self->error_message('Sample does not have any nomenclature! '.$sample->{sample}->name);
            return;
        }
    }

    $self->_imports($imports);

    return 1;
}

sub _load_instrument_data {
    my $self = shift;

    my $imports = $self->_imports;

    my %instrument_data = map { $_->original_data_path, $_ } Genome::InstrumentData::Imported->get(
        original_data_path => [ map { $_->{source_files} } @$imports ],
        '-hint' => [qw/ attributes /],
    );

    for my $import ( @$imports ) {
        my $instrument_data = $instrument_data{ $import->{source_files} };
        $import->{instrument_data} = $instrument_data;
        $import->{instrument_data_file} = eval{
            my $attribute = $instrument_data->attributes(attribute_label => 'bam_path');
            if ( not $attribute ) {
                $attribute = $instrument_data->attributes(attribute_label => 'archive_path');
            }
            return $attribute->attribute_value if $attribute;
        };

    }

    $self->_imports($imports);

    return 1;
}

sub _load_models {
    my $self = shift;

    my $model_params = $self->_model_params;
    return 1 if not $model_params;

    my $imports = $self->_imports;
    my $model_class = $self->_model_class;
    for my $sample ( @$imports ) {
        # only get/create model once inst data has been created
        next if not $sample->{instrument_data};
        # and has data file
        next if not $sample->{instrument_data_file} or not -s $sample->{instrument_data_file};

        # Only get model for these params and inst data
        $model_params->{subject} = $sample->{sample};
        $model_params->{'instrument_data.id'} = $sample->{instrument_data}->id;

        my $model = $model_class->get(%$model_params);
        if ( not $model ) {
            $model = $model_class->create(%$model_params);
            if ( not $model ) {
                $self->error_message('Failed to create model for sample! '.$sample->{sample_name});
                return;
            }
            $model->add_instrument_data( $sample->{instrument_data} );
        }

        $sample->{model} = $model;
        $sample->{build} = $model->latest_build;
    }

    $self->_imports($imports);

    return 1;
}

sub _load_sample_statuses {
    my $self = shift;

    my $sample_job_statuses = $self->_load_sample_job_statuses;
    return if not $sample_job_statuses;

    my $get_status_for_sample = sub{
        my $sample = shift;
        return 'sample_needed' if not $sample->{sample};
        return 'import_'.$sample->{job_status} if $sample->{job_status};
        return 'import_needed' if not $sample->{instrument_data};
        return 'import_failed' if not defined $sample->{instrument_data_file} or not -s $sample->{instrument_data_file};
        return 'model_needed' if not $sample->{model};
        return 'build_needed' if not $sample->{build};
        return 'build_'.lc($sample->{build}->status);
    };

    my $imports = $self->_imports;
    for my $sample ( @$imports ) {
        $sample->{job_status} = $sample_job_statuses->{ $sample->{sample_name} };
        $sample->{status} = $get_status_for_sample->($sample);
    }

    $self->_imports($imports);

    return 1;
}

sub _load_sample_job_statuses {
    my $self = shift;

    my $job_list_cmd = $self->_import_list_command;
    $job_list_cmd .= ' 2>/dev/null |';
    my $fh = IO::File->new($job_list_cmd);
    if ( not $fh ) {
        $self->error_message('Failed to execute import list command! '.$job_list_cmd);
        return;
    }
    
    my $name_column = $self->_import_list_sample_name_column;
    my $status_column = $self->_import_list_status_column;
    my %sample_job_statuses;
    while ( my $line = $fh->getline ) {
        chomp $line;
        my @tokens = split(/\s+/, $line);
        $sample_job_statuses{ $tokens[$name_column] } = lc $tokens[$status_column];
    }

    return \%sample_job_statuses;
}

sub _launch_imports {
    my $self = shift;

    return 1 if not $self->import_launch_config;

    my $imports = $self->_imports;
    for my $sample ( @$imports ) {
        next if not grep { $sample->{status} ne $_ } (qw/ import_needed import_failed /);

        if ( $sample->{status} eq 'import_failed' and $sample->{instrument_data} ) {
            $sample->{instrument_data}->delete;
            $sample->{instrument_data} = undef;
        }

        my $cmd = $self->_resolve_instrument_data_import_command_for_sample($sample);
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            $self->error_message($@) if $@;
            $self->error_message('Failed to launch instrument data import command for sample! '.$sample->{sample_name});
            return;
        }

        $sample->{status} = 'import_pend';
    }

    return 1;
}

sub _resolve_instrument_data_import_command_for_sample {
    my ($self, $sample) = @_;

    Carp::confess('No samples to resolved instrument data import command!') if not $sample;

    my $cmd_format = $self->_import_command_format;
    my $substitutions = $self->instrument_data_import_command_substitutions;
    for my $name ( keys %$substitutions ) {
        my $value = $sample->{$name};
        if ( not defined $value ) {
            $self->error_message("No value for attribute ($name) specified in launch import command. ".$self->instrument_data_import_command_format);
            return;
        }
        my $pattern = $substitutions->{$name};
        $cmd_format =~ s/$pattern/$value/g;
    }

    my $cmd .= sprintf(
        $cmd_format,
        $sample->{source_files},
        $sample->{sample}->nomenclature,
        ( 
            $sample->{instrument_data_attributes}
            ? ' --instrument-data-properties '.join(',', @{$sample->{instrument_data_attributes}})
            : ''
        ),
    );

    return $cmd;
}

sub _output_status {
    my $self = shift;

    my %totals;
    my $status = join("\t", (qw/ name status inst_data model build /))."\n";
    for my $sample ( sort { $a->{sample_name} cmp $b->{sample_name} } @{$self->_imports} ) {
        $totals{total}++;
        $totals{ $sample->{status} }++;
        $totals{build}++ if $sample->{status} =~ /^build/;
        my ($model_id, $build_id) = (qw/ NA NA /);
        if ( $sample->{model} ) {
            $model_id = $sample->{model}->id;
            if ( $sample->{build} and $sample->{status} ne 'build_requested' ) {
                $build_id = $sample->{build}->id;
            }
        }
        $status .= join(
            "\t",
            $sample->{sample_name},
            $sample->{status}, 
            ( $sample->{instrument_data} ? $sample->{instrument_data}->id : 'NA' ),
            $model_id,
            $build_id,
        )."\n";
    }

    print STDOUT "$status";
    print STDERR "\nSummary:\n".join("\n", map { sprintf('%-16s %s', $_, $totals{$_}) } sort { $a cmp $b } keys %totals)."\n";

    return 1;
}

1;

