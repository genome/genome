package Genome::InstrumentData::Command::Import::Manager;

use strict;
use warnings;

use Genome;

use IO::File;

class Genome::InstrumentData::Command::Import::Manager {
    is => 'Command::V2',
    doc => 'Manage importing sequence files into GMS',
    has => [
        source_files_tsv => {
            is => 'Text',
            doc => <<DOC
TAB separated file containing sample names, source files and instrument data attributes to import.
 
Required Columns
 sample_name     Name of the sample. The sample and extlibs library must exist.
 source_files    Source files [bam, fastq, sra, etc] to import. Separate files by comma (,). 

Additional Columns
 The columns are to specify attributes for instrument data. Attributes are skipped if they are empty for a particular instruemnt data.

Example

This will look to run/check imports for 2 samples. Sample 1 has 2 source files [bams] to import. The second sample, Sample-02, has 2 fastqs to import. They will be downloaded, unzipped and converted to bam. In addition, the flow_cell_id, lane and index_sequence will be added to the respective instrument data as attributes.

sample_name source_files    flow_cell_id    lane    index_sequence
Sample-01   sample-1.1.bam  XAXAXA 1   AATTGG
Sample-01   sample-1.2.bam  XAXAXA 1   TTAACC
Sample-02   http://fastqs.org/sample-2.fwd.fastq.gz,http://fastqs.org/sample-2.rev.fastq.gz  XYYYYX  2   GAACTT
DOC

        },
    ],
    has_optional => [
        launch_config => {
            is => 'Text',
            doc => <<DOC
The command to run to list running imports. Give the command, job name column number and status column number, separated by semicolons (;). The command will be run and parsed, expecting to find the job name and status. This is then used to determine the next course of action for each source file(s).

Example for LSF
 This will run the bjobs command, looking for the job name in column 7 and status in column 3.

 bjobs -w -g /me/mygroup;7;3
            
DOC
        },
        list_config => {
            is => 'Text',
            doc => <<DOC
Launch imports [if needed] using this command. Insert '%{job_name}' into the command so the manager can monitor status.
 
The import commands will be printed to the screen [on STDERR] if:
 launch config is not given
 launch config does not have a '%{job_name}' in it
 list config is not given

Example for LSF
 Launch the job into group /me/mygroup and logging to /users/me/logs/%{job_name}

 bsub -J %{job_name} -g /me/mygroup -oo /users/me/logs/%{job_name}

DOC
        },
        model_params => {
            is => 'Text',
            is_many => 1,
            doc => 'Parameters to use to create a model for the created instrument data. Give as comma separated key value pairs: param1=value1,param2=value2,... The model class will be derived fromm the processing profile.',
        },
    ],
    has_optional_transient => [
        _imports => { is => 'Array', },
        _list_command => { is => 'Text', },
        _list_job_name_column => { is => 'Text', },
        _list_status_column => { is => 'Text', },
        _launch_command_format => { is => 'Text', },
        _launch_command_has_job_name => { is => 'Boolean', default_value => 0, },
        _launch_command_substitutions => { 
            is => 'Hash', 
            default_value => { map { $_ => qr/%{$_}/ } (qw/ job_name sample_name /), },
        },
        _model_class => { is => 'Text' },
        _model_params => { is => 'Hash' },
    ],
};

sub help_detail {
    return <<HELP;
Manage importing sequence files into GMS.
HELP
}

sub __errors__ {
    my $self = shift;

    my @errors = $self->SUPER::__errors__;
    return @errors if @errors;

    my $import_cmd_error = $self->_resolve_launch_command;
    if ( $import_cmd_error ) {
        push @errors, UR::Object::Tag->create(
            type => 'invalid',
            properties => [qw/ launch_config /],
            desc => $import_cmd_error,
        );
        return @errors;
    }

    my $list_config = $self->list_config;
    if ( $list_config ) {
        my %list_config;
        @list_config{qw/ command job_name_column status_column /} = split(';', $list_config);
        for my $attr ( keys %list_config ) {
            if ( not defined $list_config{$attr} ) {
                push @errors, UR::Object::Tag->create(
                    type => 'invalid',
                    properties => [qw/ list_config /],
                    desc => "Missing $attr in $list_config",
                );
                return @errors;
            }
            $list_config{$attr}-- if $attr =~ /col/;
            my $method = '_list_'.$attr;
            $self->$method( $list_config{$attr} );
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
    if ( not -f $source_files_tsv or not -s $source_files_tsv ) {
        return 'Invalid source files tsv! '.$source_files_tsv;
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
            push @{$import->{instrument_data_attributes}}, $attr."='".$value."'";
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

sub _resolve_launch_command {
    my $self = shift;

    my $cmd_format = $self->launch_config;
    if ( $cmd_format ) {
        my $substitutions = $self->_launch_command_substitutions;
        my $required_substitution = 'job_name';
        my $substitution = $substitutions->{$required_substitution};
        if ( $cmd_format =~ /$substitution/ ) {
            $self->_launch_command_has_job_name(1);
        }
        $cmd_format .= ' ';
    }

    $cmd_format .= 'genome instrument-data import basic --sample name=%{sample_name} --source-files %s --import-source-name %s%s',
    $self->_launch_command_format($cmd_format);

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
    my %sample_names_seen;
    for my $import ( @$imports ) {
        $import->{sample} = Genome::Sample->get(name => $import->{sample_name});
        next if not $import->{sample};
        $import->{sample}->{nomenclature} = $import->{sample}->nomenclature // 'WUGC'; #FIXME
        my $sample_name = $import->{sample}->name;
        $sample_names_seen{$sample_name}++;
        $import->{job_name} = $sample_name;
        $import->{job_name} .= '.'.$sample_names_seen{$sample_name} if $sample_names_seen{$sample_name} > 1;
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
    for my $import ( @$imports ) {
        # only get/create model once inst data has been created
        next if not $import->{instrument_data};
        # and has data file
        next if not $import->{instrument_data_file} or not -s $import->{instrument_data_file};

        # Only get model for these params and inst data
        $model_params->{subject} = $import->{sample};
        $model_params->{'instrument_data.id'} = $import->{instrument_data}->id;

        my $model = $model_class->get(%$model_params);
        if ( not $model ) {
            $model = $model_class->create(%$model_params);
            if ( not $model ) {
                $self->error_message('Failed to create model for sample! '.$import->{sample_name});
                return;
            }
            $model->add_instrument_data( $import->{instrument_data} );
        }

        $import->{model} = $model;
        $import->{build} = $model->latest_build;
    }

    $self->_imports($imports);

    return 1;
}

sub _load_sample_statuses {
    my $self = shift;

    my $sample_job_statuses = $self->_load_sample_job_statuses;
    return if not $sample_job_statuses;

    my $get_status_for_import = sub{
        my $import = shift;
        return 'sample_needed' if not $import->{sample};
        return 'import_'.$import->{job_status} if $import->{job_status};
        return 'import_needed' if not $import->{instrument_data};
        return 'import_needed' if not defined $import->{instrument_data_file} or not -s $import->{instrument_data_file};
        return 'model_needed' if not $import->{model};
        return 'build_needed' if not $import->{build};
        return 'build_'.lc($import->{build}->status);
    };

    my $imports = $self->_imports;
    for my $import ( @$imports ) {
        $import->{job_status} = $sample_job_statuses->{ $import->{job_name} } if $import->{job_name};
        $import->{status} = $get_status_for_import->($import);
    }

    $self->_imports($imports);

    return 1;
}

sub _load_sample_job_statuses {
    my $self = shift;

    my $job_list_cmd = $self->_list_command;
    $job_list_cmd .= ' 2>/dev/null |';
    my $fh = IO::File->new($job_list_cmd);
    if ( not $fh ) {
        $self->error_message('Failed to execute import list command! '.$job_list_cmd);
        return;
    }
    
    my $name_column = $self->_list_job_name_column;
    my $status_column = $self->_list_status_column;
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

    my $launch_sub;
    if ( not $self->list_config or not $self->_launch_command_has_job_name or not $self->launch_config ) {
        $launch_sub = sub{
            my $import = shift;
            my $cmd = $self->_resolve_launch_command_for_import($import);
            print STDOUT "$cmd\n";
        };
    }
    else {
        $launch_sub = sub{
            my $import = shift;
            my $cmd = $self->_resolve_launch_command_for_import($import);
            my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
            if ( not $rv ) {
                $self->error_message($@) if $@;
                $self->error_message('Failed to launch instrument data import command!');
                return;
            }
            $import->{status} = 'import_pend';
        }
    }

    my $imports = $self->_imports;
    for my $import ( @$imports ) {
        next if $import->{status} ne 'import_needed';
        $launch_sub->($import) or return;
    }

    return 1;
}

sub _resolve_launch_command_for_import {
    my ($self, $import) = @_;

    Carp::confess('No import to resolve launch command!') if not $import;

    my $cmd_format = $self->_launch_command_format;
    my $substitutions = $self->_launch_command_substitutions;
    for my $name ( keys %$substitutions ) {
        my $value = $import->{$name};
        next if not defined $value;
        my $pattern = $substitutions->{$name};
        $cmd_format =~ s/$pattern/$value/g;
    }

    my $cmd .= sprintf(
        $cmd_format,
        $import->{source_files},
        $import->{sample}->{nomenclature},
        ( 
            @{$import->{instrument_data_attributes}}
            ? ' --instrument-data-properties '.join(',', @{$import->{instrument_data_attributes}})
            : ''
        ),
    );

    return $cmd;
}

sub _output_status {
    my $self = shift;

    my %totals;
    my $status = join("\t", (qw/ name status inst_data model build /))."\n";
    for my $import ( sort { $a->{sample_name} cmp $b->{sample_name} } @{$self->_imports} ) {
        $totals{total}++;
        $totals{ $import->{status} }++;
        $totals{build}++ if $import->{status} =~ /^build/;
        my ($model_id, $build_id) = (qw/ NA NA /);
        if ( $import->{model} ) {
            $model_id = $import->{model}->id;
            if ( $import->{build} and $import->{status} ne 'build_requested' ) {
                $build_id = $import->{build}->id;
            }
        }
        $status .= join(
            "\t",
            $import->{sample_name},
            $import->{status}, 
            ( $import->{instrument_data} ? $import->{instrument_data}->id : 'NA' ),
            $model_id,
            $build_id,
        )."\n";
    }

    print STDERR "$status";
    print STDERR "\nSummary:\n".join("\n", map { sprintf('%-16s %s', $_, $totals{$_}) } sort { $a cmp $b } keys %totals)."\n";

    return 1;
}

1;

