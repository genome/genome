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
 The columns are to specify attributes for instrument data. Attributes are skipped if they are empty for a particular instrument data.

Example

This will look to run/check imports for 2 samples. Sample 1 has 2 source files [bams] to import. The second sample, Sample-02, has 2 fastqs to import. They will be downloaded, unzipped and converted to bam. In addition, the flow_cell_id, lane and index_sequence will be added to the respective instrument data as attributes.

sample_name source_files    flow_cell_id    lane    index_sequence
Sample-01   sample-1.1.bam  XAXAXA 1   AATTGG
Sample-01   sample-1.2.bam  XAXAXA 1   TTAACC
Sample-02   http://fastqs.org/sample-2.fwd.fastq.gz,http://fastqs.org/sample-2.rev.fastq.gz  XYYYYX  2   GAACTT
DOC

        },
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'Analysis project to assign to the created instrument data.',
        },
    ],
    has_optional => [
        launch_config => {
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
        list_config => {
            is => 'Text',
            doc => <<DOC
The command to run to list running imports. Give the command, job name column number and status column number, separated by semicolons (;). The command will be run and parsed, expecting to find the job name and status. This is then used to determine the next course of action for each source file(s).

Example for LSF
 This will run the bjobs command, looking for the job name in column 7 and status in column 3.

 bjobs -w -g /me/mygroup;7;3

DOC
        },
        show_import_commands => {
            is => 'Boolean',
            doc => 'Show the import commands for source files that need to be imported *instead* of executing them.',
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
    ],
};

sub help_detail {
    return <<HELP;
Outputs
There are 3 potential outputs:

 WHAT            SENT TO  DESCRIPTION    
 Status          STDOUT   One for each sample/source files set.
 Import command  STDERR   Command to import the source files. Printed if the --show-import-commands option is indicated.
 Stats           STDERR   Summary stats for statuses.

Status Output
 The status output is formatted so that the columns are lined up. A line is output for each sample name and source file set.

 Example:

 sample_name # status  inst_data
 Sample-01   1 success 6c43929f065943a09f8ccc769e42c41d
 Sample-01   2 run     NA
 Sample-02   1 needed  NA

 Column definitions:

 sample_name  Name of the sample.
 #            Since the sample name may be used more than once, this is the iteration as it appears in the source files tsv.
 status       Import status - no_sample, no_library, needed, pend, run, success.
 inst_data    Instrument data id.

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

    $cmd_format .= 'genome instrument-data import basic --sample name=%{sample_name} --analysis-project id=%s --source-files %s --import-source-name %s%s',
    $self->_launch_command_format($cmd_format);

    return;
}

sub execute {
    my $self = shift;

    my $load_samples = $self->_load_samples;
    return if not $load_samples;

    my $load_instrument_data = $self->_load_instrument_data;
    return if not $load_instrument_data;

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
        # sample name, number and job name
        my $sample_name = $import->{sample_name};
        $import->{sample_number} = ++$sample_names_seen{$sample_name};
        $import->{job_name} = $import->{sample_name}.'.'.$import->{sample_number};
        # genome sample
        $import->{sample} = Genome::Sample->get(name => $sample_name);
        next if not $import->{sample};
        # nomenclature
        $import->{sample}->{nomenclature} = $import->{sample}->nomenclature // 'WUGC'; #FIXME
        # library
        my @libraries = Genome::Library->get( # TODO base on library
            name => $sample_name.'-extlibs',
            sample => $import->{sample},
        );
        $import->{libraries} = \@libraries if @libraries; 
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

sub _load_sample_statuses {
    my $self = shift;

    my $sample_job_statuses = $self->_load_sample_job_statuses;
    return if not $sample_job_statuses;

    my $get_status_for_import = sub{
        my $import = shift;
        return 'no_sample' if not $import->{sample};
        return 'no_library' if not $import->{libraries};
        return 'too_many_libraries' if @{$import->{libraries}} > 1;
        return $import->{job_status} if $import->{job_status};
        return 'needed' if not $import->{instrument_data};
        return 'needed' if not defined $import->{instrument_data_file} or not -s $import->{instrument_data_file};
        return 'success';
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

    if ( not $self->show_import_commands ) {
        # not printing commands, can they be launched?
        if ( not $self->launch_config ) {
            $self->warning_message('Cannot launch jobs because there is no launch config!');
            return 1;
        }
        elsif ( not $self->list_config ) { 
            $self->warning_message('Can not launch jobs because there is no list config!');
            return 1;
        }
        elsif ( not $self->_launch_command_has_job_name ) {
            $self->warning_message('Cannot launch jobs because there is no %{job_name} in launch config!');
            return 1;
        }
    }

    my $launch_sub = sub{
        my $import = shift;
        my $cmd = $self->_resolve_launch_command_for_import($import);
        if ( $self->show_import_commands ) {
            print STDERR "$cmd\n";
        }
        else {
            my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
            if ( not $rv ) {
                $self->error_message($@) if $@;
                $self->error_message('Failed to launch instrument data import command!');
                return;
            }
            $import->{status} = 'pend';
        }
    };

    my $imports = $self->_imports;
    for my $import ( @$imports ) {
        next if $import->{status} ne 'needed';
        $launch_sub->($import) or return;
    }
    print STDERR "\n";

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
        $self->analysis_project->id,
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

    my @status = ( ['sample_name'], ['#'], ['status'], ['inst_data'], );
    my ($i, @row, %totals);
    for my $import ( sort { $a->{sample_name} cmp $b->{sample_name} } @{$self->_imports} ) {
        $totals{total}++;
        $totals{ $import->{status} }++;
        @row = (
            $import->{sample_name},
            $import->{sample_number}, 
            $import->{status}, 
            ( $import->{instrument_data} ? $import->{instrument_data}->id : 'NA' ),
        );
        for ( $i = 0; $i < @row; $i++ ) {
            push @{$status[$i]}, $row[$i];
        }
    }

    my @column_formats;
    for ( $i = 0; $i < @status; $i++ ) {
        my ($column_width) = sort { $b <=> $a } map { length($_) } @{$status[$i]};
        $column_formats[$i] = '%-'.$column_width.'s';
    }

    my $format = join(' ', @column_formats)."\n";
    my $status;
    my $rownum = @{$status[0]};
    for ( $i = 0; $i < $rownum; $i++ ) {
        $status .= sprintf($format, map { $_->[$i] } @status);
    }

    printf STDOUT $status;
    print STDERR "\nSummary:\n".join("\n", map { sprintf('%-16s %s', $_, $totals{$_}) } sort { $a cmp $b } keys %totals)."\n";

    return 1;
}

1;

