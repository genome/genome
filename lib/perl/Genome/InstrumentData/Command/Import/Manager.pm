package Genome::InstrumentData::Command::Import::Manager;

use strict;
use warnings;

use Genome;

use IO::File;

class Genome::InstrumentData::Command::Import::Manager {
    is => 'Command::V2',
    doc => 'Manage importing sequence files into GMS',
    has => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'Analysis project to assign to the created instrument data.',
        },
        source_files_tsv => {
            is => 'Text',
            doc => <<DOC
TAB separated file containing library names, source files and instrument data attributes to import.
 
Required Columns
 library_name    Name of the library. The library must exist.
 source_files    Source files [bam, fastq, sra, etc] to import. Separate files by comma (,). 

Additional Columns
 The columns are to specify attributes for instrument data. Attributes are skipped if they are empty for a particular instrument data.

Example

This will look to run/check imports for 2 libraires. Library 1 has 2 source files [bams] to import. The second library, Sample-02-extlibs, has 2 fastqs to import. They will be downloaded, unzipped and converted to bam. In addition, the flow_cell_id, lane and index_sequence will be added to the respective instrument data as attributes.

library_name        source_files    flow_cell_id    lane    index_sequence
Sample-01-extlibs   sample-1.1.bam  XAXAXA          1       AATTGG
Sample-01-extlibs   sample-1.2.bam  XAXAXA          1       TTAACC
Sample-02-extlibs   http://fastqs.org/sample-2.fwd.fastq.gz,http://fastqs.org/sample-2.rev.fastq.gz  XYYYYX  2   GAACTT

DOC
        },
    ],
    has_optional => [
        launch_config => {
            is => 'Text',
            doc => <<DOC
Launch imports [if needed] using this command. Insert '%{job_name}' place holder into the command so the manager can monitor status.
 
The import commands will be printed to the screen [on STDERR] if:
 launch config is not given
 launch config does not have a '%{job_name}' in it
 list config is not given

The temp space required can also be filled in be using the tmp family of placeholders. 
 %{gtmp}  gigabytes
 %{mtmp}  megabytes
 %{kbtmp} kilobytes

Example for LSF
 Launch the job into group /me/mygroup and logging to /users/me/logs/%{job_name}

 bsub -J %{job_name} -g /me/mygroup -oo /users/me/logs/%{job_name} -M 16000000 -R 'select [mem>16000 & gtmp>%{gtmp}] rsuage[mem=16000,gtmp=%{gtmp}]'

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
            default_value => { map { $_ => qr/%{$_}/ } (qw/ job_name library_name gtmp mtmp kbtmp /), },
        },
    ],
};

sub help_detail {
    return <<HELP;
Outputs
There are 3 potential outputs:

 WHAT            SENT TO  DESCRIPTION    
 Status          STDOUT   One for each library/source files set.
 Import command  STDERR   Command to import the source files. Printed if the --show-import-commands option is indicated.
 Stats           STDERR   Summary stats for statuses.

Status Output
 The status output is formatted so that the columns are lined up. A line is output for each library name and source file set.

 Example:

 library_name      job_name status  inst_data
 Sample-01-extlibs ahsgdj   success 6c43929f065943a09f8ccc769e42c41d
 Sample-01-extlibs jfuexm   run     NA
 Sample-02-extlibs oodkmq   needed  NA

 Column definitions:

 library_name Name of the library.
 #            Since the library name may be used more than once, this is the iteration as it appears in the source files tsv.
 status       Import status - no_library, needed, pend, run, success.
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
        return 'Failed to open source files tsv! '.$source_files_tsv;
    }

    my %headers_not_found = ( library_name => 1, source_files => 1, );
    for my $header ( @{$info_reader->headers} ) {
        delete $headers_not_found{$header};
    }

    if ( %headers_not_found ) {
        return 'No '.join(' ', map { '"'.$_.'"' } keys %headers_not_found).' column in source files tsv! '.$self->source_files_tsv;
    }

    my (@imports, %seen);
    while ( my $hash = $info_reader->next ) {
        my $id = substr(Genome::Sys->md5sum_data( join('', map { $hash->{$_} } sort keys %$hash) ), 0, 6);
        my $source_files = delete $hash->{source_files};
        if ( $seen{$id} ) {
            $self->error_message('Duplicate entry on line '.$info_reader->line_number.'!');
            return;
        }
        $seen{$id}++;
        my $library_name = delete $hash->{library_name};
        if ( not $library_name ) {
            $self->error_message('No library name in source files tsv on line '.$info_reader->line_number.'!');
            return;
        }
        my $import = {
            library_name => $library_name,
            source_files => $source_files,
            job_name => $id,
        };
        push @imports, $import;
        my %instrument_data_properties;
        for my $name ( sort keys %$hash ) {
            my $value = $hash->{$name};
            next if not defined $value or $value eq '';
            if ( defined $instrument_data_properties{$name} ) {
                $self->error_message('');
                return;
            }
            $instrument_data_properties{$name} = $value;
        }
        $import->{instrument_data_properties} = \%instrument_data_properties;
    }

    $self->_imports(\@imports);

    return $info_reader->error_message;

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

    $cmd_format .= "genome instrument-data import basic --library name=%{library_name} --source-files %s --import-source-name '%s'%s%s%s",
    $self->_launch_command_format($cmd_format);

    return;
}

sub execute {
    my $self = shift;

    my $source_files_ok = $self->_check_source_files_and_set_kb_required_for_processing;
    return if not $source_files_ok;
 
    my $load_libraries = $self->_load_libraries;
    return if not $load_libraries;

    my $load_instrument_data = $self->_load_instrument_data;
    return if not $load_instrument_data;

    my $load_statuses = $self->_load_statuses;
    return if not $load_statuses;

    my $launch_imports = $self->_launch_imports;
    return if not $launch_imports;

    return $self->_output_status;
}

sub _check_source_files_and_set_kb_required_for_processing {
    my $self = shift;

    my $imports = $self->_imports;
    my %library_names_seen;
    for my $import ( @$imports ) {
        # get disk space required [checks if source files exist]
        my $disk_space_required_in_kb = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->kilobytes_required_for_processing_of_source_files( split(',', $import->{source_files}) );
        return if Genome::InstrumentData::Command::Import::WorkFlow::Helpers->error_message;
        $disk_space_required_in_kb = 1048576 if $disk_space_required_in_kb < 1048576; # 1 Gb 
        $import->{gtmp} = sprintf('%.0f', $disk_space_required_in_kb / 1048576);
        $import->{mtmp} = sprintf('%.0f', $disk_space_required_in_kb / 1024);
        $import->{kbtmp} = $disk_space_required_in_kb;
     }
    $self->_imports($imports);

    return if $self->error_message;
    return 1;
}

sub _load_libraries {
    my $self = shift;

    my $imports = $self->_imports;
    my %library_names_seen;
    for my $import ( @$imports ) {
        # library name, number and job name
        my $library_name = $import->{library_name};
        $import->{library_number} = ++$library_names_seen{$library_name};
        # genome library - get as array in case there are many with the same name
        my @libraries = Genome::Library->get(name => $library_name);
        $import->{library} = \@libraries if @libraries;
        next if not $import->{library}; # error will be displayed later
    }
    $self->_imports($imports);

    return if $self->error_message;
    return 1;
}

sub _load_instrument_data {
    my $self = shift;

    my $imports = $self->_imports;

    # Get instdata by source files
    my @instrument_data = Genome::InstrumentData::Imported->get(
        original_data_path => [ map { $_->{source_files} } @$imports ],
        '-hint' => [qw/ attributes /],
    );

    # Create map w/ original data path files and downsample ratio as ids
    my %instrument_data;
    for my $instrument_data ( @instrument_data ) {
        my $original_data_path = $instrument_data->original_data_path;
        my $downsample_ratio_attr = $instrument_data->attributes(attribute_label => 'downsample_ratio');
        my $id = $original_data_path;
        $id .= sprintf('%f', $downsample_ratio_attr->attribute_value) if $downsample_ratio_attr;
        push @{$instrument_data{$id}}, $instrument_data;
    }

    for my $import ( @$imports ) {
        my $lookup_id = $import->{source_files};
        $lookup_id .= sprintf('%f', $import->{instrument_data_properties}->{downsample_ratio}) if $import->{instrument_data_properties}->{downsample_ratio};
        $import->{instrument_data} = $instrument_data{$lookup_id};
    }

    $self->_imports($imports);

    return 1;
}

sub _load_statuses {
    my $self = shift;

    my $job_statuses = $self->_load_job_statuses;
    return if not $job_statuses;

    my $get_status_for_import = sub{
        my $import = shift;
        return 'no_library' if not $import->{library};
        return 'too_many_libraries' if @{$import->{library}} > 1;
        return $import->{job_status} if $import->{job_status};
        return 'needed' if not $import->{instrument_data};
        return 'success';
    };

    my $imports = $self->_imports;
    for my $import ( @$imports ) {
        $import->{job_status} = $job_statuses->{ $import->{job_name} };
        $import->{status} = $get_status_for_import->($import);
    }

    $self->_imports($imports);

    return 1;
}

sub _load_job_statuses {
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
    my %job_statuses;
    while ( my $line = $fh->getline ) {
        chomp $line;
        my @tokens = split(/\s+/, $line);
        $job_statuses{ $tokens[$name_column] } = lc $tokens[$status_column];
    }

    return \%job_statuses;
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

    my $instrument_data_properties = $import->{instrument_data_properties};
    my $cmd .= sprintf(
        $cmd_format,
        $import->{source_files},
        ( $import->{library}->[0]->sample->nomenclature // 'WUGC' ), #FIXME nomenclature
        ( 
            %{$import->{instrument_data_properties}}
            ? ' --instrument-data-properties '.  join(',', map { $_."='".$instrument_data_properties->{$_}."'"; } sort keys %$instrument_data_properties)
            : ''
        ),
        $self->analysis_project ? " --analysis-project id=".$self->analysis_project->id : '',
        ( 
            defined $instrument_data_properties->{downsample_ratio}
            ? " --downsample-ratio ".$instrument_data_properties->{downsample_ratio} 
            : '' 
        ),
    );

    return $cmd;
}

sub _output_status {
    my $self = shift;

    my @status = ( ['library_name'], ['job_name'], ['status'], ['inst_data'], );
    my ($i, @row, %totals);
    for my $import ( sort { $a->{library_name} cmp $b->{library_name} } @{$self->_imports} ) {
        $totals{total}++;
        $totals{ $import->{status} }++;
        @row = (
            $import->{library_name}, 
            $import->{job_name}, 
            $import->{status}, 
            ( $import->{instrument_data} ? join(' ', map { $_->id } @{$import->{instrument_data}}) : 'NA' ),
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

