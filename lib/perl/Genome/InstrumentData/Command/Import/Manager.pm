package Genome::InstrumentData::Command::Import::Manager;

use strict;
use warnings;

use Genome;

use Genome::InstrumentData::Command::Import::WorkFlow::SourceFiles;
use Genome::InstrumentData::Command::Import::CsvParser;

use IO::File;
use Params::Validate ':types';

class Genome::InstrumentData::Command::Import::Manager {
    is => 'Command::V2',
    doc => 'Manage importing sequence files into GMS',
    has => [
        analysis_project => {
            is => 'Genome::Config::AnalysisProject',
            doc => 'Analysis project to assign to the created instrument data.',
        },
        file => {
            is => 'Text',
            doc => Genome::InstrumentData::Command::Import::CsvParser->csv_help,
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
Output

 WHAT            SENT TO  DESCRIPTION    
 Status          STDOUT   One for each library/source files set.
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

sub execute {
    my $self = shift;

    $self->_resolve_launch_command;
    $self->_resolve_list_config;
    $self->_load_file;
    $self->_check_source_files_and_set_kb_required_for_processing;
    $self->_load_instrument_data;
    $self->_load_statuses;
    $self->_launch_imports;
    $self->_output_status;

    return 1
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

    return 1;
}

sub _resolve_list_config {
    my $self = shift;

    my $list_config = $self->list_config;
    return 1 if not $list_config;

    my %list_config;
    @list_config{qw/ command job_name_column status_column /} = split(';', $list_config);
    for my $attr ( sort keys %list_config ) {
        if ( not defined $list_config{$attr} ) {
            die $self->error_message("Missing %s in list config: %s", $attr, $list_config);
        }
        $list_config{$attr}-- if $attr =~ /col/;
        my $method = '_list_'.$attr;
        $self->$method( $list_config{$attr} );
    }

    return 1;
}

sub _load_file {
    my $self = shift;

    my $parser = Genome::InstrumentData::Command::Import::CsvParser->create(file => $self->file);
    my (@imports, %seen);
    while ( my $import = $parser->next ) {
        my $string = join(' ', $import->{library}->{name}, map { $import->{instdata}->{$_} } keys %{$import->{instdata}});
        $string .= $import->{instdata}->{downsample_ratio} if $import->{instdata}->{downsample_ratio};
        my $id = substr(Genome::Sys->md5sum_data($string), 0, 6);
        if ( $seen{$id} ) {
            die $self->error_message("Duplicate source file/library combination! $string");
        }
        $seen{$id}++;
        $import->{job_name} = $id;
        $import->{library_name} = $import->{library}->{name};
        my @libraries = Genome::Library->get(name => $import->{library}->{name});
        $import->{library_cnt} = scalar @libraries;
        push @imports, $import;
    }
    $self->_imports(\@imports);

    return 1;
}

sub _check_source_files_and_set_kb_required_for_processing {
    my $self = shift;

    for my $import ( @{$self->_imports} ) {
        # get disk space required [checks if source files exist]
        my $disk_space_required_in_kb = Genome::InstrumentData::Command::Import::WorkFlow::SourceFiles->create(
            paths => [ split(',', $import->{instdata}->{source_files}) ],
        )->kilobytes_required_for_processing;
        $disk_space_required_in_kb = 1048576 if $disk_space_required_in_kb < 1048576; # 1 Gb 
        $import->{gtmp} = sprintf('%.0f', $disk_space_required_in_kb / 1048576);
        $import->{mtmp} = sprintf('%.0f', $disk_space_required_in_kb / 1024);
        $import->{kbtmp} = $disk_space_required_in_kb;
    }

    return 1;
}

sub _load_instrument_data {
    my $self = shift;

    # Get instdata by source files
    my @instrument_data = Genome::InstrumentData::Imported->get(
        original_data_path => [ map { $_->{instdata}->{source_files} } @{$self->_imports} ],
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

    for my $import ( @{$self->_imports} ) {
        my $lookup_id = $import->{instdata}->{source_files};
        $lookup_id .= sprintf('%f', $import->{instdata}->{downsample_ratio}) if $import->{instdata}->{downsample_ratio};
        $import->{instrument_data} = $instrument_data{$lookup_id};
    }

    return 1;
}

sub _load_statuses {
    my $self = shift;

    my $job_statuses = $self->_load_job_statuses;
    return if not $job_statuses;

    my $get_status_for_import = sub{
        my $import = shift;
        return 'no_library' if $import->{library_cnt} == 0;
        return 'too_many_libraries' if $import->{library_cnt} > 1;
        return $import->{job_status} if $import->{job_status};
        return 'needed' if not $import->{instrument_data};
        return 'success';
    };

    my $imports = $self->_imports;
    for my $import ( @$imports ) {
        $import->{job_status} = $job_statuses->{ $import->{job_name} };
        $import->{status} = $get_status_for_import->($import);
    }

    return 1;
}

sub _load_job_statuses {
    my $self = shift;

    my $job_list_cmd = $self->_list_command;
    $job_list_cmd .= ' 2>/dev/null |';
    my $fh = IO::File->new($job_list_cmd);
    die $self->error_message('Failed to execute import list command! '.$job_list_cmd) if not $fh;
    
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

    my $launch_sub = sub{
        my $import = shift;
        my $cmd = $self->_resolve_launch_command_for_import($import);
        my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
        if ( not $rv ) {
            $self->error_message($@) if $@;
            die $self->error_message('Failed to launch instrument data import command!');
        }
        $import->{status} = 'pend';
    };

    my $imports = $self->_imports;
    for my $import ( @$imports ) {
        next if $import->{status} ne 'needed';
        $launch_sub->($import) or return;
    }

    return 1;
}

sub _resolve_launch_command_for_import {
    my ($self, $import) = Params::Validate::validate_pos(@_, {type => OBJECT}, {type => HASHREF});

    my $cmd_format = $self->_launch_command_format;
    my $substitutions = $self->_launch_command_substitutions;
    for my $name ( keys %$substitutions ) {
        my $value = $import->{$name};
        next if not defined $value;
        my $pattern = $substitutions->{$name};
        $cmd_format =~ s/$pattern/$value/g;
    }

    my %instrument_data_properties = %{$import->{instdata}};
    my $source_files = delete $instrument_data_properties{source_files};
    my $downsample_ratio = delete $instrument_data_properties{downsample_ratio};
    my $cmd .= sprintf(
        $cmd_format,
        $source_files,
        ( $import->{sample}->{nomenclature} ),
        ( 
            %instrument_data_properties
            ? ' --instrument-data-properties '.  join(',', map { $_."='".$instrument_data_properties{$_}."'"; } sort keys %instrument_data_properties)
            : ''
        ),
        $self->analysis_project ? " --analysis-project id=".$self->analysis_project->id : '',
        ( 
            defined $downsample_ratio
            ? " --downsample-ratio ".$downsample_ratio
            : '' 
        ),
    );

    return $cmd;
}

sub _output_status {
    my $self = shift;

    my @status = ( ['library_name'], ['job_name'], ['status'], ['inst_data'], );
    my ($i, @row, %totals);
    for my $import ( sort { $a->{library}->{name} cmp $b->{library}->{name} } @{$self->_imports} ) {
        $totals{total}++;
        $totals{ $import->{status} }++;
        @row = (
            $import->{library}->{name}, 
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

