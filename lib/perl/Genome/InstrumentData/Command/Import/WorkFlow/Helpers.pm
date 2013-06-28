package Genome::InstrumentData::Command::Import::WorkFlow::Helpers;

use strict;
use warnings;

use Genome;

require Carp;
require File::Copy;

class Genome::InstrumentData::Command::Import::WorkFlow::Helpers { 
    is => 'UR::Singleton',
};

#<INST DATA INFO>#
sub local_source_files_for_instrument_data {
    my ($self, $instrument_data) = @_;

    Carp::confess('No instrument data to get local source files!') if not $instrument_data;

    my $directory = $instrument_data->data_directory;
    Carp::confess('No instrument data directory to get local source files!') if not $directory;

    my @local_source_files;
    for my $source_file ( split(',', $instrument_data->original_data_path) ) {
        my $source_file_basename = File::Basename::basename($source_file);
        $source_file_basename =~ s/\.gz$//;
        push @local_source_files, $directory.'/'.$source_file_basename;
    }

    return @local_source_files;
}
#<>#

#<WORKFLOW>#
sub add_operation_to_workflow {
    my ($self, $workflow, $name) = @_;

    my $command_class_name = 'Genome::InstrumentData::Command::Import::WorkFlow::'.join('', map { ucfirst } split(' ', $name));
    my $operation_type = Workflow::OperationType::Command->create(command_class_name => $command_class_name);
    if ( not $operation_type ) {
        $self->error_message("Failed to create work flow operation for $name");
        return;
    }

    my $operation = $workflow->add_operation(
        name => $name,
        operation_type => $operation_type,
    );

    return $operation;
}
#<>#

#<MOVE>#
sub move_file {
    my ($self, $from, $to) = @_;

    $self->status_message('Move file...');
    my $from_sz = -s $from;
    $self->status_message("From: $from");
    $self->status_message("To: $to");
    my $move_ok = File::Copy::move($from, $to);
    if ( not $move_ok ) {
        $self->error_message('Move failed!');
        return;
    }
    my $to_sz = -s $to;
    if ( not $to_sz or $to_sz != $from_sz ) {
        $self->error_message("Move succeeded, but destination size is diffeerent from original! $to_sz vs $from_sz");
        return;
    }

    $self->status_message('Move file...done');
    return 1;
}
#<>#

#<SAMTOOLS>#
sub run_flagstat {
    my ($self, $bam_path, $flagstat_path) = @_;
    $self->status_message('Run flagstat...');

    $flagstat_path ||= $bam_path.'.flagstat';
    $self->status_message("Bam path: $bam_path");
    $self->status_message("Flagstat path: $flagstat_path");
    my $cmd = "samtools flagstat $bam_path > $flagstat_path";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $flagstat_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run flagstat!');
        return;
    }

    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_path);
    $self->status_message('Flagstat output:');
    $self->status_message( join("\n", map { ' '.$_.': '.$flagstat->{$_} } sort keys %$flagstat) );
    if ( not $flagstat->{total_reads} > 0 ) {
        $self->error_message('Flagstat determined that there are no reads in bam! '.$bam_path);
        return;
    }

    $self->status_message('Run flagstat...done');
    return $flagstat;
}
#<>#

1;

