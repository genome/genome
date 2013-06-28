package Genome::InstrumentData::Command::Import::WorkFlow::Helpers;

use strict;
use warnings;

use Genome;

require Carp;

class Genome::InstrumentData::Command::Import::WorkFlow::Helpers { 
    is => 'UR::Singleton',
};

#<INST DATA INFO>#
sub local_file_for_source_file {
    my ($self, $source_file) = @_;

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

sub run_flagstat {
    my ($self, $bam_file, $flagstat_file) = @_;
    $self->status_message('Run flagstat...');

    $flagstat_file ||= $bam_file.'.flagstat';
    $self->status_message("Flagstat file: $flagstat_file");
    my $cmd = "samtools flagstat $bam_file > $flagstat_file";
    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $flagstat_file ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run flagstat!');
        return;
    }

    my $flagstat = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_file);
    $self->status_message('Flagstat output:');
    $self->status_message( join("\n", map { ' '.$_.': '.$flagstat->{$_} } sort keys %$flagstat) );
    if ( not $flagstat->{total_reads} > 0 ) {
        $self->error_message('Flagstat determined that there are no reads in bam! '.$bam_file);
        return;
    }

    $self->status_message('Run flagstat...done');
    return $flagstat;
}

1;

