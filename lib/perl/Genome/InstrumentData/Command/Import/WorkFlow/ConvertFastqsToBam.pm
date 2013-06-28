package Genome::InstrumentData::Command::Import::WorkFlow::ConvertFastqsToBam;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require List::Util;
require List::MoreUtils;

class Genome::InstrumentData::Command::Import::WorkFlow::ConvertFastqsToBam { 
    is => 'Command::V2',
    has_input => [
        instrument_data => { 
            is => 'Genome::InstrumentData',
            is_output => 1,
            doc => 'Instrument data.',
        },
    ],
    has_output => [ 
        bam_path => {
            calculate_from => [qw/ instrument_data /],
            calculate => q( return $instrument_data->data_directory.'/tmp/unsorted.bam'; ),
            doc => 'The path of the bam.',
        },
    ],
};

sub execute {
    my $self = shift;
    $self->status_message('Convert fastqs to bam...');

    my $convert_ok = $self->_convert_fastqs_to_bam;
    return if not $convert_ok;

    #my $flagstat = Genome::InstrumentData::Command::Helpers->run_flagstat();

    $self->status_message('Convert fastqs to bam...done');
    return 1;
}

sub _convert_fastqs_to_bam {
    my $self = shift;
    $self->status_message('Run picard fastq to sam...');

    my $instrument_data = $self->instrument_data;
    my @fastqs = Genome::InstrumentData::Command::Import::WorkFlow::Helpers->local_source_files_for_instrument_data($instrument_data);
    if ( not @fastqs ) {
        $self->error_message('No local source files for instrument data! '.$instrument_data->id);
        return;
    }

    $self->status_message("Fastq 1: $fastqs[0]");
    my $bam_path = $self->bam_path;
    my $cmd = "gmt picard fastq-to-sam --fastq $fastqs[0] --output $bam_path --quality-format Standard --sample-name ".$instrument_data->sample->name.' --read-group-name '.$instrument_data->library->name;
    if ( $fastqs[1] ) {
        $self->status_message("Fastq 2: $fastqs[1]") if $fastqs[1];
        $cmd .= ' --fastq2 '.$fastqs[1]
    }
    $self->status_message("Bam path: $bam_path");

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    if ( not $rv or not -s $bam_path ) {
        $self->error_message($@) if $@;
        $self->error_message('Failed to run picard fastq to sam!');
        return;
    }

    $self->status_message('Run picard fastq to sam...done');
    return 1;
}

1;

