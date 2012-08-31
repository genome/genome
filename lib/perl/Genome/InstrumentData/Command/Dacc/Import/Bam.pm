package Genome::InstrumentData::Command::Dacc::Import::Bam;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
require File::Basename;
require File::Path;

class Genome::InstrumentData::Command::Dacc::Import::Bam {
    is  => 'Genome::InstrumentData::Command::Dacc::Import',
    has => [
        format => {
            is => 'Text',
            is_constant => 1,
            value => 'bam',
        },
    ],
};

sub help_brief {
    return 'Import the dl\'d BAMs from the DACC';
}

sub help_detail {
    return help_brief();
}

sub _execute {
    my $self = shift;

    $self->status_message('Move BAM...');
    
    my $dl_directory = $self->_dl_directory;
    my @bams = glob($dl_directory.'/*.bam');
    if ( @bams > 1 ) { 
        $self->error_message("Expected only one BAM in download directory ($dl_directory), but found ".scalar(@bams));
        return;
    }

    my $instrument_data = $self->_main_instrument_data;
    my $destination_file = $instrument_data->destination_file;
    if ( @bams ) {
        $instrument_data->original_data_path($bams[0]);
        if ( not UR::Context->commit ) { # must save this now
            $self->error_message('Cannot commit instrument data original path.');
            return;
        }
        my $move = $self->_move_file($bams[0], $destination_file);
        return if not $move;
    }
    if ( not -e $destination_file ){
        $self->error_message("BAM not found in download directory: $dl_directory");
        return;
    }

    $self->status_message('Move BAM...OK');

    $self->status_message('Update instrument data...');

    my %metrics = $self->_load_bam_metric_file;
    return if not %metrics;
    if ( $self->sra_sample_id ne $metrics{id} ) {
        $self->error_message('SRA ID does not match the one in the BAM metrics file: '.$self->sra_sample_id.' v. '.$metrics{id});
        return;
    }
    if ( not defined $metrics{reads} ) {
        $self->error_message('No read count in BAM metric file');
        return;
    }
    $instrument_data->read_count($metrics{reads});
    $instrument_data->description($self->sra_sample_id.' BAM of mapped reads from the DACC');

    $self->status_message('Update instrument data...OK');

    return 1;
}

sub _load_bam_metric_file {
    my $self = shift;

    $self->status_message('Load BAM metrics...');

    my ($metric_file) = glob($self->_dl_directory.'/*_metric.txt');
    if ( not $metric_file ) {
        $self->error_message('Cannot find metric file in dl directory: '.$self->_dl_directory);
        return;
    }

    my $fh = eval{ Genome::Sys->open_file_for_reading($metric_file); };
    if ( not $fh ) {
        $self->error_message("Cannot open metric file ($metric_file): $@");
        return;
    }

    my $line = $fh->getline;
    if ( not $line ) {
       $self->error_message("No data in metric file: $metric_file");
       return;
    }
    chomp $line;

    my %metrics;
    @metrics{qw/ id reads clusters mapped /} = split(/\s+/, $line);

    $self->status_message('Load BAM metrics...OK');

    return %metrics;

}

1;

