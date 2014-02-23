package Genome::InstrumentData::Command::CalculateReadCount;

use strict;
use warnings;

use Carp "confess";

use Genome;

class Genome::InstrumentData::Command::CalculateReadCount {
    is => 'Command::V2',
    has => {
        instrument_data => {
            is => 'Genome::InstrumentData',
            doc => "Instrument data that has 'bam_path' and will accept the 'read_count' attribute.",
            shell_args_position => 1,
        },
    },
};

sub help_brief {
    return 'A command to calculate the number of reads in the bam file associated with this instrument-data.';
}

sub help_detail {
    return <<EOS
    This tool will use samtools flagstat to determine the number of reads in a bam file.
    The number of reads will then be stored as an attribute with attribute_label="read_count" on the
    instrument-data.
EOS
}

sub execute {
    my $self = shift;

    my $instrument_data = $self->instrument_data;
    # does the instrument data already have this done?
    my @read_count_attrs = $instrument_data->attributes(attribute_label=>"read_count"); 
    if (scalar (@read_count_attrs)) {
        $self->status_message("Found existing read_count attribute was already defined.");
        return 1;
    }

    my $id = $instrument_data->id;
    # find the bam on the instrument_data object
    my $bam_path = $instrument_data->bam_path;
    confess "Couldn't find bam for instrument_data: $id\n" unless $bam_path;
    $self->debug_message("Found bam_path: $bam_path");

    # check to see if a flagstat file already exists
    my $flagstat_file = $bam_path . ".flagstat";
    if( not -e $flagstat_file ){
        $flagstat_file = Genome::Sys->create_temp_file_path();
        # get the flagstat object
        my $flagstat_object = Genome::Model::Tools::Sam::Flagstat->create(bam_file => $bam_path, output_file =>$flagstat_file);
        $self->debug_message("Created flagstat object to create output file: $flagstat_file");
        # determine the read-count of the bam
        unless($flagstat_object->execute()){
            confess "Failed to run flagstat on bam file at: $bam_path";
        }
        $self->debug_message("Ran flagstat and about to parse output file...");
    } else {
        $self->debug_message("Found existing flagstat and about to parse output file...");
    }

    my $flag_data = Genome::Model::Tools::Sam::Flagstat->parse_file_into_hashref($flagstat_file);
    my $read_count = $flag_data->{"total_reads"};
    $self->debug_message("Found $read_count for the bam file at: $bam_path");

    # create the read_count attribute on the instrument_data
    my $new_attribute = $instrument_data->add_attribute(
        attribute_label => 'read_count',
        attribute_value => $read_count, 
    );
    if (not $new_attribute){
        confess "Failed to create new read_count attribute on instrument: $id";
    } else {
        $self->debug_message("Everything complete.");
        return 1;
    }
}

1;
