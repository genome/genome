package Genome::Model::Tools::CgHub::Query;

use strict;
use warnings;

use Genome;

class Genome::Model::Tools::CgHub::Query {
    is => "Genome::Model::Tools::CgHub::Base",
    has_input => {
        uuid => {
            is => "Text",
            is_output => 1,
        },
    },
    has_optional_output => {
        xml_file => {
            is => "Text",
            doc => 'Save XML output to this file instead of formatted results.',
        },
        output_file => {
            doc => 'Save standard output/error to this file. If given the xml file, the will only contain the status of the query.',
        },
    },
    doc => 'Query CG Hub for metadata.',
};

sub _build_command {
    my $self = shift;

    my $cmd = "cgquery";
    my $xml_file = $self->xml_file;
    $cmd .= " -o $xml_file" if $xml_file;
    $cmd .= " analysis_id=".$self->uuid;

    return $cmd;
}

sub _verify_success {
    my $self = shift;

    my $output_file = $self->output_file;
    if ( not -s $output_file ) {
        die $self->error_message('CG Query failed to generate output file!');
    }

    my $output_fh = Genome::Sys->open_file_for_reading($output_file);
    my $number_of_hits = 0;
    while ( my $line = $output_fh->getline ) {
        if ( $line =~ /Matching Objects\s+:\s+(\d+)/ ) {
            $number_of_hits = $1;
            last;
        }
    }
    $output_fh->close;

    if ( $number_of_hits == 0 ) {
        die $self->error_message("Failed to find matching objects for '".$self->uuid."' on CG Hub!");
    }

    return 1;
}

1;

