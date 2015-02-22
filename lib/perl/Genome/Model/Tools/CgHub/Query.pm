package Genome::Model::Tools::CgHub::Query;

use strict;
use warnings;

use Genome;

use File::Temp;

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
            doc => 'Save standard output/error to this file. If given the xml file, this file will only contain the status of the query.',
        },
    },
    doc => 'Query CG Hub for metadata.',
};

sub _build_command {
    my $self = shift;

    my @cmd_parts = ( "cgquery", "'analysis_id=".$self->uuid."'" );
    if ( $self->xml_file ) {
        push @cmd_parts, "-o", $self->xml_file;
    }

    return join(' ', @cmd_parts);
}

sub _verify_success {
    my $self = shift;

    my $output_file = $self->output_file;
    if ( not -s $output_file ) {
        $self->error_message('CG Query failed to generate output file!');
        return;
    }

    my $output_fh = Genome::Sys->open_file_for_reading($output_file);
    my $number_of_hits_regexp = qr/Matching Objects\s+:\s+(\d+)/;
    my $number_of_hits_line = List::Util::first { m/$number_of_hits_regexp/ } $output_fh->getlines;
    if ( not $number_of_hits_line ) {
        $self->error_message('Failed to find matching objects line in query output!');
        return;
    }

    my ($number_of_hits) = $number_of_hits_line =~ m/$number_of_hits_regexp/;
    if ( not $number_of_hits or $number_of_hits == 0 ) {
        $self->error_message('No matching objects found for uuid! '.$self->uuid);
        return;
    }

    return 1;
}

1;

