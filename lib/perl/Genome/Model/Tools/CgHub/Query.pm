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
            doc => 'Save standard output/error to this file. If given the xml file, this fie will only contain the status of the query.',
        },
    },
    has_transient => {
        _was_output_file_given => { is => 'Boolean', },
    },
    doc => 'Query CG Hub for metadata.',
};

sub _build_command {
    my $self = shift;

    my @cmd_parts = ( "cgquery", "analysis_id=".$self->uuid );
    if ( $self->xml_file ) {
        push @cmd_parts, "-o", $self->xml_file;
    }

    my $output_file = $self->output_file;
    if ( not $output_file ) {
        $output_file = $self->output_file( File::Temp::tempfile(CLEANP => 1) );
        $self->_was_output_file_given(1);
    }
    push @cmd_parts, '>', $output_file;

    return join(' ', @cmd_parts);
}

sub _verify_success {
    my $self = shift;

    my $output_file = $self->output_file;
    if ( not -s $output_file ) {
        die $self->error_message('CG Query failed to generate output file!');
    }

    my $print_output;
    if ( $self->_was_output_file_given ) {
        my $fh = Genome::Sys->open_file_for_reading($output_file);
        $print_output = sub{ $fh->print($_); };
    }
    else {
        $print_output = sub{ print $_; };
    }

    my $output_fh = Genome::Sys->open_file_for_reading($output_file);
    my $number_of_hits = 0;
    while ( my $line = $output_fh->getline ) {
        $print_output->($line);
        if ( $line =~ /Matching Objects\s+:\s+(\d+)/ ) {
            $number_of_hits = $1;
        }
    }
    $output_fh->close;

    if ( $number_of_hits == 0 ) {
        die $self->error_message("Failed to find matching objects for '".$self->uuid."' on CG Hub!");
    }

    return 1;
}

1;

