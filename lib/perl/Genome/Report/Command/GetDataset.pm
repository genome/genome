package Genome::Report::Command::GetDataset;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use IO::Handle;

class Genome::Report::Command::GetDataset {
    is => 'Genome::Report::Command',
    has => [ 
        dataset_name => {
            is => 'Text',
            doc => 'Name of the dataset to get.',
        },
    ],
    has_optional => [
        output_file => {
            is => 'Text',
            doc => 'File name to save output.  If none given, prints to STDOUT',
        },
        output_type => {
            is => 'Text', 
            default_value => default_output_type(),
            doc => 'Get the dataset as this format: '.join(', ', output_types()).'.  Deafult is '.default_output_type(),
        },
        _output => {
            is_input => 0,
            doc => 'Output object',
        },

    ],
};

#< Output Types >#
sub output_types {
    return (qw/ csv xml /);
}

sub default_output_type {
    return (output_types)[0];
}

#< Helps >#
sub help_brief {
    'Get a dataset from a report, then save or print';
}

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}

#< Command >#
sub execute {
    my $self = shift;

    # report
    my $report = $self->report;
    unless ( $report ) {
        $self->error_message(
            sprintf('Report name (%s) not found for build (%s)', $self->report_name, $self->build_id)
        );
        return;
    }

    # dataset
    unless ( $self->dataset_name ) {
        $self->error_message("No dataset name given to create");
        return;
    }
    my $dataset = $report->get_dataset( $self->dataset_name );
    unless ( $dataset ) {
        $self->error_message( sprintf('Dataset (%s) not found in report (%s)', $self->dataset_name, $report->name) );
        return;
    }

    # output handler
    if ( $self->output_file ) { 
        my $fh = Genome::Sys->open_file_for_writing($self->output_file);
        unless ( $fh ) {
            return;
        }
        $fh->autoflush;
        $self->_output($fh);
    }
    else { # STDOUT
        my $handle = IO::Handle->new();
        $handle->fdopen(fileno(STDOUT), "w");
        $handle->autoflush;
        $self->_output($handle);
    }

    # generate and print
    my $method = '_'.$self->output_type;
    return $self->_output->print( $self->$method($dataset) );
}

sub _xml {
    return $_[1]->to_xml_string."\n";
}

sub _csv {
    return $_[1]->to_separated_value_string(separator => ',');
}

1;

#$HeadURL$
#$Id$
