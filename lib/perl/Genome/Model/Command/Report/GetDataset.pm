package Genome::Model::Command::Report::GetDataset;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';
use IO::Handle;

class Genome::Model::Command::Report::GetDataset {
    is => 'Genome::Model::Command::Report',
    has => [ 
    dataset_name => {
        is => 'Text',
        doc => 'Name of the dataset to get.',
    },
    _output => {
        is_input => 0,
        doc => 'Output object',
    },
    ],
    has_optional => [
    output_file => {
        is => 'Text',
        doc => 'File name to save output.  If none given, prints to STDOUT',
    },
    output_type => {
        is => 'Text', 
        default_value => __PACKAGE__->_default_output_type,
        doc => 'Get the dataset as this format: '.join(', ', __PACKAGE__->_output_types).'.  Deafult is '.__PACKAGE__->_default_output_type,
    },
    ],
};

#< Output Types >#
sub _output_types {
    return (qw/ csv xml /);
}

sub _default_output_type {
    return (_output_types)[0];
}

#< Helps >#
sub help_brief {
    'Get a dataset from a report';
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
sub create {
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    unless ( $self->report ) {
        $self->error_message(
            sprintf('Report name (%s) not found for build (%s)', $self->report_name, $self->build_id)
        );
        $self->delete;
        return;
    }

    unless ( $self->dataset_name ) {
        $self->error_message("No dataset name given to create");
        $self->delete;
        return;
    }

    if ( $self->output_file ) { 
        my $fh = Genome::Sys->open_file_for_writing($self->output_file);
        unless ( $fh ) {
            $self->delete;
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

    return $self;
}

sub execute {
    my $self = shift;
    
    my $report = $self->report
        or return;

    my ($dataset) = $report->get_dataset_nodes_for_name( $self->dataset_name )
       or return; # TODO handle mutliple datasets with the same name

    my $method = '_'.$self->output_type;
    return $self->_output->print( $self->$method($dataset) );
}

sub _xml {
    return $_[1]->toString."\n";
}

sub _csv {
    my ($self, $dataset) = @_;

    my @rows = grep { $_->nodeType == 1 } $dataset->findnodes('*');
    # headers
    my $csv = join(',', map { $_->nodeName } grep { $_->nodeType == 1 } $rows[0]->findnodes('*'))."\n";
    # data
    for my $row ( @rows ) {
        $csv .= join(',', map { $_->to_literal } grep { $_->nodeType == 1 } $row->getChildnodes)."\n";
    }

    return $csv;
}

1;

#$HeadURL$
#$Id$
