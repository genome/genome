package Genome::Report::FromSeparatedValueFile;
#:adukes check

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Report::FromSeparatedValueFile {
    is => 'Genome::Report::Generator',
    has => [
    name => {
        is => 'Text',
        doc => 'Name to give the report.  Will usually have a default/calculated value',
    },
    description => {
        is => 'Text',
        doc => 'Description to give the report.  Will usually have a default/calculated value',
    },
    svr => {
        is => 'Genome::Utility::IO::SeparatedValueReader',
        doc => 'Separated value reader that file to import',
    },
    ],
};

#< Create >#
sub create { 
    my $class = shift;

    my $self = $class->SUPER::create(@_)
        or return;

    unless ( $self->name ) {
        $self->error_message("Name is required to create");
        $self->delete;
        return;
    }

    unless ( $self->description ) {
        $self->error_message("Description is required to create");
        $self->delete;
        return;
    }
    
    unless ( $self->svr ) {
        $self->error_message('Separated value reader (svr) is required to create');
        $self->delete;
        return;
    }

    return $self;
}

#< Generate >#
sub _add_to_report_xml {
    my $self = shift;

    my $svr = $self->svr;
    my @rows;
    while ( my $ref = $svr->next ) {
        push @rows, [ map { $ref->{$_} } @{$svr->headers} ];
    }

    unless ( @rows ) {
        $self->error_message("No data found in separated value reader: ".$self->orignal_input);
        return;
    }

    $self->_add_dataset(
        name => 'dataset',
        headers => $svr->headers,
        row_name => 'row',
        rows => \@rows,
    );
    
    return 1;
}

1;

#$HeadURL$
#$Id$
