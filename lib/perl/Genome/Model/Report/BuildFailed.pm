package Genome::Model::Report::BuildFailed;

use strict;
use warnings;

use Genome;
use Workflow;

use Data::Dumper 'Dumper';
require Text::Wrap;

class Genome::Model::Report::BuildFailed {
    is => 'Genome::Model::Report::BuildEventBase',
    has => [
    errors => {
        is => 'Array',
        doc => 'The errors from the build.',
    },
    ],
};

sub create {
    my ($class, %params) = @_;

    my $errors = delete $params{errors};
    my $self = $class->SUPER::create(%params)
        or return;

    unless ( $errors and @$errors ) {
        $self->error_message('No errors given to generate build fail report');
        $self->delete;
        return;
    }

    $self->errors($errors);

    return $self;
}

sub _add_to_report_xml {
    my $self = shift;

    $self->_add_build_event_dataset
        or return; # bad
    
    #< Errors >#
    my @rows;
    my @headers = (qw/
        build_event_id stage_event_id stage step_event_id step error error_wrapped error_log
        /);
    for my $error ( @{$self->errors} ) {
        push @rows, [ map { $error->$_ } @headers ];
    }

    $self->_add_dataset(
        name => 'errors',
        row_name => 'error',
        headers => [ map { s#\_#\-#g; $_ } @headers ],
        rows => \@rows,
    ) or return; # bad

    return 1;
}

1;

#$HeadURL$
#$Id$
