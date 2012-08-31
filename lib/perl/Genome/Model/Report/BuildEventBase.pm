package Genome::Model::Report::BuildEventBase;

use strict;
use warnings;

use Genome;

use Data::Dumper 'Dumper';

class Genome::Model::Report::BuildEventBase {
    is => 'Genome::Model::Report',
    is_abstract => 1,
    has => [
    name => {
        calculate => q| 
        return 'Build '.ucfirst($self->event_name),
        |,
        is_constant => 1,
    },
    description => {
        calculate => q| 
        return sprintf(
            'Build %s for Model (%s)',
            ucfirst($self->event_name),
            $self->model_name,
        );
        |,
        is_constant => 1,
    },
    event_name => {
        calculate => q| 
        my $class = $self->class;
        $class =~ s#Genome::Model::Report::Build##;
        return $class;
        |,
        is_constant => 1,
    },
    ],
};

sub _add_build_event_dataset {
    my $self = shift;

    my $build_event = $self->build->build_event;
    $self->_add_dataset(
        name => 'build-events',
        row_name => 'build-event',
        headers => [qw/ id status date-scheduled date-completed /],
        rows => [ [ map { $build_event->$_ } (qw/ id event_status date_scheduled date_completed /) ] ],
    ) or return;

    return 1;
}

1;

#$HeadURL$
#$Id$
