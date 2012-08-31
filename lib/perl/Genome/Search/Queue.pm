package Genome::Search::Queue;

use Carp;
use Genome;

class Genome::Search::Queue {
    id_generator => '-uuid',
    id_by => [
        id => {
            is => 'Text',
        },
    ],
    has => [
        subject_class => {
            is => 'Text',
            doc => 'Class of the subject to be indexed by search.',
        },
        subject_id => {
            is => 'Text',
            doc => 'ID of the subject to be indexed by search.',
        },
        timestamp => {
            is => 'Time',
            doc => 'Timestamp of first request. Automatically added if not provided.',
        },
        priority => {
            is => 'Number',
            doc => 'Priority describes the order to process. (0 = high, 1 = normal, 2-9 = low)',
            default_value => 1,
        },
    ],
    data_source => 'Genome::DataSource::GMSchema',
    table_name => 'SEARCH_INDEX_QUEUE',
};

sub __display_name__ {
    my $self = shift;
    return $self->id . ' (' . $self->subject_class . ' ' . $self->subject_id . ')';
}

sub queue_iterator {
    my $class = shift;
    return $class->create_iterator(-order_by => ['priority', 'timestamp']);
}

sub create {
    my $class = shift;

    my $bx = $class->define_boolexpr(@_);

    my $subject_class = $bx->value_for('subject_class');
    unless ($subject_class) {
        Carp::croak "subject_class not specified, cannot check if it is indexable";
    }
    unless (Genome::Search->is_indexable($subject_class)) {
        Carp::croak "subject_class must be indexable in order to add to IndexQueue";
    }

    unless ($bx->specifies_value_for('timestamp')) {
        $bx = $bx->add_filter('timestamp' => UR::Context->now);
    }

    $index_queue = $class->SUPER::create($bx);

    return $index_queue;
}

sub default_priority {
    my $class = shift;
    my $meta = $class->__meta__;
    my $property = $meta->property('priority');
    return $property->{default_value};
}

1;
