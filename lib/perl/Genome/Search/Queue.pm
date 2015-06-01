package Genome::Search::Queue;

use Carp;
use Genome;

use strict;
use warnings;

use Params::Validate qw(:types);

class Genome::Search::Queue {
    table_name => 'web.search_index_queue',
    id_by => [
        id => { is => 'Text', len => 32 },
    ],
    has => [
        subject_class => {
            is => 'Text',
            len => 255,
            doc => 'Class of the subject to be indexed by search.',
        },
        subject_id => {
            is => 'Text',
            len => 256,
            doc => 'ID of the subject to be indexed by search.',
        },
        timestamp => {
            is => 'DateTime',
            len => 11,
            default_value => UR::Context->now,
            doc => 'Timestamp of first request. Automatically added if not provided.',
        },
        priority => {
            is => 'Integer',
            len => 1,
            default_value => 1,
            is_optional => 1,
            doc => 'Priority describes the order to process. (0 = high, 1 = normal, 2-9 = low)',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    id_generator => '-uuid',
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
        Carp::croak "subject_class ($subject_class) must be indexable in order to add to IndexQueue";
    }

    unless ($bx->specifies_value_for('timestamp')) {
        $bx = $bx->add_filter('timestamp' => UR::Context->now);
    }

    my $index_queue = $class->SUPER::create($bx);

    return $index_queue;
}

sub default_priority {
    my $class = shift;
    my $meta = $class->__meta__;
    my $property = $meta->property('priority');
    return $property->{default_value};
}

sub create_dedup_iterator {
    return $_[0]->create_iterator(
        -group_by => [qw(subject_class subject_id)],
    );
}

sub dedup_set {
    my ($class, $set) = Params::Validate::validate_pos(@_,
        { type => SCALAR },
        { isa => 'Genome::Search::Queue::Set',
          callbacks => {
              'specifies subject_id' => sub { shift->rule->value_for('subject_id') },
              'specifies subject_class' => sub { shift->rule->value_for('subject_class') },
          },
        },
    );

    my $m_iter = $set->member_iterator;
    my $q = $m_iter->next;
    my $duplicate_count = 0;
    while (my $q = $m_iter->next) {
        $duplicate_count++;
        $q->delete;
    }

    return $duplicate_count;
}

sub dedup {
    my $class = shift;

    my $max = 500;

    my $lw = UR::Context->object_cache_size_lowwater();
    my $lw_guard = Scope::Guard->new( sub { UR::Context->object_cache_size_lowwater($lw) } );
    UR::Context->object_cache_size_lowwater($UR::Context::all_objects_cache_size);

    my $hw = UR::Context->object_cache_size_highwater();
    my $hw_guard = Scope::Guard->new( sub { UR::Context->object_cache_size_highwater($hw) } );
    UR::Context->object_cache_size_highwater(UR::Context->object_cache_size_lowwater + 20 * $max);

    my $delete_count = 0;
    my $commit_and_prune = sub {
        Genome::Logger->infof("Committing the removal of %d duplicates...\n", $delete_count);
        $delete_count = 0;
        UR::Context->commit();
        UR::Context->prune_object_cache();
    };
    my $iter = $class->create_dedup_iterator();
    while (my $s = $iter->next) {
        if ($s->count > 1) {
            my $duplicate_count = $class->dedup_set($s);
            $delete_count += $duplicate_count;
            Genome::Logger->infof("Deleted %d duplicates.\n", $duplicate_count);

            if ($delete_count > $max) {
                $commit_and_prune->();
            }
        }
    }

    $commit_and_prune->();

    return 1;
}

1;
