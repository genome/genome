package Genome::Sys::Command::Search::Index;

use Genome;

class Genome::Sys::Command::Search::Index {
    is => ['Genome::Role::Logger', 'Command'],
    has => [
        action => {
            is => 'Text',
            default => 'add',
            valid_values => ['add', 'delete'],
        },
        subject_text => {
            is => 'Text',
            shell_args_position => 1,
        },
        confirm => {
            is => 'Boolean',
            default => 1,
        },
        max_changes_per_commit => {
            is => 'Number',
            default => 50,
        },
        loop_sleep => {
            is => 'Number',
            default => 10,
        },
        tie_stderr => {
            is => 'Boolean',
            default => 1,
            doc => '(warning) globally tie STDERR to this logger',
        },
        commited_changes => {
            is => 'Number',
            default => 0,
        },
    ],
};

sub execute {
    my $self = shift;

    # Manually init logger so it ties STDERR.
    $self->log_dispatch();

    if ($self->subject_text ne 'list') {
        my $confirmed = $self->prompt_for_confirmation() if $self->confirm;
        if ($self->confirm && !$confirmed) {
            $self->info('Aborting.');
            return;
        }
    }

    if ($self->subject_text eq 'all') {
        $self->index_all;
    }
    elsif ($self->subject_text eq 'queued') {
        $self->index_queued;
    }
    elsif ($self->subject_text eq 'daemon') {
        $self->daemon;
    }
    elsif ($self->subject_text eq 'list') {
        $self->list;
    }
    else {
        die "Not able to modify specific items at this time";
    }

    return 1;
}

sub prompt_for_confirmation {
    my $self = shift;

    my $solr_server = $ENV{GENOME_SYS_SERVICES_SOLR};
    print "Are you sure you want to rebuild the index for the search server at $solr_server? ";
    my $response = <STDIN>;
    chomp $response;
    $response = lc($response);

    return ($response =~ /^(y|yes)$/);
}

sub index_all {
    my $self = shift;

    my $action = $self->action;

    my @classes_to_index = $self->indexable_classes;
    for my $class (@classes_to_index) {
        $self->info("Scanning $class...");
        my @subjects = $class->get();
        for my $subject (@subjects) {
            my $subject_class = $subject->class;
            my $subject_id = $subject->id;
            $self->modify_index($action, $subject_class, $subject_id);
        }
    }

    return 1;
}

my $signaled_to_quit;
sub daemon {
    my $self = shift;

    local $SIG{INT} = sub { print STDERR "\nDaemon will exit as soon as possible.\n"; $signaled_to_quit = 1 };
    local $SIG{TERM} = sub { print STDERR "\nDaemon will exit as soon as possible.\n"; $signaled_to_quit = 1 };

    while (!$signaled_to_quit) {
        my $pid = fork();
        if (not defined $pid) {
            die "Failed to fork";
            exit if $signaled_to_quit;
        }
        elsif ($pid == 0) {
            # CHILD
            local $SIG{INT} = sub { $signaled_to_quit = 1 };
            local $SIG{TERM} = sub { $signaled_to_quit = 1 };

            if ($signaled_to_quit) {
                $self->inf("CHILD($$): signaled to quit");
                exit;
            }

            my $max = $self->max_changes_per_commit;
            $self->info("CHILD($$): Processing index queue max=$max)");
            eval { $self->index_queued(max_changes_count => $max); };
            if ($@) { $self->info("CHILD($$): ahh shit: $@"); }

            if ($signaled_to_quit) {
                $self->info("CHILD($$): signaled to quit");
                exit;
            }

            $self->info("CHILD($$): Commiting...");
            UR::Context->commit;

            exit;
        }
        else {
            # PARENT
            waitpid($pid, 0);

            last if $signaled_to_quit;
            my $loop_sleep = $self->loop_sleep;
            $self->info("PARENT($$): Sleeping for $loop_sleep seconds...");
            sleep $loop_sleep;
        }
    }

    $self->info("Exiting...");
    return 1;
}

sub list {
    my $self = shift;

    my $index_queue_iterator = Genome::Search::Queue->queue_iterator();

    print join("\t", 'PRIORITY', 'TIMESTAMP', 'SUBJECT_CLASS', 'SUBJECT_ID') . "\n";
    print join("\t", '--------', '---------', '-------------', '----------') . "\n";
    while (my $index_queue_item = $index_queue_iterator->next) {
        print join("\t",
            $index_queue_item->priority,
            $index_queue_item->timestamp,
            $index_queue_item->subject_class,
            $index_queue_item->subject_id,
        ) . "\n";
    }

    return 1;
}

sub index_queued {
    my $self = shift;
    my %params = @_;

    my $max_changes_count = delete $params{max_changes_count};

    # TODO Should optimize this by grouping by subject id and class and removing all related rows
    my $index_queue_iterator = Genome::Search::Queue->queue_iterator();

    my $subject_seen = {};
    my $modified_count = 0;
    while (
        !$signaled_to_quit
        && (!defined($max_changes_count) || $modified_count < $max_changes_count)
        && (my $index_queue_item = $index_queue_iterator->next)
    ) {
        my $subject_class = $index_queue_item->subject_class;
        my $subject_id = $index_queue_item->subject_id;
        last if $signaled_to_quit;

        # if we've already seen this subject during this iterator then we do not need to re-index it
        if ($subject_seen->{$subject_class}->{$subject_id}) {
            $index_queue_item->delete();
        }
        else {
            my $action;
            if (not $subject_class->can('get')) {
                $self->warning_message("Class ($subject_class) cannot 'get'. Deleting item (ID: $subject_id) from queue.");
                $action = 'delete';
            } else {
                $action = ($subject_class->get($subject_id) ? 'add' : 'delete');
            }
            if($action eq 'add' and $subject_class->isa('UR::Object::Set')) {
                my $set = $subject_class->get($subject_id);
                unless ($set->members) {
                    $action = 'delete';
                }
            }
            last if $signaled_to_quit;
            if ($self->modify_index($action, $subject_class, $subject_id)) {
                $subject_seen->{$subject_class}->{$subject_id}++;
                $index_queue_item->delete();
                $modified_count++;
            } else {
                # Move it to the back of the line.
                $index_queue_item->timestamp(UR::Context->now);
            }
        }
    }

    return 1;
}

sub modify_index {
    my ($self, $action, $subject_class, $subject_id) = @_;

    my $display_name = "(Class: $subject_class, ID: $subject_id)";

    my ($rv, $error);
    if ($action eq 'add') {
        $rv = eval {
            my $subject = $subject_class->get($subject_id);
            return if $signaled_to_quit;
            unless ($subject) { die "Could not get object $display_name" };
            Genome::Search->add($subject);
        };
        $error = $@;
    }
    elsif ($action eq 'delete') {
        $rv = eval { Genome::Search->delete_by_class_and_id($subject_class, $subject_id) };
        $error = $@;
    }

    if ($error) {
        $self->error($error);
    }

    if ($rv) {
        my $display_action = ($action eq 'add' ? 'Added' : 'Deleted');
        $self->info("$display_action $display_name");
    }
    else {
        my $display_action = ($action eq 'add' ? 'Failed to add' : 'Failed to delete');
        $self->info("$display_action $display_name\n$@");
    }

    return $rv;
}

sub indexable_classes {
    my $self = shift;

    my @searchable_classes = Genome::Search->searchable_classes();

    my @classes_to_add;
    for my $class (@searchable_classes) {
        eval "use $class";
        my $use_errors = $@;
        if ($use_errors) {
            $self->debug("Class ($class) in searchable_classes is not usable ($use_errors).");
            next;
        }

        my $class_is_indexable = Genome::Search->is_indexable($class);
        if (!$class_is_indexable) {
            $self->debug("Class ($class) in searchable_classes is not indexable.");
            next;
        }

        push @classes_to_add, $class;
    }

    return @classes_to_add;
}

sub delete_view_and_aspects {
    my ($self, $view) = @_;
    for my $aspect ($view->aspects) {
        if (my $delegate_view = $aspect->delegate_view) {
            $self->delete_view_and_aspects($delegate_view);
        }
        $aspect->delete;
    }
    $view->delete;
}

1;
