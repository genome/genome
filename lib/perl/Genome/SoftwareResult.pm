package Genome::SoftwareResult;

use strict;
use warnings;

use Genome;
use Genome::Sys::LockProxy qw();
use Digest::MD5 qw(md5_hex);
use File::Basename qw(fileparse);
use Data::Dumper;
use Date::Manip;
use List::MoreUtils qw(uniq);
use File::Grep qw(fgrep);

use Carp;

use JSON;

use Genome::Utility::Instrumentation;

require Genome::Model::Build;

class Genome::SoftwareResult {
    is_abstract => 1,
    table_name => 'result.software_result',
    subclass_description_preprocessor => 'Genome::SoftwareResult::_expand_param_and_input_properties',
    subclassify_by => 'subclass_name',
    id_generator => '-uuid',
    id_by => [
        id => { is => 'Text', len => 32 },
    ],
    attributes_have => [
        is_param => { is => 'Boolean', is_optional=>'1' },
        is_input => { is => 'Boolean', is_optional=>'1' },
        is_metric => { is => 'Boolean', is_optional=>'1' }
    ],
    has => [
        module_version      => { is => 'Text', len => 64, column_name => 'VERSION', is_optional => 1 },
        subclass_name       => { is => 'Text', len => 255, column_name => 'CLASS_NAME' },
        lookup_hash         => { is => 'Text', len => 32, column_name => 'LOOKUP_HASH', is_optional => 1 },
        output_dir          => { is => 'Text', len => 1000, column_name => 'OUTPUTS_PATH', is_optional => 1 },
        test_name           => { is_param => 1, is_delegated => 1, is_mutable => 1, via => 'params', to => 'value_id', where => ['name' => 'test_name'], is => 'Text', doc => 'Assigns a testing tag to the result.  These will not be used in default processing', is_optional => 1 },
        _lock_proxy          => { is_optional => 1, is_transient => 1 },
    ],
    has_many_optional => [
        params              => { is => 'Genome::SoftwareResult::Param', reverse_as => 'software_result'},
        inputs              => { is => 'Genome::SoftwareResult::Input', reverse_as => 'software_result'},
        metrics             => { is => 'Genome::SoftwareResult::Metric', reverse_as => 'software_result'},
        users               => { is => 'Genome::SoftwareResult::User', reverse_as => 'software_result'},
        disk_allocations    => { is => 'Genome::Disk::Allocation', reverse_as => 'owner'},
        builds              => { is => 'Genome::Model::Build', via => 'users', to => 'user', where => ["user_class_name like" => "Genome::Model::Build%"] },
        build_ids           => { via => 'builds', to => 'id', is_deprecated => 1 },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
    doc => 'base class for managed data sets, with database tracking for params, inputs, metrics, and disk',
};

UR::Observer->register_callback(
    subject_class_name => 'Genome::SoftwareResult',
    aspect => 'subclass_loaded',
    callback => sub {
        my ($classname, $aspect, $subclassname) = @_;

        my $subclass = $subclassname->__meta__;
        unless ($subclass) {
            die "Failed to get object Type for $subclassname";
        }

        # classes that have table_name will complain if no column exists in DB for a property
        unless ($subclass->table_name) {
            my @properties = grep { $_->class_name->isa(__PACKAGE__) }
                            $subclass->properties();
            # try to define what "ambiguous" means...
            my @ambiguous_properties = grep { !(
                    $_->property_name eq 'subclass_name'
                    || $_->is_param
                    || $_->is_input
                    || $_->is_calculated
                    || $_->is_delegated
                    || $_->class_name ne $subclassname
                    || $_->is_transient
                )} @properties;

            if (@ambiguous_properties) {
                die sprintf("ambiguous properties: %s\n",
                        join(", ", map { $_->property_name } @ambiguous_properties));
            }
        }
    },
);

our %LOCKS;

sub __display_name__ {
    my $self = shift;
    return $self->class . ' (' . $self->id . ')';
}

sub _faster_get {
    my $class = shift;
    my @args = @_;

    my $statsd_prefix = "software_result_get.full_time.";
    my $statsd_class_suffix = "$class";
    $statsd_class_suffix =~ s/::/_/g;
    my @objects;

    Genome::Utility::Instrumentation::timer(
        $statsd_prefix . 'total',
        $statsd_prefix . $statsd_class_suffix,
        sub {
            my $lookup_hash = $class->calculate_lookup_hash_from_arguments(@args);
            # NOTE we do this so that get can be noisy when called directly.
            @objects = $class->SUPER::get(lookup_hash => $lookup_hash);
            if (@objects > 1) {
                die $class->error_message(
                    'Multiple results found for lookup hash %s!',
                    $lookup_hash
                );
            }

            if (@objects) {
                my $calculated_lookup_hash = $objects[0]->calculate_lookup_hash();
                my $lookup_hash_property = $objects[0]->lookup_hash;
                if ($calculated_lookup_hash ne $lookup_hash_property) {
                    die $class->error_message(
                        "SoftwareResult lookup_hash (%s) does not match it's calculated lookup_hash (%s).  Cannot trust that the correct SoftwareResult was retrieved.",
                        $lookup_hash_property,
                        $calculated_lookup_hash
                    );
                }
            }
        }
    );

    return $objects[0];
}

# Override get to be noisy.  This should generally not be called.
sub get {
    my $class = shift;

    # We only want to apply this restriction if we're asking for a scalar.
    if (defined wantarray && not wantarray) {
        # UR::Context::query has special logic for when it is called on an object
        unless (ref $class) {
            if (not $class->_is_id_only_query(@_)) {
                # XXX This should be a warning, but warnings are not yet sent to logstash
                $class->error_message(
                    "Calling get on SoftwareResult (unless getting by id) is slow and possibly incorrect.\n"
                    . Carp::longmess);
            }
        }
    }

    return $class->SUPER::get(@_);
}

sub _is_id_only_query {
    my $class = shift;

    my ($bx, @extra) = $class->define_boolexpr(@_);
    return if @extra;
    return $bx->is_id_only;
}

# You must specify enough parameters to uniquely identify an object to get a result.
# If two users specify different sets of parameters that uniquely identify the result,
# they will create different locks.
sub get_with_lock {
    my $class = shift;

    return Genome::SoftwareResult::User->with_registered_users(
        $class->_preprocess_params_for_callback(@_),
        callback => sub {
            my $params_processed = $class->_gather_params_for_get_or_create($class->_preprocess_params_for_get_or_create(@_));

            my %is_input = %{$params_processed->{inputs}};
            my %is_param = %{$params_processed->{params}};

            # Only try with lock if object does not exist since locking causes
            # a performance hit. It is assumed that if an object is found it is
            # complete. If this is a bad assumption then we need to add a
            # status to SoftwareResults.
            my $lock;
            my $result = $class->_faster_get(@_);
            unless ($result) {
                my $subclass = $params_processed->{subclass};
                my $lookup_hash = $subclass->calculate_lookup_hash_from_arguments(@_);

                unless ($lock = $subclass->_lock($lookup_hash)) {
                    die "Failed to get a lock for " . Dumper(@_);
                }

                UR::Context->current->reload($subclass, lookup_hash => $lookup_hash);

                eval {
                    $result = $subclass->_faster_get(@_);
                };
                my $error = $@;

                if ($error) {
                    $subclass->_release_lock_or_die($lock, "Failed to unlock during get_with_lock.");
                    die $subclass->error_message('Failed in get! %s', $error);
                }
            }


            if ($result && $lock) {
                $result->_lock_proxy($lock);

                $result->debug_message("Cleaning up lock %s...", $lock->resource);
                unless ($result->_unlock) {
                    $result->error_message("Failed to unlock after getting software result");
                    die "Failed to unlock after getting software result";
                }
                $result->debug_message("Cleanup completed for lock %s.", $lock->resource);
            } elsif ($lock) {
                $class->_release_lock_or_die($lock, "Failed to unlock after not finding software result.");
            }

            $result->_auto_unarchive if $result;
            return ($result, 0);
        }
    );

}

sub get_or_create {
    my $class = shift;

    return Genome::SoftwareResult::User->with_registered_users(
        $class->_preprocess_params_for_callback(@_),
        callback => sub {
            my $params_processed = $class->_gather_params_for_get_or_create($class->_preprocess_params_for_get_or_create(@_));

            my %is_input = %{$params_processed->{inputs}};
            my %is_param = %{$params_processed->{params}};

            my $result = $class->_faster_get(@_);
            my $newly_created = 0;

            unless ($result) {
                $result = $class->create(@_);
                if ($result) {
                    $newly_created = 1;
                } else {
                    # see if the reason we failed was b/c the objects were created while we were locking...
                    $result = $class->_faster_get(@_);
                }
            }

            unless ($result) {
                die $class->error_message(
                    'Could not create a %s for params %s even after trying!',
                    $class,
                    Data::Dumper::Dumper(\@_),
                );
            }

            $result->_auto_unarchive();
            return ($result, $newly_created);
        });
}

#hook for subclasses
sub _preprocess_params_for_callback {
    my $class = shift;

    return @_;
}

sub create {
    my $class = shift;

    if ($class eq __PACKAGE__ || $class->__meta__->is_abstract) {
        # this class is abstract, and the super-class re-calls the constructor from the correct subclass
        return $class->SUPER::create(@_);
    }

    my $params_processed = $class->_gather_params_for_get_or_create($class->_preprocess_params_for_get_or_create(@_));

    my %is_input = %{$params_processed->{inputs}};
    my %is_param = %{$params_processed->{params}};

    my $previously_existing = $class->_faster_get(@_);

    if (defined($previously_existing)) {
        $class->error_message("Attempt to create an $class but it looks like we already have one with those params " . Dumper(\@_));
        return;
    }

    my $lock;
    my $lookup_hash = $class->calculate_lookup_hash_from_arguments(@_);
    unless ($lock = $class->_lock($lookup_hash)) {
        die "Failed to get a lock for " . Dumper(@_);
    }

    # we might have had to wait on the lock, in which case someone else was probably creating that entity
    # do a "reload" here to force another trip back to the database to see if a software result was created
    # while we were waiting on the lock.
    ($previously_existing) = UR::Context->current->reload($class,
        lookup_hash => $class->calculate_lookup_hash_from_arguments(@_));

    if (defined($previously_existing)) {
        $class->error_message("Attempt to create an $class but it looks like we already have one with those params " . Dumper(\@_));
        $class->_release_lock_or_die($lock, "Failed to release lock in create before committing SoftwareResult.");
        return;
    }

    # We need to update the indirect mutable accessor logic for non-nullable
    # hang-offs to delete the entry instead of setting it to null.  Otherwise
    # we get SOFTWARE_RESULT_PARAM entries with a NULL, and unsavable PARAM_VALUE.
    # also remove empty strings because that's equivalent to a NULL to the database

    # Do the same for inputs (e.g. alignment results have nullable segment values for instrument data, which are treated as inputs)
    my @param_remove = grep { not (defined $is_param{$_}) || $is_param{$_} eq "" } keys %is_param;
    my @input_remove = grep { not (defined $is_input{$_}) || $is_input{$_} eq "" } keys %is_input;
    my $bx = $class->define_boolexpr($class->_preprocess_params_for_get_or_create(@_));
    for my $i (@param_remove, @input_remove) {
        $bx = $bx->remove_filter($i);
    }

    my $self = $class->SUPER::create($bx);
    unless ($self) {
        $class->_release_lock_or_die($lock,"Failed to unlock during create after committing SoftwareResult.");
        return;
    }

    $self->_lock_proxy($lock);

    $self->_add_unlock_observers;

    if (my $output_dir = $self->output_dir) {
        if (-d $output_dir) {
            my @files = glob("$output_dir/*");
            if (@files) {
                $self->delete;
                die "Found files in output directory $output_dir!:\n\t"
                    . join("\n\t", @files);
            }
            else {
                $self->debug_message("No files in $output_dir.");
            }
        }
        else {
            $self->debug_message("Creating output directory $output_dir...");
            eval {
                Genome::Sys->create_directory($output_dir)
            };
            if (my $error = $@) {
                $self->delete;
                die $error;
            }
        }
    }

    $self->module_version($self->resolve_module_version) unless defined $self->module_version;
    $self->subclass_name($class);
    $self->lookup_hash($self->calculate_lookup_hash());
    return $self;
}

sub _add_unlock_observers {
    my $self = shift;

    my $unlock_callback = sub {
        $self->_unlock;
    };
    my @observers =
        $self->add_observer(aspect => 'commit', callback => $unlock_callback, once => 1),
        $self->add_observer(aspect => 'delete', callback => $unlock_callback);

    return @observers;
}

sub _gather_params_for_get_or_create {
    my $class = shift;

    my ($bx, @extra) = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, @_);
    die sprintf('got extra parameters: [%s]', join(',', @extra)) if @extra;

    my %params = $bx->params_list;
    my %is_input;
    my %is_param;
    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        if ($meta->{is_input} && exists $params{$key}) {
            $is_input{$key} = $params{$key};
        } elsif ($meta->{is_param} && exists $params{$key}) {
            $is_param{$key} = $params{$key};
        }

    }

    my %software_result_params = (subclass_name=>$class);

    return {
        software_result_params => \%software_result_params,
        subclass => $class,
        inputs=>\%is_input,
        params=>\%is_param,
    };
}

sub _preprocess_params_for_get_or_create {
    my $class = shift;
    if(scalar @_ eq 1) {
        return @_; #don't process a UR::BoolExpr or plain ID
    }

    my %params = @_;

    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);

        for my $t ('input', 'param') {
            if ($meta->{'is_' . $t} && $meta->is_many) {
                my $value_list = delete $params{$key};
                if((defined $value_list) && (scalar @$value_list)) {
                    my @values = sort map { Scalar::Util::blessed($_)? $_->id : $_ } @$value_list;
                    my $t_params = $params{$t . 's'} || [];
                    for my $i (0..$#values) {
                        my $value = $values[$i];
                        push @$t_params, {'name' => $key . '-' . $i, 'value_id' => $value};
                    }
                    $params{$t . 's'} = $t_params;

                    $params{$key . '_count'} = scalar @values;
                    $params{$key . '_md5'} = md5_hex( join(':', @values));
                } else {
                    $params{$key . '_count'} = 0;
                    $params{$key . '_md5'} = undef;
                }
            }
        }
    }
    return %params;
}

my $recalculate_lookup_hash_callback = sub {
    my ($object, $aspect) = @_;
    return unless $object;
    return unless $object->software_result;
    $object->software_result->recalculate_lookup_hash();
};
for my $name (qw(Param Input)) {
    my $classname = join('::', qw(Genome SoftwareResult), $name);
    for my $aspect (qw(create value_id)) {
        UR::Observer->register_callback(subject_class_name => $classname, aspect => $aspect, callback => $recalculate_lookup_hash_callback);
    }
}

sub calculate_query {
    my $self = shift;

    # Pre-fetch things, since we're going loop through them.
    $self->inputs;
    $self->params;

    my @query;

    my $class_object = $self->__meta__;
    for my $key ($self->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        next unless $meta->{is_input} or $meta->{is_param};
        next if $key =~ /(.+?)_(?:md5|count)$/ and $class_object->property_meta_for_name($1); #TODO remove these params completely!

        if($meta->is_many) {
            push @query,
                $key, [$self->$key];
        } else {
            push @query,
                $key, $self->$key;
        }
    }

    return @query;
}

sub calculate_lookup_hash_from_arguments {
    my $class = shift;

    my %processed_params = $class->_process_params_for_lookup_hash(@_);
    return $class->_generate_lookup_hash(\%processed_params);
}

sub calculate_lookup_hash {
    my $self = shift;

    my @query = $self->calculate_query;
    return $self->calculate_lookup_hash_from_arguments(@query);
}

sub recalculate_lookup_hash {
    my $self = shift;
    return $self->lookup_hash($self->calculate_lookup_hash);
}

sub _process_params_for_lookup_hash {
    my $class = shift;
    my %initial_params;

    # Handle the case of a boolean expression (used by _faster_get)
    if (1 == scalar(@_)) {
        %initial_params = $_[0]->params_list;
    } else {
        %initial_params = @_;
    }

    $class->_modify_params_for_lookup_hash(\%initial_params);

    my ($bx, @extra) = $class->define_boolexpr(%initial_params);

    die sprintf('got extra parameters: [%s]', join(',', @extra)) if @extra;

    my %params = $bx->params_list;

    my $class_object = $class->__meta__;
    for my $key ($class->property_names) {
        my $meta = $class_object->property_meta_for_name($key);
        unless ($meta->{is_input} or $meta->{is_param}) {
            delete $params{$key};
            next;
        }

        if ($meta->is_transient) {
            delete $params{$key};
            next;
        }

        if (defined($meta->default_value) and not exists $params{$key}) {
            $params{$key} = $meta->default_value;
        }

        unless ($meta->is_optional or $meta->is_many or exists $params{$key}) {
            confess('incomplete object specification: missing ' . $key);
        }

        for my $t ('input', 'param') {
            if ($meta->{'is_' . $t} && $meta->is_many) {
                my $value_list = delete $params{$key};
                if((defined $value_list) && (scalar @$value_list)) {
                    my @values = sort map { _resolve_object_id($_) } @$value_list;
                    $params{$key} = \@values;
                }
            } else {
                $params{$key} = _resolve_object_id($params{$key});
            }
        }

        if(not defined $params{$key} or $params{$key} eq '') {
            delete $params{$key};
        }

        next unless exists $params{$key};

        if(defined($meta->data_type) and
                ($meta->data_type eq 'Boolean' or
                 $meta->data_type eq 'UR::Value::Boolean') and $params{$key} eq 0){
            delete $params{$key};
        }
    }

    return %params;
}

sub _modify_params_for_lookup_hash {
    # overridden in some subclasses
}

sub _resolve_object_id {
    my $object = shift;

    return Scalar::Util::blessed($object) ? $object->id : $object;
}

sub _generate_lookup_hash {
    my $class = shift;
    my $hash_to_encode = shift;

    my $json = JSON->new();
    $json->canonical([1]);
    my $result = $json->encode($hash_to_encode);

    return Genome::Sys->md5sum_data($result);
}

sub remove_test_name {
    my $self = shift;
    my $param = $self->params(name => 'test_name');
    return $param->delete;
}

sub resolve_module_version {
    my $class = shift;
    my $revision = Genome::Sys->snapshot_revision;
    # the revision may be a series of long paths
    # truncate the revision if it is too long to store, but put a * at the end so we can tell this is the case
    my $pmeta = $class->__meta__->property("module_version");
    my $len = $pmeta->data_length - 1;
    if (length($revision) > $len) {
        $revision = substr($revision,0,$len) . '*';
    }
    return $revision;
}

sub _prepare_output_directory {
    my $self = shift;

    my ($allocation) = $self->disk_allocations;
    if ( $allocation ) {
        my $absolute_path = $allocation->absolute_path;
        $self->output_dir($absolute_path) if $self->output_dir and $self->output_dir ne $absolute_path;
        return $self->output_dir;
    }

    my $subdir = $self->resolve_allocation_subdirectory;
    unless ( $subdir ) {
        die $self->error_message("failed to resolve subdirectory for output data.  cannot proceed.");
    }

    my %allocation_create_parameters = (
        disk_group_name => $self->resolve_allocation_disk_group_name,
        allocation_path => $subdir,
        kilobytes_requested => $self->resolve_allocation_kilobytes_requested,
        owner_class_name => $self->class,
        owner_id => $self->id
    );
    $allocation = Genome::Disk::Allocation->allocate(%allocation_create_parameters);
    unless ($allocation) {
        die $self->error_message("Failed to get disk allocation with params:\n". Data::Dumper::Dumper(%allocation_create_parameters));
    }

    my $output_dir = $allocation->absolute_path;
    unless (-d $output_dir) {
        die $self->error_message("Allocation path $output_dir doesn't exist!");
    }

    $self->output_dir($output_dir);
    return $output_dir;
}

sub _auto_unarchive {
    my $self = shift;

    my (@allocations) = grep { $_->is_archived } $self->disk_allocations;
    return unless @allocations;

    $self->status_message('Result %s is archived... unarchiving now.', $self->id);

    my $user = Genome::Sys->username;
    my ($package, undef, $line) = caller;
    my $reason = sprintf('Automatically unarchiving due to request from %s line %s for software result %s run by %s', $package, $line, $self->id, $user);

    my $build = Genome::Model::Build::get_build_id();
    $reason .= sprintf(' for build %s', $build) if $build;

    for my $allocation (@allocations) {
        $allocation->unarchive(reason => $reason);
    }

    UR::Context->reload($self);

    return 1;
}

sub _expand_param_and_input_properties {
    my ($class, $desc) = @_;
    for my $t ('input','param','metric') {
        for my $prop_name (sort keys %{ $desc->{has} }) {
            my $prop_desc = $desc->{has}{$prop_name};

            if (exists $prop_desc->{'is_'.$t} and $prop_desc->{'is_'.$t}) {
                my $is_many = ($t ne 'metric' and exists $prop_desc->{'is_many'} and $prop_desc->{'is_many'});

                my $name_name;
                if ($t eq 'metric') {
                    $prop_desc->{'to'} = 'metric_value';
                    $name_name = 'metric_name';
                }
                else {
                    # TODO This logic was borrowed in the Model.pm's _resolve_to_for_prop_desc so
                    # when this is refactored, that should also be updated.
                    if (exists $prop_desc->{'data_type'} and $prop_desc->{'data_type'}) {
                        my $prop_class = UR::Object::Property->_convert_data_type_for_source_class_to_final_class(
                            $prop_desc->{'data_type'},
                            $class
                        );
                        if ($prop_class->isa("UR::Value")) {
                            $prop_desc->{'to'} = 'value_id';
                        } else {
                            $prop_desc->{'to'} = 'value_obj';
                        }
                    }
                    else {
                        $prop_desc->{'to'} = 'value_id';
                    }
                    $name_name = 'name';
                }

                $prop_desc->{'is_delegated'} = 1;

                if($is_many) {
                    $prop_desc->{'where'} = [
                        $name_name . ' like' => $prop_name . '-%',
                    ];
                }
                else {
                    $prop_desc->{'where'} = [
                        $name_name => $prop_name
                    ];
                }


                $prop_desc->{'is_mutable'} = 1;
                $prop_desc->{'via'} = $t.'s';

                if($is_many) {
                    my $md5_name = $prop_name . '_md5';
                    unless(exists $desc->{has}{$md5_name}) {
                        my $md5_prop = {};
                        $md5_prop->{'is'} = 'Text';
                        $md5_prop->{'is_param'} = 1;
                        $md5_prop->{'is_delegated'} = 1;
                        $md5_prop->{'via'} = 'params';
                        $md5_prop->{'to'} = 'value_id';
                        $md5_prop->{'where'} = [ 'name' => $md5_name ];
                        $md5_prop->{'doc'} = 'MD5 sum of the sorted list of values for ' . $prop_name;
                        $md5_prop->{'is_mutable'} = 1;
                        $md5_prop->{'is_optional'} = 1;

                        $md5_prop->{'property_name'} = $md5_name;
                        $md5_prop->{'class_name'} = $desc->{class_name};
                        $desc->{has}{$md5_name} = $md5_prop;
                    }

                    my $count_name = $prop_name . '_count';
                    unless(exists $desc->{has}{$count_name}) {
                        my $count_prop = {};
                        $count_prop->{'is'} = 'Number';
                        $count_prop->{'is_param'} = 1;
                        $count_prop->{'is_delegated'} = 1;
                        $count_prop->{'via'} = 'params';
                        $count_prop->{'to'} = 'value_id';
                        $count_prop->{'where'} = [ 'name' => $count_name ];
                        $count_prop->{'doc'} = 'number of values for ' . $prop_name;
                        $count_prop->{'is_mutable'} = 1;
                        $count_prop->{'is_optional'} = 1;

                        $count_prop->{'property_name'} = $count_name;
                        $count_prop->{'class_name'} = $desc->{class_name};
                        $desc->{has}{$count_name} = $count_prop;
                    }
                }
            }
        }
    }
    return $desc;
}

sub delete {
    my $self = shift;

    my $class_name = $self->class;
    my @users = $self->users;
    my @active_users = grep{$_->active} @users;
    if (@active_users) {
        my $name = $self->__display_name__;
        die "Refusing to delete $class_name $name as it still has active users:\n\t"
            .join("\n\t", map { $_->user_class_name . "\t" . $_->user_id } @active_users);
    }

    my @to_nuke = ($self->params, $self->inputs, $self->metrics, @users);

    #If we use any other results, unregister ourselves as users
    push @to_nuke, Genome::SoftwareResult::User->get(user_class_name => $class_name, user_id => $self->id);

    for (@to_nuke) {
        unless($_->delete) {
            die "Failed to delete: " . Data::Dumper::Dumper($_);
        }
    }

    #creating an anonymous sub to delete allocations when commit happens
    my $id = $self->id;
    my $upon_delete_callback = sub {
        print "Now Deleting Allocation with owner_id = $id\n";
        my $allocation = Genome::Disk::Allocation->get(owner_id=>$id, owner_class_name=>$class_name);
        if ($allocation) {
            $allocation->deallocate;
        }
    };

    #hook our anonymous sub into the commit callback
    Genome::Sys::CommitAction->create( on_commit => $upon_delete_callback );

    return $self->SUPER::delete(@_);
}

sub _lock {
    my $class = shift;
    my $lookup_hash = _validate_lookup_hash(shift);
    my $wait = shift || 1;

    my $resource_lock_name = $class->_resolve_lock_name($lookup_hash);

    # Maybe we should memoize on scope and resource but all SoftwareResult are
    # in the 'site' scope which isn't currently expected to change.
    my $lock = $LOCKS{$resource_lock_name}{object};
    my $count = $LOCKS{$resource_lock_name}{count};
    return $lock if $count;

    $lock = Genome::Sys::LockProxy->new(
        resource => $resource_lock_name,
        scope => 'site',
    )->lock(
        max_try => 2,
    );
    if (!$lock && $wait) {
        if (Genome::Config::get("software_result_async_locking")) {
            my $exit_code = Genome::Config::get(
                "software_result_async_locking_exit_code");
            $class->status_message("SoftwareResult asynchronous locking " .
                "is enabled... exiting now with exit_code (%s).", $exit_code);
            exit $exit_code;
        } else {
            $class->debug_message("This data set is still being processed by " .
                "its creator.  Waiting for existing data lock...");
            $lock = Genome::Sys::LockProxy->new(
                resource => $resource_lock_name,
                scope => 'site',
            )->lock(
                wait_announce_interval => 600,
            );
            unless ($lock) {
                $class->error_message("Failed to get existing data lock!");
                die($class->error_message);
            }
        }
    }

    $LOCKS{$resource_lock_name}{object} = $lock;
    $LOCKS{$resource_lock_name}{count}++;

    return $lock;
}

sub _unlock_resource {
    my $class = shift;
    my $resource_lock = shift;

    my $resource_lock_name = $resource_lock->resource;

    $class->debug_message("Cleaning up lock $resource_lock_name...");

    if (!exists $LOCKS{$resource_lock_name}{count})  {
        $class->error_message("Attempt to unlock $resource_lock_name but this was never locked!");
        die $class->error_message;
    }

    $LOCKS{$resource_lock_name}{count}--;

    return 1 if ($LOCKS{$resource_lock_name}{count} >= 1);

    unless ($resource_lock->unlock()) {
        $class->error_message("Couldn't unlock $resource_lock_name.  error message was " . $class->error_message);
        die $class->error_message;
    }

    delete $LOCKS{$resource_lock_name};

    $class->debug_message("Cleanup completed for lock $resource_lock_name.");
    return 1;
}

sub _unlock {
    my $self = shift;

    my $resource_lock = $self->_lock_proxy;
    $self->_unlock_resource($resource_lock);
}

sub _resolve_lock_name {
    my $class = shift;
    my $lookup_hash = _validate_lookup_hash(shift);

    my $class_string = $class->class;
    return "genome/$class_string/" .  $lookup_hash;
}

sub _validate_lookup_hash {
    my $lookup_hash = shift;
    unless ($lookup_hash =~ /^[[:xdigit:]]{32}$/) {
        croak "invalid lookup_hash: $lookup_hash";
    }
    return $lookup_hash;
}

# override _resolve_lock_name (for testing) to append username and time
# This override is used to prevent lock collisions when tests are being run concurrently on the same machine.
if ($ENV{UR_DBI_NO_COMMIT}) {
    warn 'Overriding Genome::SoftwareResult::_resolve_lock_name since UR_DBI_NO_COMMIT is on.' . "\n";
    my $suffix = Genome::Sys->username . '_' . time;
    my $original_resolve_lock_name_sub = \&Genome::SoftwareResult::_resolve_lock_name;
    require Sub::Install;
    Sub::Install::reinstall_sub({
        into => 'Genome::SoftwareResult',
        as => '_resolve_lock_name',
        code => sub {
            my $lock_name = &$original_resolve_lock_name_sub(@_);
            $lock_name .= "_$suffix" unless $lock_name =~ /$suffix/;
            return $lock_name;
        },
    });
}

sub metric_names {
    my $class = shift;
    my $meta = $class->__meta__;
    my @properties = grep { $_->{is_metric} } $meta->_legacy_properties();
    my @names = map { $_->property_name } @properties;
    return @names;
}

sub metrics_hash {
    my $self = shift;
    my @names = $self->metric_names;
    my %hash = map { $self->name } @names;
    return %hash;
}

sub generate_expected_metrics {
    my $self = shift;
    my @names = @_;
    unless (@names) {
        @names = $self->metric_names;
    }

    # pre-load all metrics
    my @existing_metrics = $self->metrics;

    for my $name (@names) {
        my $metric = $self->metric(name => $name);
        if ($metric) {
            $self->debug_message(
                $self->display_name . " has metric "
                . $metric->name
                . " with value "
                . $metric->value
            );
            next;
        }
        my $method = "_calculate_$name";
        unless ($self->can($method)) {
            $self->error_message("No method $method found!");
            die $self->error_message;
        }
        $self->debug_message(
            $self->display_name . " is generating a value for metric "
            . $metric->name
            . "..."
        );
        my $value = $self->$method();
        unless (defined($value)) {
            $self->error_message(
                $self->display_name . " has metric "
                . $metric->name
                . " FAILED TO CALCULATE A DEFINED VALUE"
            );
            next;
        }
        $self->$metric($value);
        $self->debug_message(
            $self->display_name . " has metric "
            . $metric->name
            . " with value "
            . $metric->value
        );
    }
}

sub _available_cpu_count {
    my $self = shift;

    # Not running on LSF, allow only one CPU
    if (!exists $ENV{LSB_MCPU_HOSTS}) {
        return 1;
    }

    my $mval = $ENV{LSB_MCPU_HOSTS};
    my @c = split /\s+/, $mval;

    if (scalar @c != 2) {
        $self->error_message("LSB_MCPU_HOSTS environment variable doesn't specify just one host and one CPU count. (value is '$mval').  Is the span[hosts=1] value set in your resource request?");
        die $self->error_message;
    }

    if ($mval =~ m/(\.*?) (\d+)/) {
        return $2;
    } else {
        $self->error_message("Couldn't parse the LSB_MCPU_HOSTS environment variable (value is '$mval'). ");
        die $self->error_message;
    }

}

sub _resolve_param_value_from_text_by_name_or_id {
    my $class = shift;
    my $param_arg = shift;

    #First try default behaviour of looking up by name or id
    my @results = Command::V2->_resolve_param_value_from_text_by_name_or_id($class, $param_arg);

    #If that didn't work, and the argument is a filename, see if it's part of our output directory.
    if(!@results and -f $param_arg) {
        my $abs_path = Genome::Sys->abs_path($param_arg);
        my (undef, $dir) = fileparse($abs_path);
        $dir =~ s!/$!!; #remove trailing slash!
        @results = Genome::SoftwareResult->get(output_dir => $dir);
    }

    return @results;
}

sub _release_lock_or_die {
    my ($class, $lock, $error_message) = @_;

    $class->debug_message("Cleaning up lock %s...", $lock->resource);

    unless ($lock->unlock()) {
        $class->error_message($error_message);
        die $error_message;
    }
    delete $LOCKS{$lock->resource};

    $class->debug_message("Cleanup completed for lock %s.", $lock->resource);
}

# children are things that use this
sub children {
    my $self = shift;
    return uniq map { $_->user } $self->users;
}

# parents are things that this uses
sub parents {
    my $self = shift;
    return uniq map { $_->software_result }
        Genome::SoftwareResult::User->get(user => $self);
}

# ancestors are recursive parents, sort order is not guaranteed
sub ancestors {
    my $self = shift;
    my @parents = $self->parents;
    if (@parents) {
        return uniq(map { $_->ancestors }
            grep { $_->isa('Genome::SoftwareResult') } @parents), @parents;
    } else {
        return;
    }
}

# descendents are recursive children, sort order is not guaranteed
sub descendents {
    my $self = shift;
    my @children = $self->children;
    if (@children) {
        return @children, uniq map { $_->descendents }
            grep { $_->isa('Genome::SoftwareResult') } @children;
    } else {
        return;
    }
}

sub best_guess_date {
    my $self = shift;
    my ($earliest_time) = sort map { $_->creation_time }
        $self->disk_allocations;
    return $earliest_time;
}

sub best_guess_date_numeric {
    return UnixDate(shift->best_guess_date, "%s");
}

sub expunge_results_containing_object {
    my $class = shift;
    my $object = shift;
    my $reason = shift;

    my @params_and_inputs = (Genome::SoftwareResult::Param->get(value_id => $object->id), Genome::SoftwareResult::Input->get(value_id => $object->id));
    my @results_to_expunge = Genome::SoftwareResult->get([map $_->software_result_id, @params_and_inputs]);

    map $_->expunge($reason), @results_to_expunge;
}

sub expunge {
    my $self = shift;
    my $reason = shift;

    for my $child (grep { $_->isa('Genome::SoftwareResult') } $self->children) {
        $child->expunge($reason);
    }

    if($self->disk_allocation) {
        eval {
            $self->disk_allocation->purge($reason);
        };
        if(my $error = $@) {
            $self->warning_message($@);
        }
    }

    #disk allocation activity may have updated object so set this afterwards.
    $self->test_name($reason);

    return 1;
}

1;
