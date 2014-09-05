package Genome::SoftwareResult;

use strict;
use warnings;

use Genome;
use Digest::MD5 qw(md5_hex);
use Cwd;
use File::Basename qw(fileparse);
use Data::Dumper;
use Date::Manip;
use List::MoreUtils qw(uniq);
use File::Grep qw(fgrep);

use Carp;

use JSON;

use Genome::Utility::Instrumentation;

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
        #inputs_bx           => { is => 'UR::BoolExpr', id_by => 'inputs_id', is_optional => 1 },
        inputs_id           => { is => 'Text', len => 4000, column_name => 'INPUTS_ID', implied_by => 'inputs_bx', is_optional => 1 },
        #params_bx           => { is => 'UR::BoolExpr', id_by => 'params_id', is_optional => 1 },
        params_id           => { is => 'Text', len => 4000, column_name => 'PARAMS_ID', implied_by => 'params_bx', is_optional => 1 },
        output_dir          => { is => 'Text', len => 1000, column_name => 'OUTPUTS_PATH', is_optional => 1 },
        test_name           => { is_param => 1, is_delegated => 1, is_mutable => 1, via => 'params', to => 'value_id', where => ['name' => 'test_name'], is => 'Text', doc => 'Assigns a testing tag to the result.  These will not be used in default processing', is_optional => 1 },
        _lock_name          => { is_optional => 1, is_transient => 1 },
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

Genome::SoftwareResult->add_observer(
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

    my $statsd_prefix = "software_result_get.";
    my $statsd_class_suffix = "$class";
    $statsd_class_suffix =~ s/::/_/g;

    my $start_time = Time::HiRes::time();

    my $lookup_hash = $class->calculate_lookup_hash_from_arguments(@_);

    # NOTE we do this so that get can be noisy when called directly.
    my @objects = $class->SUPER::get(lookup_hash => $lookup_hash);

    my $final_time = Time::HiRes::time();
    my $full_time = 1000 * ($final_time - $start_time);
    Genome::Utility::Instrumentation::timing($statsd_prefix . "full_time.total",
            $full_time);
    Genome::Utility::Instrumentation::timing(
        $statsd_prefix . "full_time." . $statsd_class_suffix, $full_time);

    return @objects;
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

    return if (@_ > 1);

    if (@_ == 1) {
        if (Scalar::Util::blessed($_[0])) {
            if ($_[0]->isa('UR::BoolExpr') && $_[0]->is_id_only) {
                return 1;
            }
            return;
        }

        return 1;
    }
    return;
}

# You must specify enough parameters to uniquely identify an object to get a result.
# If two users specify different sets of parameters that uniquely identify the result,
# they will create different locks.
sub get_with_lock {
    my $class = shift;

    my $params_processed = $class->_gather_params_for_get_or_create($class->_preprocess_params_for_get_or_create(@_));

    my %is_input = %{$params_processed->{inputs}};
    my %is_param = %{$params_processed->{params}};

    # Only try with lock if object does not exist since locking causes
    # a performance hit. It is assumed that if an object is found it is
    # complete. If this is a bad assumption then we need to add a
    # status to SoftwareResults.
    my $lock;
    my @objects = $class->_faster_get(@_);
    unless (@objects) {
        my $subclass = $params_processed->{subclass};

        my $lookup_hash = $subclass->calculate_lookup_hash_from_arguments(@_);
        unless ($lock = $subclass->_lock($lookup_hash)) {
            die "Failed to get a lock for " . Dumper(@_);
        }

        UR::Context->current->reload($class,
            lookup_hash => $class->calculate_lookup_hash_from_arguments(@_));

        eval {
            @objects = $class->_faster_get(@_);
        };
        my $error = $@;

        if ($error) {
            $class->error_message('Failed in get! ' . $error);
            $class->_release_lock_or_die($lock, "Failed to unlock during get_with_lock.");
            die $class->error_message;
        }
    }

    if (@objects > 1) {
        $class->error_message("Multiple results returned for SoftwareResult ${class}::get_with_lock.  To avoid this, call get_with_lock with enough parameters to uniquely identify a SoftwareResult.");
        $class->error_message("Parameters used for the get: " . Data::Dumper::Dumper %is_input . Data::Dumper::Dumper %is_param);
        $class->error_message("Objects gotten: " . Data::Dumper::Dumper @objects);
        $class->_release_lock_or_die($lock, "Failed to unlock during get_with_lock with multiple results.") if $lock;
        die $class->error_message;
    }

    my $result = $objects[0];

    if ($result) {
        my $calculated_lookup_hash = $result->calculate_lookup_hash();
        my $lookup_hash = $result->lookup_hash;
        if ($calculated_lookup_hash ne $lookup_hash) {
            my $m = sprintf(q{SoftwareResult lookup_hash (%s) does not match it's calculated lookup_hash (%s).  Cannot trust that the correct SoftwareResult was retrieved.}, $lookup_hash, $calculated_lookup_hash);
            # Really we would just call get(%$params) but that might undermine the whole lookup_hash optimization.
            die $class->error_message($m);
        }
    }

    if ($result && $lock) {
        $result->_lock_name($lock);

        $result->debug_message("Cleaning up lock $lock...");
        unless ($result->_unlock) {
            $result->error_message("Failed to unlock after getting software result");
            die "Failed to unlock after getting software result";
        }
        $result->debug_message("Cleanup completed for lock $lock.");
    } elsif ($lock) {
        $class->_release_lock_or_die($lock, "Failed to unlock after not finding software result.");
    }

    $result->_auto_unarchive if $result;
    return $result;
}

sub get_or_create {
    my $class = shift;

    my @objects = $class->_faster_get(@_);

    unless (@objects) {
        @objects = $class->create(@_);
        unless (@objects) {
            # see if the reason we failed was b/c the objects were created while we were locking...
            @objects = $class->_faster_get(@_);
            unless (@objects) {
                $class->error_message("Could not create a $class for params " . Data::Dumper::Dumper(\@_) . " even after trying!");
                confess $class->error_message();
            }
        }
    }

    if (@objects > 1) {
        if(wantarray) {
            map $_->_auto_unarchive, @objects;
            return @objects;
        }
        my @ids = map { $_->id } @objects;
        die "Multiple matches for $class but get or create was called in scalar context!  Found ids: @ids";
    } else {
        $objects[0]->_auto_unarchive if $objects[0];
        return $objects[0];
    }
}

sub create {
    my $class = shift;

    if ($class eq __PACKAGE__ || $class->__meta__->is_abstract) {
        # this class is abstract, and the super-class re-calls the constructor from the correct subclass
        return $class->SUPER::create(@_);
    }

    my @previously_existing = $class->_faster_get(@_);
    if (@previously_existing > 0) {
        $class->error_message("Attempt to create an $class but it looks like we already have one with those params " . Dumper(\@_));
        return;
    }

    my $lookup_hash = $class->calculate_lookup_hash_from_arguments(@_);
    my $lock = $class->_get_lock_or_die($lookup_hash);

    # we might have had to wait on the lock, in which case someone else was probably creating that entity
    # do a "reload" here to force another trip back to the database to see if a software result was created
    # while we were waiting on the lock.
    (@previously_existing) = UR::Context->current->reload($class,
        lookup_hash => $lookup_hash);

    if (@previously_existing > 0) {
        $class->error_message("Attempt to create an $class but it looks like we already have one with those params " . Dumper(\@_));
        $class->_release_lock_or_die($lock, "Failed to release lock in create before committing SoftwareResult.");
        return;
    }

    my $bx = $class->_get_safe_boolexpr(@_);
    my $self = $class->SUPER::create($bx);
    unless ($self) {
        $class->_release_lock_or_die($lock,"Failed to unlock during create after committing SoftwareResult.");
        return;
    }

    $self->_lock_name($lock);
    $self->_set_unlock_callbacks();

    if ($self->output_dir) {
        if (-d $self->output_dir) {
            $self->_validate_output_dir_is_empty;
        } else {
            $self->_create_output_dir_or_die;
        }
    }

    $self->module_version($self->resolve_module_version) unless defined $self->module_version;
    $self->subclass_name($class);
    $self->lookup_hash($lookup_hash);
    return $self;
}

sub _get_lock_or_die {
    my $class = shift;
    my $lookup_hash = shift;

    my $lock;
    if (my $lock = $class->_lock($lookup_hash)) {
        return $lock;
    } else {
        die "Failed to get a lock for class ($class) and lookup_hash ($lookup_hash)";
    }
}

# TODO update the indirect mutable accessor logic for non-nullable
# hang-offs to delete the entry instead of setting it to null.
#
# Otherwise we get SOFTWARE_RESULT_PARAM entries with a NULL, and unsavable
# PARAM_VALUE.  also remove empty strings because that's equivalent to a
# NULL to the database
#
# This function generates a boolexpr that filters out such params and inputs.
sub _get_safe_boolexpr {
    my $class = shift;

    my $params_processed = $class->_gather_params_for_get_or_create($class->_preprocess_params_for_get_or_create(@_));

    my %is_input = %{$params_processed->{inputs}};
    my %is_param = %{$params_processed->{params}};

    my @param_remove = grep { not (defined $is_param{$_}) || $is_param{$_} eq "" } keys %is_param;

    # Do the same for inputs (e.g. alignment results have nullable segment values for instrument data, which are treated as inputs)
    my @input_remove = grep { not (defined $is_input{$_}) || $is_input{$_} eq "" } keys %is_input;

    my $bx = $class->define_boolexpr($class->_preprocess_params_for_get_or_create(@_));
    for my $i (@param_remove, @input_remove) {
        $bx = $bx->remove_filter($i);
    }
    return $bx;
}

# Ensure that we release the lock we have obtained if the SR object gets
# committed or deleted.
sub _set_unlock_callbacks {
    my $self = shift;

    my $unlock_callback = sub {
        $self->_unlock;
    };
    $self->create_subscription(method=>'commit', callback=>$unlock_callback);
    $self->create_subscription(method=>'delete', callback=>$unlock_callback);
    return;
}


sub _validate_output_dir_is_empty {
    my $self = shift;

    my $output_dir = $self->output_dir;
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

sub _create_output_dir_or_die {
    my $self = shift;

    $self->debug_message("Creating output directory (%s)...", $self->output_dir);
    eval {
        Genome::Sys->create_directory($self->output_dir)
    };
    if ($@) {
        $self->delete;
        die $@;
    }
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

    #my $inputs_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_input);
    #my $params_bx = UR::BoolExpr->resolve_normalized_rule_for_class_and_params($class, %is_param);

    my %software_result_params = (#software_version=>$params_bx->value_for('aligner_version'),
        #params_id=>$params_bx->id,
        #inputs_id=>$inputs_bx->id,
        subclass_name=>$class
    );

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
    local $@;
    $object->software_result->recalculate_lookup_hash();
};
for my $name (qw(Param Input)) {
    my $classname = join('::', qw(Genome SoftwareResult), $name);
    for my $aspect (qw(create value_id)) {
        $classname->add_observer(aspect => $aspect, callback => $recalculate_lookup_hash_callback);
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

    my $build = $ENV{GENOME_BUILD_ID};
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
        while (my ($prop_name, $prop_desc) = each(%{ $desc->{has} })) {
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
    my $observer;
    my $upon_delete_callback = sub {
        print "Now Deleting Allocation with owner_id = $id\n";
        $observer->delete if $observer;
        my $allocation = Genome::Disk::Allocation->get(owner_id=>$id, owner_class_name=>$class_name);
        if ($allocation) {
            $allocation->deallocate;
        }
    };

    #hook our anonymous sub into the commit callback
    $observer = $class_name->ghost_class->add_observer(aspect=>'commit', callback=>$upon_delete_callback);

    return $self->SUPER::delete(@_);
}

sub _lock {
    my $class = shift;
    my $lookup_hash = _validate_lookup_hash(shift);
    my $wait = shift || 1;

    my $resource_lock_name = $class->_resolve_lock_name($lookup_hash);

    # if we're already locked, just increment the lock count
    $LOCKS{$resource_lock_name} += 1;
    return $resource_lock_name if ($LOCKS{$resource_lock_name} > 1);

    my $lock = Genome::Sys->lock_resource(resource_lock => $resource_lock_name, max_try => 2);
    if (!$lock && $wait) {
        $class->debug_message("This data set is still being processed by its creator.  Waiting for existing data lock...");
        $lock = Genome::Sys->lock_resource(resource_lock => $resource_lock_name, wait_announce_interval => 600);
        unless ($lock) {
            $class->error_message("Failed to get existing data lock!");
            die($class->error_message);
        }
    }

    return $lock;
}

sub _unlock_resource {
    my $class = shift;
    my $resource_lock_name = shift;

    $class->debug_message("Cleaning up lock $resource_lock_name...");

    if (!exists $LOCKS{$resource_lock_name})  {
        $class->error_message("Attempt to unlock $resource_lock_name but this was never locked!");
        die $class->error_message;
    }
    $LOCKS{$resource_lock_name} -= 1;

    return 1 if ($LOCKS{$resource_lock_name} >= 1);

    unless (Genome::Sys->unlock_resource(resource_lock=>$resource_lock_name)) {
        $class->error_message("Couldn't unlock $resource_lock_name.  error message was " . $class->error_message);
        die $class->error_message;
    }

    delete $LOCKS{$resource_lock_name};
    $class->debug_message("Cleanup completed for lock $resource_lock_name.");
    return 1;
}

sub _unlock {
    my $self = shift;

    my $resource_lock_name = $self->_lock_name;
    $self->_unlock_resource($resource_lock_name);
}

sub _resolve_lock_name {
    my $class = shift;
    my $lookup_hash = _validate_lookup_hash(shift);

    my $class_string = $class->class;
    return $ENV{GENOME_LOCK_DIR} . "/genome/$class_string/" .  $lookup_hash;
}

sub _validate_lookup_hash {
    my $lookup_hash = shift;
    unless (length($lookup_hash) == 32 && $lookup_hash =~ /^[\d\w]+$/) {
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
        my $abs_path = Cwd::abs_path($param_arg);
        my (undef, $dir) = fileparse($abs_path);
        $dir =~ s!/$!!; #remove trailing slash!
        @results = Genome::SoftwareResult->get(output_dir => $dir);
    }

    return @results;
}

sub _release_lock_or_die {
    my ($class, $lock, $error_message) = @_;

    $class->debug_message("Cleaning up lock $lock...");

    unless (Genome::Sys->unlock_resource(resource_lock=>$lock)) {
        $class->error_message($error_message);
        die $error_message;
    }
    delete $LOCKS{$lock};

    $class->debug_message("Cleanup completed for lock $lock.");
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

sub creation_grid_job {
    my $self = shift;
    my ($creation_build, $lsf_job_id) = $self->creation_build_and_lsf_job_id;
    unless ($creation_build && $lsf_job_id) {
        $self->error_message('Failed to find a grid job since the creation build and/or event is lost for software result: '. $self->id);
        return;
    }
    my $project_name = 'build' . $creation_build->id;
    my $grid_job = Genome::Site::TGI::GridJobsFinished->get(
        project_name => $project_name,
        job_id => $lsf_job_id,
    );
    unless ($grid_job) {
        $self->error_message('Failed to find grid job with project '. $project_name .' and job_id '. $lsf_job_id);
        return;
    }
    return $grid_job;
}

sub creation_build_and_lsf_job_id {
    my $self = shift;

    # Find the disk allocation ID to look for in the LSF error logs
    my @allocations = $self->disk_allocations;
    if (scalar(@allocations) != 1) {
        $self->error_message('Please add logic to handle multiple allocations for a software result!');
        return;
    }
    my $allocation = shift @allocations;
    my $allocation_id = $allocation->id;

    # Generate an ordered list of builds, search the oldest builds first
    my @ordered_builds = $self->builds(-order_by => 'date_scheduled');
    unless ( scalar(@ordered_builds) ) {
        # There are no direct links to builds, walk through 1st generation children to find builds
        my @children = $self->children;
        for my $child (@children) {
            # TODO: Should reorder after this push, but this will work as a first pass
            push @ordered_builds, $child->builds(-order_by => 'date_scheduled');
        }
        unless ( scalar(@ordered_builds) ) {
            $self->error_message('Failed to find any builds for software result: '. $self->id);
            return;
        }
    }
    # The output directory of the software result should be in the LSF error log
    my $output_dir = $self->output_dir;
    $self->debug_message('For software result '. $self->id .', searching for LSF error log containing output directory: '. $output_dir);

    # Iterate over every build linked to this software result and grep all event error logs
    my $creation_build;
    my $creation_lsf_job_id;
    for ( my $i=0; $i < scalar(@ordered_builds); $i++ ) {
        my $build = $ordered_builds[$i];
        $self->debug_message('Evaluating build '. $build->id .' as the potential creation build.');
        my @events = $build->events;
        for my $event (@events) {
            my $error_log_file = $event->error_log_file();
            $self->debug_message('Grep through error log: '. $error_log_file);
            # This requires that one line in the error log ends with the output_dir
            if (fgrep { /Allocation \($allocation_id\) created at $output_dir$/ } $error_log_file) {
                $self->debug_message('Found disk allocation id '. $allocation_id .' and output_dir '. $output_dir .' in error log '. $error_log_file);
                unless ($creation_build) {
                    $creation_build = $build;
                } else {
                    $self->error_message('Found multiple potential creation builds: '. $creation_build->id .' and '. $build->id);
                    return;
                }
                unless ($creation_lsf_job_id) {
                    $creation_lsf_job_id = $event->lsf_job_id;
                } else {
                    $self->error_message('Found multiple potential creation lsf job ids: '. $creation_lsf_job_id .' and '. $event->lsf_job_id);
                    return;
                }
            }
        }
        unless ($creation_build && $creation_lsf_job_id) {
            my @instances = ($build->newest_workflow_instance, $build->child_workflow_instances);
            my @ie = map { $_->current } @instances;
            for my $ie (@ie) {
                my $error_log_file = $ie->stderr;
                # TODO: Remove redundany with above, refactor to subroutine? or get common event/instance interface
                # This requires that one line in the error log ends with the output_dir
                if (fgrep { /Allocation \($allocation_id\) created at $output_dir$/ } $error_log_file) {
                    $self->status_message('Found disk allocation id '. $allocation_id .' and output_dir '. $output_dir .' in error log '. $error_log_file);
                    unless ($creation_build) {
                        $creation_build = $build;
                    } else {
                        $self->error_message('Found multiple potential creation builds: '. $creation_build->id .' and '. $build->id);
                        return;
                    }
                    unless ($creation_lsf_job_id) {
                        $creation_lsf_job_id = $ie->dispatch_identifier;
                    } else {
                        $self->error_message('Found multiple potential creation events: '. $creation_lsf_job_id .' and '. $ie->dispatch_identifier);
                        return;
                    }
                }
            }
        }
    }
    # Return the build and LSF job id most likely to originally have created this software result
    return ($creation_build, $creation_lsf_job_id);
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

    $self->test_name($reason);

    if($self->disk_allocation) {
        eval {
            $self->disk_allocation->purge($reason);
        };
        if(my $error = $@) {
            $self->warning_message($@);
        }
    }

    return 1;
}

1;
