package Genome::Model::Event;

use strict;
use warnings;
use File::Path;
use YAML;

use Genome;
class Genome::Model::Event {
    is => 'Genome::Model::Command::BaseDeprecated',
    table_name => 'model.event',
    is_abstract => 1,
    first_sub_classification_method_name => '_resolve_subclass_name',
    sub_classification_method_name => '_resolve_subclass_name',
    subclass_description_preprocessor => 'Genome::Model::Event::_preprocess_subclass_description',
    type_name => 'genome model event',
    id_generator => '-uuid',
    id_by => [
        genome_model_event_id => {
            is => 'Text',
            len => 32,
        },
    ],
    has => [
        model => {
            is => 'Genome::Model',
            id_by => 'model_id',
            constraint_name => 'GME_GM_FK',
        },
        event_type => {
            is => 'VARCHAR2',
            len => 255,
        },
        event_status => {
            is => 'VARCHAR2',
            len => 32,
            is_optional => 1,
        },
        user_name => {
            is => 'VARCHAR2',
            len => 64,
            is_optional => 1,
        },
        build_id => {
            is => 'NUMBER',
            is_optional => 1,
        },
        run_id => {
            is => 'NUMBER',
            len => 11,
            is_optional => 1,
        },
        workflow_instance_id => {
            is => 'Number',
            len => 11,
            is_optional => 1,
        },
    ],
    has_optional => [
        instrument_data_id => {
            is => 'VARCHAR2',
            len => 100,
        },
        instrument_data => {
            is => 'Genome::InstrumentData',
            id_by => 'instrument_data_id',
            doc => 'The id of the instrument data on which to operate',
        },
        instrument_data_input => {
            calculate_from => [ 'instrument_data_id', 'model' ],
            calculate => q(
                return $model->input_for_instrument_data_id($instrument_data_id);
            ),
        },
        ref_seq_id => {
            is => 'VARCHAR2',
            len => 64,
        },
        parent_event => {
            is => 'Genome::Model::Event',
            id_by => 'parent_event_id',
            constraint_name => 'GME_PAEID_FK',
        },
        prior_event => {
            is => 'Genome::Model::Event',
            id_by => 'prior_event_id',
            constraint_name => 'GME_PPEID_FK',
        },
        date_completed => { is => 'TIMESTAMP' },
        date_scheduled => { is => 'TIMESTAMP' },
        lsf_job_id => {
            is => 'VARCHAR2',
            len => 64,
        },
        retry_count => {
            is => 'NUMBER',
            len => 3,
        },
        status_detail => {
            is => 'VARCHAR2',
            len => 200,
        },
        build => {
            is => 'Genome::Model::Build',
            id_by => 'build_id',
        },
        should_calculate => {
            calculate_from => 'event_status',
            calculate => q(
                                 if ($event_status eq 'Failed' or $event_status eq 'Crashed') {
                                     return 0;
                                 }
                                 return 1;
                             ),
            doc => 'a flag to determine metric calculations',
        },
        build_directory => {
            calculate_from => 'build',
            calculate => q( return $build->data_directory ),
            doc => 'the directory where this step should put data',
        },
    ],
    has_many_optional => [
        sibling_events => {
            via => 'parent_event',
            to => 'child_events',
            is_many => 1,
        },
        child_events => {
            is => 'Genome::Model::Event',
            reverse_as => 'parent_event',
        },
        inputs => {
            is => 'Genome::Model::Event::Input',
            reverse_as => 'event',
        },
        outputs => {
            is => 'Genome::Model::Event::Output',
            reverse_as => 'event',
        },
        metrics => {
            is => 'Genome::Model::Event::Metric',
            reverse_as => 'event',
        },
        metric_names => {
            via => 'metrics',
            to => 'name',
        },
        instrument_data_segment_id_param => {
            is => 'Genome::Model::Event::Input',
            reverse_as => 'event',
            where => [ name => 'instrument_data_segment_id' ],
        },
        instrument_data_segment_type_param => {
            is => 'Genome::Model::Event::Input',
            reverse_as => 'event',
            where => [ name => 'instrument_data_segment_type' ],
        },
        instrument_data_segment_id => {
            via => 'instrument_data_segment_id_param',
            to => 'value',
        },
        instrument_data_segment_type => {
            via => 'instrument_data_segment_type_param',
            to => 'value',
        },
    ],
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',
};

sub command_name {
    my $self = shift;
    my $command_name = $self->SUPER::command_name(@_);
    $command_name =~ s/event // unless $command_name eq 'genome model event generic';
    return $command_name;
}

sub _preprocess_subclass_description {
    my ($class,$desc) = @_;
    my $has = $desc->{has};
    for my $attribute (values %$has) {
        #print "event has $attribute " . Data::Dumper::Dumper($attribute) ,"\n";
        # via => 'inputs', to => 'value', where => [ name => 'foobar' ]
        if ($attribute->{is_input}) {
            delete($attribute->{is_input});
            $attribute->{via} = 'inputs';
            $attribute->{to} = 'value';
            $attribute->{where} = [ name => $attribute->{property_name} ];
        } elsif ($attribute->{is_output}) {
            delete($attribute->{is_output});
            $attribute->{via} = 'outputs';
            $attribute->{to} = 'value';
            $attribute->{where} = [ name => $attribute->{property_name} ];
        } elsif ($attribute->{is_param}) {
            die($attribute->{property_name} .' is a param and parameters are not handled correctly by events');
        }
        $has->{$attribute->{property_name}} = $attribute;
    }
    $desc->{has} = $has;
    return $desc;
}

sub create {
    my $class = shift;
    my $self = $class->SUPER::create(@_);
    return unless $self;
    unless ($self->event_type) {
        $self->event_type($self->command_name);
    }
    unless ($self->date_scheduled) {
        $self->date_scheduled(UR::Context->current->now);
    }
    unless ($self->user_name) {
        my $user = getpwuid($<);
        $self->user_name($user);
    }
    return $self;
}

sub revert {
    my $self = shift;
    # don't roll back inputs.
    for my $obj (grep {! $_->isa('Genome::Model::Event::Input') } $self->get_all_objects) {
        if ($obj->isa('Genome::Model::Event')) {
            if ($obj->parent_event_id eq $self->id) {
                # Remove foreign keys
                $obj->parent_event_id(undef);
                #UR::Context->_sync_databases;
            }
        }
        $self->warning_message('deleting '. $obj->class .' with id '. $obj->id);
         #this delete is a general UR delete not the cute delete just below - Jim & Chris
        $obj->delete;
    }
    return 1;
}

sub get_all_objects {
    my $self = shift;
    my @inputs = $self->inputs;
    my @outputs = $self->outputs;
    my @metrics = $self->metrics;
    return sort {$a->id cmp $b->id} (@inputs, @outputs, @metrics);
}

sub yaml_string {
    my $self = shift;
    my $string = YAML::Dump($self);
    for my $object ($self->get_all_objects) {
        $string .= YAML::Dump($object);
    }
    return $string;
}

sub delete {
    my $self = shift;

    # get rid of associated metrics
    my @metrics = Genome::Model::Event::Metric->get(event_id=>$self->id);
    for (@metrics) {
        $_->delete;
    }

    my @inputs = Genome::Model::Event::Input->get(event_id=>$self->id);
    for (@inputs) {
        $_->delete;
    }

    #removing constraints
    my @next_events = $self->next_events;
    for my $next_event (@next_events) {
        $next_event->prior_event_id(undef);
        #UR::Context->_sync_databases;
        print Data::Dumper::Dumper($next_event);
        $self->warning_message('Setting event '.$next_event->id.' to undef to remove prior_event_id constraint.');
    }

    if ($self->{db_committed}) {
        $self->warning_message('deleting ' . $self->class . ' with id  ' . $self->id);
    }
    $self->revert;
    $self->SUPER::delete();
    return 1;
}

#TODO: fix the way this method is called an either inherit from G:U:F
#      or call class method directly like below
sub shellcmd {
    my $self = shift;
    return Genome::Sys->shellcmd(@_);
}

#< Logging >#
sub resolve_log_directory { return log_directory(@_); }
sub log_directory {
    my $self = shift;
    my $log_directory = sprintf('%s/logs/',
                                $self->build_directory,
                            );
    if (defined $self->ref_seq_id) {
        $log_directory .= $self->ref_seq_id;
    }
    return $log_directory
}

sub create_log_directory {
    my $self = shift;

    my $log_dir = $self->log_directory;
    Genome::Sys->create_directory($log_dir)
        or die "Can't create log directory ($log_dir) for event: ".$self->id;

    return 1
}

sub _log_file {
    my ($self, $ext) = @_;

    unless ( $self->build ) {
        $self->error_message("Can't get log file because there is no build.");
        return;
    }

    my $old = sprintf(
        '%s/%s.%s',
        $self->log_directory,
        $self->genome_model_event_id,
        $ext,
    );

    if (-e $old) {
        return $old;
    }
    else {
        my $new = $self->event_type;
        my $log_directory = $self->log_directory;
        $log_directory =~ s|/$||;
        $new =~ s/^genome.model.//;
        $new =~ s/[ -]/_/g;
        if ($new eq 'build') {
            $new = 'workflow-server';
        }
        else {
            $new =~ s/^build//;
            $new .= '.' . $self->id;
            $new =~ s/^_//g;
        }
        $new = $log_directory . '/' . $new . '.' . $ext;
        return $new;
    }
}

sub error_log_file {
    return _log_file(@_, 'err');
}

sub output_log_file {
    return _log_file(@_, 'out');
}
#<>#

sub check_for_existence {
    my ($self,$path,$attempts) = @_;

    unless (defined $attempts) {
        $attempts = 5;
    }

    my $try = 0;
    my $found = 0;
    while (!$found && $try < $attempts) {
        $found = -e $path;
        sleep(1);
        $try++;
        if ($found) {
            $self->status_message("existence check passed: $path");
            return $found;
        }
    }
    return;
}

sub create_file {
    my ($self, $output_name, $path) = @_;
    if (!$path) {
        die "Output $output_name opened without a specified path!"
    }
    elsif (my @existing = $self->outputs(name => $output_name)) {
        if ($output_name and $existing[0]->value ne $path) {
            die "Input $output_name already exists with value " . $existing[0]->value
            . ".  Cannot open with value $path";
        }
        die "Attempting to re-create already created file! $output_name: $path";
    }
    elsif (-e $path) {
        die "File $path already exists!  Canot create output $output_name: $path\n";
    }
    else {
        my $fh = IO::File->new('>'.$path);
        die "Failed to make file $path! $?" unless $fh;
        unless ($self->add_output(name => $output_name, value => $path)) {
            die "Error adding ouput $output_name $path!";
        }

    return $fh;
    }
}

sub open_file {
    my ($self, $input_name, $path) = @_;
    my @existing = $self->inputs(name => $input_name);
    if (@existing and $input_name and $existing[0]->value ne $input_name) {
        die "Input $input_name already exists with value " . $existing[0]->value
            . ".  Cannot open with value $path";
    }
    if (not @existing) {
        if (!$path) {
            die "Input $input_name opened without a specified path, and no path has been set yet!"
        }
        else {
            $self->add_input(name => $input_name, value => $path);
        }
    }
    my $fh = IO::File->new($path);
    die "Failed to open file $path! $?" unless $fh;
    if ( (!-p $path) and (-z $path) ) {
        warn "warning: opening zero-length file $path for input $input_name";
    }
    return $fh;
}

sub _shell_args_property_meta {
    # exclude this class' commands from shell arguments
    return grep {
            (
                $_->class_name eq __PACKAGE__
                    ? ($_->property_name eq 'model_id' ? 1 : 0)
                    : 1
            )
        } shift->SUPER::_shell_args_property_meta(@_);
}

# This is called by the infrastructure to appropriately classify abstract events
# according to their event type because of the "sub_classification_method_name" setting
# in the class definiton...
# TODO: replace with cleaner calculated property.
sub _resolve_subclass_name {
    my $class = shift;

    if (ref($_[0]) and $_[0]->isa(__PACKAGE__)) {
        my $event_type = $_[0]->event_type;
        return $class->_resolve_subclass_name_for_event_type($event_type);
    }
    elsif (my $event_type = $class->define_boolexpr(@_)->value_for('event_type')) {
        return $class->_resolve_subclass_name_for_event_type($event_type);
    }
    else {
        # this uses the model
        return $class->_get_sub_command_class_name(@_);
    }
}

# This is called by some legacy code.
sub class_for_event_type {
    my $self = shift;
    return $self->_resolve_subclass_name_for_event_type($self->event_type);
}

# This is called by both of the above.
sub _resolve_subclass_name_for_event_type {
    my ($class,$event_type) = @_;
    my @command_parts = split(' ',$event_type);
    my $genome_model = shift @command_parts;
    if ($genome_model eq 'genome'){
        #TODO, this is to accomodate eddie's refactoring, will need to switch once command name is redeployed
        $genome_model.= '-'.shift @command_parts;
    }
    if ($genome_model !~ m/genome-model/) {
        $class->error_message("Malformed event-type $event_type.  Expected it to begin with 'genome-model'");
        return;
    }

    foreach ( @command_parts ) {
        my @sub_parts = map { ucfirst } split('-');
        $_ = join('',@sub_parts);
    }

    my $class_name;
    my $loaded_class = 0;

    foreach my $type ('Event','Command') {
        $class_name = join('::', 'Genome::Model::' . $type, @command_parts);

        # TEMP
        $class_name =~ s/AddReads::/Build::ReferenceAlignment::/;

        if (eval {$class_name->class()}) {
            if ($class_name->isa('Genome::Model::Event')) {
                $loaded_class = 1;
                last;
            }
        }
    }

    unless ($loaded_class) {
        $class_name = 'Genome::Model::Event::Generic';
    }

#    print "resolved $class_name\n";

    return $class_name;
}

sub __errors__ {
    my ($self) = shift;

    my @tags = $self->SUPER::__errors__(@_);
    if ($self->model_id) {
        unless (Genome::Model->get($self->model_id)) {
            push @tags, UR::Object::Tag->create(
                                            type => 'invalid',
                                            properties => ['model_id'],
                                            desc => "There is no model with id ". $self->model_id,
                                        );
        }
    }
    # Once expunged instrument data events are abandoned, this will cause the DB commit to fail
    #if ($self->instrument_data_id && !Genome::InstrumentData->get($self->instrument_data_id)) {
    #    push @tags, UR::Object::Tag->create(
    #                                        type => 'invalid',
    #                                        properties => ['instrument_data_id'],
    #                                        desc => "There is no instrument data with instrument_data_id ". $self->instrument_data_id,
    #                                    );
    #}
    return @tags;
}

sub desc {
    my $self = shift;
    my $desc = $self->id .' (' . $self->event_type .')';
    if (defined $self->instrument_data_id) {
        $desc .= " for read set " . $self->instrument_data->full_name;
    }
    if (defined $self->ref_seq_id) {
        $desc .= " for refseq " . $self->ref_seq_id . " on build " . $self->build_id;
    }
    return $desc;
}

#< Bsub >#
sub bsub_rusage {
    return "-R 'select[model!=Opteron250 && type==LINUX64] span[hosts=1]'";
}

sub execute_with_bsub {
    my ($self, %params) = @_;
    my $queue = $params{bsub_queue};
    my $bsub_args = $params{bsub_args};
    my $dependency_hash_ref = $params{dependency_hash_ref};
    my $model_id = $self->model_id;

    ## should check if $self isa Command??
    $queue ||= $ENV{GENOME_LSF_QUEUE_BUILD_WORKER_ALT};

    $DB::single = $DB::stopper;

    my $dependency_expression = '';
    for my $dep_type (keys %{$dependency_hash_ref}) {
        for my $dep_value (@{$$dependency_hash_ref{$dep_type}}) {
            if ($dependency_expression eq '') {
                $dependency_expression .= "$dep_type($dep_value)";
            } else {
                $dependency_expression .= " && $dep_type($dep_value)";
            }
        }
    }
    if ($self->lsf_job_name) {
        $bsub_args .= ' -J "'. $self->lsf_job_name .'" ';
    }
    if (my $bsub_rusage = $self->bsub_rusage) {
        $bsub_args .= ' ' . $bsub_rusage;
    }

    my $class = $self->class;
    my $id = $self->id;

    my $cmd = "genome-model bsub-helper";

    my $event_id = $self->genome_model_event_id;
    my $log_dir = $self->resolve_log_directory;
    unless (-d $log_dir) {
        $self->create_directory($log_dir);
    }
    my $err_log_file = $self->error_log_file;
    my $out_log_file = $self->output_log_file;
    $bsub_args .= ' -o ' . $out_log_file . ' -e ' . $err_log_file;

    my $cmdline;
    { no warnings 'uninitialized';
        $cmdline = "bsub -q $queue -H $bsub_args" .
                   ($dependency_expression && " -w '$dependency_expression'") .
                   " $cmd --model-id $model_id --event-id $event_id ";
    }
    $self->status_message("Running command: " . $cmdline);

    # Header for output and error files
    for my $log_file ( $err_log_file, $out_log_file )
    {
        $DB::single = $DB::stopper;
        if(-e $log_file && (stat($log_file))[2] != 0100664) {
            unless ( chmod(0664, $log_file) )
            {
                $self->error_message("Can't chmod log file ($log_file)");
                return;
            }
        }
        my $fh = IO::File->new(">> $log_file");
        $fh->print("\n\n########################################################\n");
        $fh->print( sprintf('Submitted at %s: %s', UR::Context->current->now, $cmdline) );
        $fh->close;
    }
    my $bsub_output = `$cmdline`;
    my $retval = $? >> 8;

    if ($retval) {
        $self->error_message("bsub returned a non-zero exit code ($retval), bailing out");
            return;
    }
    my $bsub_job_id;
    if ($bsub_output =~ m/Job <(\d+)>/) {
        $bsub_job_id = $1;
    } else {
        $self->error_message('Unable to parse bsub output, bailing out');
        $self->error_message("The output was: $bsub_output");
        return;
    }
    return $bsub_job_id;
}

#< Scheduling, Abandon >#
sub schedule {
    my $self = shift;

    # FYI: user_name is set in create
    $self->event_status("Scheduled");
    $self->date_scheduled( UR::Context->current->now );
    $self->date_completed(undef);
    $self->retry_count(0) unless defined $self->retry_count;

    return 1;
}

sub abandon {
    my $self = shift;

    if ($self->event_status eq 'Abandoned') {
        return 1;
    }
    for my $next_event ($self->next_events) {
        $next_event->abandon;
    }
    if ($self->event_status =~ /Scheduled|Running/) {
        my $user = getpwuid($<);
        unless ($user eq 'apipe-builder' || $user eq $self->user_name) {
            my $build_id = $self->build_id || '';
            my $model_id = $self->model_id || '';
            $self->error_message('Attempted to abandon '. $self->event_status .' event '. $self->id ."(model_id:$model_id, build_id:$build_id) owned by ". $self->user_name);
            Carp::confess($self->error_message);
        }
        my $lsf_job_id = $self->lsf_job_id;
        if ($lsf_job_id) {
            my $cmd = 'bkill '. $lsf_job_id . ' > /dev/null 2>&1';
            my $rv = system($cmd);
            unless ($rv == 0) {
                $self->error_message('Failed execution of command '. $cmd);
            }
        }
    }
    $self->event_status("Abandoned");
    $self->date_completed(UR::Context->current->now);
    return 1;
}

#<>#
sub next_events {
    my $self = shift;
    my @next_events = Genome::Model::Event->get(
                                                prior_event_id => $self->id,
                                            );
    return @next_events;
}


sub is_reschedulable {
    my($self) = @_;

    return 1; # was part of bsub helper, may change implementation again
}

sub max_retries {
    2;
}

sub verify_prior_event {
    my $self = shift;

    if (my $prior_event = $self->prior_event) {
        unless ($prior_event->event_status eq 'Succeeded') {
            $self->error_message('Prior event '. $self->prior_event_id .' is not Succeeded');
            return;
        }
    }

    return 1;
}


#this method is just a wrapper that tries a database call, then tries to calculate the metric and store it if its not already in the db
sub get_metric_value {
    my $self = shift;
    my $metric_name = shift;

    return unless($metric_name);

    my $metric=$self->get_metric($metric_name);
    unless ($metric) {
         return "Not Found";
    }
    return $metric->value;
}

#this method is like gimme that metric from the database or fail if it doesn't exist.
#Its public for those that want to generate a view without impyling an hour computation for unknown values
sub get_metric {
    my $self = shift;
    my @metric_names = @_;

    unless (@metric_names) {
        @metric_names = $self->metrics_for_class;
    }

    if (@metric_names == 1) {
        # A hack to make a cgi script faster.  It preloads all the metrics, but the UR
        # cache system isn't able to find them because of the in-clause.
        return Genome::Model::Event::Metric->get(name => $metric_names[0], event_id => $self->id);
    } else {
        return Genome::Model::Event::Metric->get(
                                             name => \@metric_names,
                                             event_id => $self->id,
                                         );
    }
}

sub has_all_metrics {
    my $self = shift;

    my @metric_names = $self->metrics_for_class;
    for my $metric_name (@metric_names) {
        unless ($self->get_metric($metric_name)) {
            $self->error_message("Metric $metric_name does not exist for event_id ". $self->id);
            return 0;
        }
    }
    return 1;
}


#this method is like can i have that metric? no? then i'll make one!
sub resolve_metric {
    my $self = shift;
    my @metric_names = @_;

    unless (@metric_names) {
        @metric_names = $self->metrics_for_class;
    }

    my @metrics;
    for my $metric_name (@metric_names) {
        my $metric = $self->get_metric($metric_name);

        unless ($metric) {
            $metric = $self->generate_metric($metric_name);
            unless ($metric) {
                $self->error_message("Unable to generate requested metric $metric_name for event_id ". $self->id);
                next;
            }
        }
        push @metrics, $metric;
    }
    return @metrics;
}


#this method is called by resolve metric and it dynamically figures out the calculate method to call to store a new metric
sub generate_metric {
    my $self = shift;
    my @metric_names = @_;

    unless (@metric_names) {
        @metric_names = $self->metrics_for_class;
    }

    my @metrics;
    for my $metric_name (@metric_names) {
        my $metric = $self->get_metric($metric_name);

        my $calculate_method = '_calculate_'. $metric_name;
        unless ($self->can($calculate_method)) {
            $self->error_message("Event ". $self->id ." can not $calculate_method");
            next;
        }

        my $value = $self->$calculate_method;
        unless(defined $value) {
            $self->error_message("Non Fatal Metric Error(this doesn't kill the step):Value not defined for metric $metric_name using method $calculate_method");
            next;
        }
        if ($metric) {
            $metric->value($value);
        } else {
            $metric = $self->add_metric(
                                        name    => $metric_name,
                                        value   => $value,
                                    );
        }
        unless ($metric) {
            $self->error_message("Could not create/update metric $metric_name with value $value");
            return;
        }
        push @metrics, $metric;
    }
    if($self->can('cleanup_transient_properties')) {
        $self->cleanup_transient_properties();
    }
    return @metrics;
}


sub metrics_for_class {
    my $self = shift;
    $self->error_message("Please implement me! I do not have metrics_for_class");
    return 0;
}


sub lsf_job_name {
    my $self = shift;
    my $build = $self->build;
    unless ($build) {
        $self->error_message('No build found for event('. $self->id .')');
        die;
    }
    my $build_event = $build->build_event;
    unless ($build_event) {
        $self->error_message('No build event found for build id '. $build->id);
        die
    }
    my $stage_name = $build_event->resolve_stage_name_for_class($self->class);
    unless ($stage_name) {
        $self->error_message('Failed to resolve stage name for event('.
                             $self->id .','. $self->class .') with build event ('.
                             $build_event->id .','. $build_event->class .')');
        die;
    }
    return $self->model_id .'_'. $self->build_id .'_'. $stage_name .'_'. $self->id;
}

sub lsf_dependency_condition {
    my $self = shift;
    unless ($self->lsf_job_id) {
        return;
    }
    my ($job_info,$events) = $self->lsf_state($self->lsf_job_id);
    for my $entry (@$events) {
        my ($time, $attributes) = (@$entry);
        if ($$attributes{'Dependency Condition'}) {
            return $$attributes{'Dependency Condition'};
        }
    };
    return;
}

sub lsf_job_state {
    my $self = shift;
    unless ($self->lsf_job_id) {
        return;
    }
    my ($job_info,$events) = $self->lsf_state($self->lsf_job_id);
    if ($job_info) {
        return $$job_info{Status};
    }
    return;
}

sub lsf_pending_reasons {
    my $self = shift;
    unless ($self->lsf_job_id) {
        return;
    }
    unless ($self->lsf_job_state eq 'PEND') {
        return;
    }
    my @reasons;
    my ($job_info,$events) = $self->lsf_state($self->lsf_job_id);
    for my $entry (@$events) {
        my ($time, $attributes) = (@$entry);
        if ($$attributes{'PENDING REASON'}) {
            push @reasons, $$attributes{'PENDING REASON'};
        }
    }
    return @reasons;
}

sub lsf_state {
    my ($self, $lsf_job_id) = @_;

    my $spool = `bjobs -l $lsf_job_id 2>&1`;
    return if not defined $spool;
    return if ($spool =~ /No command 'bjobs' found/);
    return if ($spool =~ /Job <$lsf_job_id> is not found/);

    # this regex nukes the indentation and line feed
    $spool =~ s/\s{22}//gm;

    my @eventlines = split(/\n/, $spool);
    shift @eventlines;  # first line is white space
    my %jobinfo = ();

    my $jobinfoline = shift @eventlines;
    if (defined $jobinfoline) {
        # sometimes the prior regex nukes the white space between Key <Value>
        $jobinfoline =~ s/(?<!\s{1})</ </g;
        # parse out a line such as
        # Key <Value>, Key <Value>, Key <Value>
        while ($jobinfoline =~ /(?:^|(?<=,\s{1}))(.+?)(?:\s+<(.*?)>)?(?=(?:$|;|,))/g) {
            $jobinfo{$1} = $2;
        }
    }
    my @pending_reasons = ();
    foreach my $line (@eventlines) {
        if ($line =~ /PENDING REASONS:/) {
            @pending_reasons = (
                                UR::Context->current->now,
                                {}
                            );
            next;
        }
        if (scalar(@pending_reasons)) {
            if ($line =~ /^\s*$/) {
                last;
            }
            $line =~ s/^\s+//;
            $pending_reasons[1]->{'PENDING REASON'} = $line;
        }
    }

    my @events = ();
    foreach my $el (@eventlines) {
        $el =~ s/(?<!\s{1})</ </g;

        my $time = substr($el,0,21,'');
        substr($time,-2,2,'');

        # see if we really got the time string
        if ($time !~ /\w{3} \w{3}\s+\d{1,2}\s+\d{1,2}:\d{2}:\d{2}/) {
            # there's stuff we dont care about at the bottom, just skip it
            next;
        }

        my @entry = (
            $time,
            {}
        );

        while ($el =~ /(?:^|(?<=,\s{1}))(.+?)(?:\s+<(.*?)>)?(?=(?:$|;|,))/g) {
            $entry[1]->{$1} = $2;
        }
        push @events, \@entry;
    }
    if (scalar(@pending_reasons)) {
        push @events, \@pending_reasons;
    }
    return (\%jobinfo, \@events);
}

1;
