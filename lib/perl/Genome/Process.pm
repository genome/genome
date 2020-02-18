package Genome::Process;

use strict;
use warnings FATAL => 'all';
use Genome;
use File::Spec;
use Params::Validate qw(validate_pos :types);
use Data::Dump qw(pp);
use Scalar::Util qw();
use Try::Tiny qw(try catch);
use JSON qw(to_json);
use List::MoreUtils qw(uniq any);
use Genome::Disk::Group::Validate::GenomeDiskGroups;
use Genome::Utility::Inputs qw(encode);
use Cwd qw(abs_path);
use File::DirCompare;

class Genome::Process {
    roles => [qw(
        Genome::Role::ObjectWithCreatedBy
        Genome::Role::Notable
        Genome::Role::SoftwareResultRequestor
        Genome::Role::Comparable
    )],
    is_abstract => 1,
    table_name => 'process.process',
    schema_name => 'GMSchema',
    data_source => 'Genome::DataSource::GMSchema',

    subclassify_by => 'subclass_name',
    id_generator => '-uuid',
    id_by => [
        id => {is => 'Text'},
    ],

    subclass_description_preprocessor => 'Genome::Process::_preprocess_for_inputs',
    attributes_have => [
        is_input => {
            is => 'Boolean',
            default => 0,
        }
    ],

    has => [
        disk_allocation => {
            is => 'Genome::Disk::Allocation',
            id_by => 'disk_allocation_id',
            doc => 'A disk allocation to store logs and other process meta-data',
        },
        status => {
            is => 'Text',
            doc => 'The current status of the process',
        },
        status_events => {
            is => 'Genome::Process::StatusEvent',
            where => [ -order_by => [ "timestamp" ] ],
            is_many => 1,
            reverse_as => 'process',
        },
        created_at => {
            is => 'Timestamp',
            is_optional => 1,
            doc => "The time that status was set to 'New'",
        },
        started_at => {
            is => 'Timestamp',
            is_optional => 1,
            doc => "The time that status was set to 'Running'",
        },
        ended_at => {
            is => 'Timestamp',
            is_optional => 1,
            doc => "The time that status was set to 'Crashed' or 'Succeeded'",
        },
        software_revision => {
            is => 'Text',
            is_mutable => 0,
            doc => 'The version of the Genome code the process was created with',
        },
        results => {
            is => 'Genome::SoftwareResult',
            via => 'result_users',
            to => 'software_result',
        },
        result_users => {
            is => 'Genome::SoftwareResult::User',
            is_many => 1,
            reverse_as => 'user',
        },
        subclass_name => {
            is => 'Text',
            is_mutable => 0,
            doc => 'Used by UR to instantiate the correct sub-class',
        },
        inputs => {
            is => 'Genome::Process::Input',
            is_many => 1,
            is_optional => 1,
            reverse_as => 'process',
        }

    ],
    doc => 'A base class to manage meta-data related to running a process (workflow)',
};

sub symlink_results {
    my $self = shift;
    $self->status_message("Process (%s) is a %s which does not symlink results.",
        $self->id, $self->class);
    return 1;
}

sub run {
    my $self = shift;
    my %p = Params::Validate::validate(@_, {
            workflow_xml => {type => SCALAR},
            workflow_inputs => {type => HASHREF},
    });

    my $transaction = UR::Context::Transaction->begin();
    $self->create_disk_allocation();

    $self->_write_workflow_file($p{workflow_xml});
    $self->_write_inputs_file($p{workflow_inputs});

    local $ENV{UR_DUMP_DEBUG_MESSAGES} = 1 unless
        exists $ENV{UR_DUMP_DEBUG_MESSAGES};
    local $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES} = 1 unless
        exists $ENV{UR_COMMAND_DUMP_DEBUG_MESSAGES};
    local $ENV{UR_DUMP_STATUS_MESSAGES} = 1 unless
        exists $ENV{UR_DUMP_STATUS_MESSAGES};
    local $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES} = 1 unless
        exists $ENV{UR_COMMAND_DUMP_STATUS_MESSAGES};

    if ($ENV{UR_DBI_NO_COMMIT}) {
        return $self->_execute_process($transaction);
    } else {
        my $backend = Genome::Config::get('workflow_builder_backend');
        if ($backend eq 'ptero') {
            $self->fatal_message('Sorry, PTero is no longer supported. Configure a different workflow_builder_backend.');
        } elsif ($backend eq 'inline') {
            $self->_execute_process($transaction);
        } elsif ($backend eq 'simple') {
            $self->_dispatch_process($transaction);
        } else {
            $self->fatal_message('Unknown backend chosen to run process: %s', $backend);
        }
    }
}

sub run_and_wait {
    my $self = shift;
    my %p = Params::Validate::validate(@_, {
            workflow_xml => {type => SCALAR},
            workflow_inputs => {type => HASHREF},
    });

    my $transaction = UR::Context::Transaction->begin();
    $self->create_disk_allocation();
    $self->_write_workflow_file($p{workflow_xml});
    $self->_write_inputs_file($p{workflow_inputs});
    return $self->_execute_process($transaction);
}

sub _write_workflow_file {
    my $self = shift;
    my $workflow_xml = shift;

    Genome::Sys->write_file($self->workflow_file, $workflow_xml);
}

sub _write_inputs_file {
    my $self = shift;
    my $inputs = shift;

    Genome::Sys->write_file($self->inputs_file,
        to_json(encode($inputs), {pretty => 1, canonical => 1}));
}

sub write_environment_file {
    my $self = shift;

    Genome::Sys->write_file($self->environment_file,
        to_json({%ENV}, {pretty => 1, canonical => 1}));
}

sub _execute_process {
    my $self = shift;
    my $transaction = shift;

    my $cmd = Genome::Process::Command::Run->create(
        process => $self,
        update_with_commit => 0
    );
    if (my $rv = $cmd->execute()) {
        if ($transaction->commit()) {
            $self->status_message("Successfully ran process (%s)",
                $self->id);
            return $rv;
        } else {
            $transaction->rollback();
            die sprintf("Failed to commit process (%s): %s",
                $self->id, $transaction->error_message || 'Reason Unknown');
        }
    } else {
        die sprintf("Failed to run process (%s) for some reason",
            $self->id);
    }
}

sub _submit_process {
    my $self = shift;
    my $transaction = shift;

    unless ($transaction->commit()) {
        $transaction->rollback();
        die sprintf("Failed to submit process (%s): %s",
            $self->id, $transaction->error_message || 'Reason Unknown');
    }

    my $cmd = Genome::Process::Command::Run->create(
        process => $self, update_with_commit => 0);
    $cmd->submit();
}

sub _dispatch_process {
    my $self = shift;
    my $transaction = shift;

    unless ($transaction->commit()) {
        $transaction->rollback();
        $self->fatal_message(
            "Failed to dispatch process (%s): %s",
            $self->id, $transaction->error_message || 'Reason Unknown'
        );
    }

    my $job_group = join('/',
        Genome::Config::get('lsf_job_group'),
        Genome::Sys->username,
        'process',
    );

    my $log_file_base = File::Spec->join($self->log_directory, 'process-%J');

    my $job_id = Genome::Sys->bsub(
        queue => Genome::Config::get('lsf_queue_build_worker'),
        cmd => [qw(genome process run), $self->id],
        hold_job => 1,
        err_file => "${log_file_base}.err",
        log_file => "${log_file_base}.out",
        project => $self->lsf_project_name,
        job_group => $job_group,
        resource_string => Genome::Config::get('lsf_resource_cwl_runner'),
    );

    unless ($job_id) {
        $self->fatal_message('Failed to dispatch process (%s): bsub did not succeed.', $self->id);
    }
    $self->lsf_job_id($job_id);
    Genome::Sys::CommitAction->create(
        on_commit => sub {
            Genome::Sys->shellcmd(cmd => ['bresume', $job_id]);
        },
    );

    return 1;
}

sub create {
    my $class = shift;

    if ($class->__meta__->is_abstract) {
        return $class->SUPER::create(@_); #UR will re-call with the concrete class.
    }

    # we do not support creating from a boolexpr because the boolexpr logic
    # will silently filter out some invalid param-values.
    my %params = @_;

    unless (exists $params{software_revision}) {
        $params{software_revision} = Genome::Sys->snapshot_revision();
    }

    my $self = $class->SUPER::create(
        $class->_preprocess_params_for_create(%params));
    return unless $self;

    $self->update_status('New');
    $self->bail_out_if_input_errors(\%params);

    return $self;
}

# because inputs are immutable, there are no opportunities to correct these
# errors, so we should die.
sub bail_out_if_input_errors {
    my $self = shift;
    my $params = shift;

    my @errors = $self->__input_errors__;
    if (@errors) {
        my $class = $self->class;
        $self->print_errors(@errors);
        $self->delete();
        die sprintf("Failed to create (%s) with params: %s",
            $class, pp($params));
    }
}

my $SET_TIMESTAMP_ON_STATUS = {
    New => 'created_at',
    Running => 'started_at',
    Crashed => 'ended_at',
    Succeeded => 'ended_at',
};

sub update_status {
    my ($self, $new_status) = validate_pos(@_,
        {type => OBJECT}, {type => SCALAR});

    my $old_status = $self->status;

    my $now = UR::Context->current->now;
    my $event = Genome::Process::StatusEvent->create(
        process => $self,
        old_status => $old_status,
        new_status => $new_status,
        timestamp => $now,
    );

    if ($event) {
        my $timestamp_accessor = $SET_TIMESTAMP_ON_STATUS->{$new_status};
        if ($timestamp_accessor) {
            $self->$timestamp_accessor($now);
        }
        $self->status($new_status);
    } else {
        die sprintf("Cannot transition Process (%s) from (%s) to (%s)",
            $self->id, $old_status, $new_status);
    }
    return $self->status;
}

sub lsf_job_id_header {
    return 'workflow lsf job_id';
}

sub lsf_job_id {
    my $self = shift;

    my $new_value = shift;
    if ($new_value) {
        my $existing = $self->notes(header_text => $self->lsf_job_id_header);
        if ($existing) {
            $existing->delete;
        }
        my $n = $self->add_note(header_text => $self->lsf_job_id_header, body_text => $new_value);
        $n->body_text($new_value); #ignore any system-generated sudo message for this note.
    } else {
        my $id_note = $self->notes(header_text => $self->lsf_job_id_header);
        if ($id_note) {
            return $id_note->body_text;
        } else {
            return;
        }
    }
}

sub workflow_name {
    my $self = shift;
    return sprintf('Genome::Process(%s)', $self->id);
}

sub lsf_project_name {
    my $self = shift;
    return $self->workflow_name;
}

sub environment_file {
    my $self = shift;
    return File::Spec->join($self->metadata_directory, 'environment.json');
}

sub inputs_file {
    my $self = shift;
    return File::Spec->join($self->metadata_directory, 'inputs.json');
}

sub workflow_file {
    my $self = shift;
    return File::Spec->join($self->metadata_directory, 'workflow.xml');
}

sub log_directory {
    my $self = shift;
    return File::Spec->join($self->metadata_directory, 'logs');
}

sub metadata_directory {
    my $self = shift;
    return $self->disk_allocation->absolute_path;
}

sub create_disk_allocation {
    my ($self, $kb_requested) = validate_pos(@_, {type => OBJECT},
        {default => 50 * 1024},
    );

    if ($self->disk_allocation) {
        $self->status_message("Process (%s) already has a disk_allocation (%s)",
            $self->id, $self->disk_allocation->id);
        return $self->disk_allocation;
    }

    my $allocation_path = File::Spec->join('model_data', 'process', $self->id);

    my %args = (
        disk_group_name => Genome::Config::get('disk_group_models'),
        allocation_path => $allocation_path,
        kilobytes_requested => $kb_requested,
        owner_class_name => $self->class,
        owner_id => $self->id,
    );
    $self->debug_message("Attempting to create a disk_allocation for process " .
        "(%s) with arguments: %s", $self->id, pp(\%args));
    my $disk_allocation = Genome::Disk::Allocation->create(%args);

    if ($disk_allocation) {
        $self->disk_allocation($disk_allocation);
        $self->status_message("Process (%s) now has disk_allocation (%s)",
            $self->id, $self->disk_allocation->id);
        $self->_ensure_disk_allocation_gets_cleaned_up();

        Genome::Sys->create_directory($self->log_directory);
        return $self->disk_allocation;
    } else {
        die sprintf("Failed to create disk allocation with " .
            "arguments: %s", pp(\%args));
    }
}

#XXX not sure why but this only works with a UR software transaction
sub _ensure_disk_allocation_gets_cleaned_up {
    my $self = shift;

    my $cleanup_closure = $self->_disk_allocation_cleanup_closure();
    my $create_disk_allocation = UR::Context::Transaction->log_change(
            $self, 'UR::Value', $self->disk_allocation_id,
            'external_change', $cleanup_closure);
    unless ($create_disk_allocation) {
        die sprintf("Couldn't log allocation (%s) created",
            $self->disk_allocation_id);
    }
}

sub _disk_allocation_cleanup_closure {
    my $self = shift;
    my $observer = shift;

    my $allocation_id = $self->disk_allocation_id;
    my $process_id = $self->id;
    my $remove_allocation = sub {
        print "Now deleting disk allocation ($allocation_id) associated " .
            "with process ($process_id)\n";
        ${$observer}->delete if $observer;
        my $allocation = Genome::Disk::Allocation->get($allocation_id);
        if ($allocation) {
            $allocation->delete;
        }
    };
    $self->debug_message("Created closure to delete disk allocation (%s) " .
        "associated with process (%s)", $allocation_id, $process_id);
    return $remove_allocation;
}


sub _preprocess_for_inputs {
    my ($class, $class_description) = @_;

    while (my ($property_name, $property_description) =
            each %{$class_description->{has}}) {
        if ($property_description->{is_input}) {
            _expand_is_input($class, $property_description, $property_name);
        }
    }
    return $class_description;
}

sub _expand_is_input {
    my $class = shift;
    my $property_description = shift;
    my $property_name = shift;

    $property_description->{'is_mutable'} = 0;
    $property_description->{'via'} = 'inputs';
    $property_description->{'is_delegated'} = 1;
    $property_description->{'to'} = 'value_id';
    $property_description->{'where'} = [
        name => $property_name,
        '-order_by' => 'array_index',
    ];
    if (my $data_type = $property_description->{data_type}) {
        my $value_class_name = UR::Object::Property->_convert_data_type_for_source_class_to_final_class(
            $data_type, $class);
        unless ($value_class_name->isa("UR::Value")) {
            $property_description->{'to'} = 'value';
        }
    }
}

sub _preprocess_params_for_create {
    my $self = shift;
    my %params = @_;

    my $input_params_list = $params{inputs} || [];
    for my $property ($self->__meta__->properties) {
        my $property_name = $property->property_name;

        if ($property->{is_input}) {
            my $value = delete $params{$property_name};
            next unless defined($value);
            if ($property->{is_many}) {
                push @$input_params_list,
                    @{_get_input_params_list($property, $value)};
            } else {
                push @$input_params_list, _get_input_params($property, $value);
            }
        }
    }
    $params{inputs} = $input_params_list;
    return %params;
}

sub _get_input_params_list {
    my ($property, $value_list) = @_;

    if (!ref $value_list || (ref $value_list && ref $value_list ne 'ARRAY')) {
        die sprintf(
            "Cannot set 'is_many' input named (%s) with non-arrayref (%s)",
            $property->property_name, $value_list);
    }

    my $input_params_list = [];
    if (defined($value_list) && scalar(@$value_list)) {
        my $index = 0;
        for my $value (@$value_list) {
            push @$input_params_list,
                _get_input_params($property, $value, $index);
            $index++;
        }
    }
    return $input_params_list;
}

sub _get_input_params {
    my ($property, $value, $index) = validate_pos(@_, 1, 1, 0);

    my $input_params = {name => $property->property_name};
    $input_params->{array_index} = $index if defined($index);

    if ($property->to eq 'value') {
        if (Scalar::Util::blessed($value) &&
                $property->is_valid_storage_for_value($value)) {
            $input_params->{value} = $value;
        } else {
            die sprintf("Cannot set input named (%s) to value (%s), " .
                "because it is not a (%s)",
                $property->property_name, $value, $property->data_type);
        }
    } else {
        if (ref($value)) {
            die sprintf("Cannot set input named (%s) to value (%s), " .
                "because it is a reference (%s)",
                $property->property_name, $value, ref($value));
        } else {
            $input_params->{value_class_name} = 'UR::Value',
            $input_params->{value_id} = $value,
        }
    }
    return $input_params;
}

sub __errors__ {
    my $self = shift;

    return (
        $self->SUPER::__errors__(),
        $self->__input_errors__,
    );
}

sub __input_errors__ {
    my $self = shift;

    my @errors;
    for my $property ($self->__meta__->properties(
            is_input=>1, is_optional=>0)) {
        my $property_name = $property->property_name;

        my @values = $self->$property_name;
        unless (scalar(@values)) {
            push @errors, UR::Object::Tag->create(
                type => 'invalid',
                properties => [$property_name],
                desc => 'no value for required input'
            );
        }
    }
    return @errors;
}

sub print_errors {
    my ($self, @errors) = @_;

    for my $error (@errors) {
        my @properties = $error->properties;
        $self->error_message("Property " .
            join(',', map { "'$_'" } @properties) .
            ': ' . $error->desc);
    }
    return;
}


sub delete {
    my $self = shift;

    for my $status_event ($self->status_events) {
        $status_event->delete;
    }

    for my $input ($self->inputs) {
        $input->delete;
    }

    my $observer = Genome::Sys::CommitAction->create(
        on_commit => $self->_disk_allocation_cleanup_closure()
    );
    if ($observer) {
        $self->status_message("Registered observer to delete disk allocation " .
            "(%s) upon commit", $self->disk_allocation_id);
    } else {
        $self->error_message("Failed to register observer to delete disk " .
            "allocation (%s), you need to delete it manually, after commiting",
            $self->disk_allocation_id);
    }

    return $self->SUPER::delete(@_);
}

sub unique_results {
    my $self = shift;
    my @results = $self->results;
    return uniq @results;
}

sub result_with_label {
    my ($self, $label) = validate_pos(@_, 1, 1);

    my @results = grep {() = $_->users(label => $label)} $self->unique_results;
    if (scalar(@results) != 1) {
        die sprintf("Found (%d) results with label (%s), but expected only one, they are: %s",
            scalar(@results), $label, pp(map {$_->id} @results));
    } else {
        return $results[0];
    }
}

sub is_cle_verified {
    my $self = shift;

    for my $result ($self->unique_results) {
        unless ($self->result_is_index_result($result) or $self->result_is_on_cle_disk_group($result)) {
            return 0;
        }
    }
    return 1;
}

sub result_is_on_cle_disk_group {
    my $self = shift;
    my $result = shift;

    my $allocation = $result->disk_allocation;
    my $disk_group_name = $allocation->disk_group_name;

    return Genome::Disk::Group::Validate::GenomeDiskGroups::is_cle_disk_group_name($disk_group_name);
}

sub result_is_index_result {
    my $self = shift;
    my $result = shift;

    return any { $result->isa($_) } (qw(Genome::Db::Ensembl::Api Genome::Db::Ensembl::VepCache));
}

sub compare_output {
    my ($self, $other_process_id) = @_;
    my $process_id = $self->id;

    unless (defined ($other_process_id)) {
        die $self->error_message('Require process ID argument!');
    }
    my $other_process = Genome::Process->get($other_process_id);
    unless ($other_process) {
        die $self->error_message('Could not get process for ID (%s)', $other_process_id);
    }

    unless ($self->class eq $other_process->class) {
        die $self->error_message('Processes (%s) and (%s) are not of the same type', $process_id, $other_process_id);
    }

    my %diffs = $self->_compare_output_files($other_process);

    return %diffs;
}

sub _compare_output_files {
    my ($self, $other_process) = @_;

    my $output_dir = Genome::Sys->create_temp_file_path();
    my $other_output_dir = Genome::Sys->create_temp_file_path();

    $self->symlink_results($output_dir);
    $other_process->symlink_results($other_output_dir);

    return $self->_compare_output_directories($output_dir, $other_output_dir, $other_process);
}

sub _compare_output_directories {
    my ($self, $output_dir, $other_output_dir, $other_process) = @_;

    my %diffs;
    File::DirCompare->compare(
        $other_output_dir,
        $output_dir,
        sub {
            my ($other_file, $file) = @_;
            my $template = 'no %s %s found for process %s';
            my $other_target = defined($other_file)? abs_path($other_file) : undef;
            my $target = defined($file) ? abs_path($file) : undef;
            if (! $target) {
                my $type = ( -f $other_target ? 'file' : 'directory' );
                my $rel_path = File::Spec->abs2rel($other_file, $other_output_dir);
                $diffs{$rel_path} = sprintf($template, $type, $rel_path, $self->id);
            } elsif (! $other_target) {
                my $type = ( -f $target ? 'file' : 'directory' );
                my $rel_path = File::Spec->abs2rel($file, $output_dir);
                $diffs{$rel_path} = sprintf($template, $type, $rel_path, $other_process->id);
            } else {
                if (-l $file || -l $other_file) {
                    if (-d $target && -d $other_target) {
                        my %additional_diffs = $self->_compare_output_directories($target, $other_target, $other_process);
                        %diffs = (%diffs, %additional_diffs);
                    }
                    elsif (-f $target && -f $other_target) {
                        if ($self->file_comparer->compare($target, $other_target)) {
                            $diffs{File::Spec->abs2rel($file, $output_dir)} = sprintf(
                                'symlinks are not the same (diff -u %s %s) (diff -u %s %s)',
                                $file, $other_file, $target, $other_target
                            );
                        }
                        else {
                            #Symlink targets are in fact the same content - do nothing
                        }
                    }
                    else {
                        $diffs{File::Spec->abs2rel($file, $output_dir)} = sprintf(
                            'paths are not the same (diff -u %s %s) (diff -u %s %s)',
                            $file, $other_file, $target, $other_target
                        );
                    }
                }
                elsif ($self->file_comparer->compare($target, $other_target)) {
                    $diffs{File::Spec->abs2rel($file, $output_dir)} = sprintf(
                        'files are not the same (diff -u %s %s)',
                        $file, $other_file
                    );
                }
                else {
                    # Files are the same
                }
            }
        },
    );
    return %diffs;
}

1;
