package Genome::Model::Tools::Picard::Base;

use strict;
use warnings;

use Genome;
use Carp qw(confess);
use File::Spec qw();

class Genome::Model::Tools::Picard::Base {
    is => 'Genome::Model::Tools::Picard',
    has_param => [
        lsf_queue => {
            default_value => Genome::Config::get('lsf_queue_build_worker_alt'),
            is_optional => 1,
            doc => 'queue to use when running in a workflow',
        },
    ],
};

# basename of the jar file (e.g., 'CleanSam.jar')
sub _jar_name {
    confess "sub _jar_name must be overridden in child class";
}

# full path to jar file, should not typically be overridden
sub _jar_path {
    my $self = shift;
    return $self->picard_path if $self->version_newer_than('1.123');
    return File::Spec->catfile($self->picard_path, $self->_jar_name);
}

# Java class as an array (e.g., qw(picard sam SortSam))
sub _java_class {
    confess "sub _java_class must be overridden in child class";
}

# Handles assembly of correct class name
sub _java_class_name {
    my $self = shift;
    my @class = $self->_java_class;
    if ($self->version_older_than('1.114')) {
        unshift @class, qw(net sf);
    }
    # Only one JAR file exists for version older than 1.124
    # Instead of the class name, only the sub-command name is necessary
    if ($self->version_newer_than('1.123')) {
        my @last_class = split(/\./,$class[-1]);
        #The last class in the array or . delimited string is the subcommand name
        return $last_class[-1];
    }
    return join '.', @class;
}

# picard params should be 'inputs' and also have picard_param_name defined
sub _picard_param_metas {
    my $class = shift;
    my @property_metas = $class->__meta__->properties(is_input => 1);
    return grep {$_->can("picard_param_name") && $_->picard_param_name} @property_metas;
}

# Handle conversion from boolean => 'true'/'false' and general picard arg
# formatting (argname=argvalue).
sub _format_picard_arg {
    my ($type, $name, $value) = @_;
    $value = $value ? 'true' : 'false' if $type eq 'Boolean';
    return sprintf('%s=%s', $name, Genome::Sys->quote_for_shell($value));
}

# given a property meta object (defining property name, and ostensibly
# picard_param_name), return the list of cmdline args that should be passed
# to picard based on the values set on the current object.
#
# this is where is_many properties get handled.
sub _picard_args_from_meta {
    my ($self, $meta) = @_;

    my $ur_name = $meta->property_name;
    # using // instead of ? : below made a lims perl test fail
    my $type = defined $meta->data_type ? $meta->data_type : 'String';
    my $picard_name = $meta->picard_param_name;

    return unless defined $self->$ur_name;

    # This works for things that are is_many or not.
    return map {_format_picard_arg($type, $picard_name, $_)} $self->$ur_name;
}

# command line args that come after the jar and class
#
# The default implementation just translates all the inputs that have
# picard_param_name defined. Boolean values get translated to true/false,
# and is many attributes get repeated, e.g.,
#
#  create_md5_file => {..., is => 'Boolean'}
#       when set to something that evalutes to false
#       yields
#  qw(CREATE_MD5_FILE=false)
#
# and
#
#  program_name => {..., is_many => 1}
#       when set to ['a', 'b', 'c']
#       yields
#  qw(PROGRAM_NAME=a PROGRAM_NAME=b PROGRAM_NAME=c)
#
# Override for custom behavior.
sub _cmdline_args {
    my $self = shift;

    # CREATE_MD5_FILE will cause errors for certain picard commands (e.g., CompareSAMs)
    # the picard default value is false.
    my @metas = grep {$_->property_name ne 'create_md5_file'} $self->_picard_param_metas;
    my @args = map {$self->_picard_args_from_meta($_)} @metas;

    push @args, 'CREATE_MD5_FILE=true' if $self->create_md5_file;
    return @args;
}

# this is your chance to make a loud sound and die if a caller asks for
# something ridiculous. params can also be edited here if that is somehow
# appropriate.
sub _validate_params {
    # example:
    # my $self = shift;
    # die unless -s $self->input_file;
    # $self->be_noisy(0) unless $self->log_file;
}

# want to pass extra stuff to Genome::Sys->shellcmd? return a hash (not ref)
# of those params here.
sub _shellcmd_extra_params {
    # example:
    # my $self = shift
    # return (
    #   input_files => [$self->input_file, $self->reference_sequence],
    #   output_files => [$self->output_file]
    #   );
}

# string to be appended after the full command line (e.g., '> foo.out')
#
# there was a bug that led to some commands being built like:
#
#   java ... foo.jar class ARGS > output_file MORE_ARGS
#
# that's bad. this sort of conflicts with the idea of log_file, but commands
# like "BamIndexStats" that have no native OUTPUT option (i.e., can only spew
# to stdout) support the existence of this function.
sub _redirects {
    my $self = shift;
    return sprintf(">> %s", $self->log_file) if $self->log_file;
}

sub build_cmdline_list {
    my $self = shift;
    my $jvm_options = $self->additional_jvm_options || '';
    my @cmdline = (
        $self->java_interpreter,
        sprintf('-Xmx%dm', int(1024 * $self->maximum_memory)),
        sprintf('-XX:MaxPermSize=%dm', $self->maximum_permgen_memory),
        split(' ', $jvm_options),
    );
    if ($self->version_newer_than('1.123')) {
        push @cmdline, (
            '-jar',
            $self->_jar_path,
        );
    } else {
        push @cmdline, (
            '-cp',
            sprintf('/usr/share/java/ant.jar:%s', $self->_jar_path),
        );
    }
    unless ($self->_java_class_name) {die;}
    push @cmdline, (
        $self->_java_class_name,
        $self->_cmdline_args,
    );
    return @cmdline;
}

sub build_cmdline_string {
    my ($self, $cmd) = @_;

    my $java_vm_cmd = join(' ', $self->build_cmdline_list);

    my $redirects = $self->_redirects;
    $java_vm_cmd = join(" ", $java_vm_cmd, $redirects) if $redirects;

    return $java_vm_cmd;
}

# override with care.
sub execute {
    my $self = shift;

    $self->enforce_minimum_version_required;
    $self->_validate_params;

    my %params = (
        cmd => $self->build_cmdline_string,
        $self->_shellcmd_extra_params
        );

    return Genome::Sys->shellcmd(%params);
}

1;
