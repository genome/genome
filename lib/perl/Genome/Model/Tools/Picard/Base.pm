package Genome::Model::Tools::Picard::Base;

use strict;
use warnings;

use Genome;
use Carp qw(confess);
use File::Spec qw();

class Genome::Model::Tools::Picard::Base {
    is => 'Genome::Model::Tools::Picard',
};

# basename of the jar file (e.g., 'CleanSam.jar')
sub _jar_name {
    confess "sub _jar_name must be overridden in child class";
}

# full path to jar file, should not typically be overridden
sub _jar_path {
    my $self = shift;
    return File::Spec->catfile($self->picard_path, $self->_jar_name);
}

# java class name within the jar to run (e.g., 'net.sf.picard.sam.CleanSam')
sub _java_class_name {
    confess "sub _java_class_name must be overridden in child class";
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
    return sprintf '%s=%s', $name, $value;
}

# given a property meta object for (defining property name, and ostensibly
# picard_param_name), return the list of cmdline args that should be passed
# to picard based on the values set on the current object.
#
# this is where is_many properties get handled.
sub _picard_args_from_meta {
    my ($self, $meta) = @_;

    my $ur_name = $meta->property_name;
    my $type = $meta->data_type // 'String';
    my $picard_name = $meta->picard_param_name;

    return unless defined $self->$ur_name;

    # This handles both is_many => true / false.
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
    return map {$self->_picard_args_from_meta($_)} $self->_picard_param_metas;
}

# this is your chance to make a loud sound and die if a caller asks for
# something ridiculous
sub _validate_params {
    # example:
    # my $self = shift;
    # $self->enforce_minimum_version('1.85');
    # die unless -s $self->input_file;
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

sub build_cmdline_string {
    my ($self, $cmd) = shift;
    my $jvm_options = $self->additional_jvm_options || '';

    my $java_vm_cmd = sprintf '%s -Xmx%dm -XX:MaxPermSize=%dm %s -cp /usr/share/java/ant.jar:%s %s %s',
        $self->java_interpreter,
        int(1024 * $self->maximum_memory),
        $self->maximum_permgen_memory,
        $jvm_options,
        $self->_jar_path,
        $self->_java_class_name,
        join(" ", $self->_cmdline_args);

    # hack to workaround premature release of 1.123 with altered classnames
    if ($java_vm_cmd =~ /picard-tools.?1\.123/) {
        $java_vm_cmd =~ s/net\.sf\.picard\./picard./;
    }

    $java_vm_cmd = join(" ", $java_vm_cmd, $self->_redirects);

    return $java_vm_cmd;
}

# override with care.
sub execute {
    my $self = shift;

    $self->_validate_params;

    my %params = (
        cmd => $self->build_cmdline_string,
        $self->_shellcmd_extra_params
        );

    return Genome::Sys->shellcmd(%params);
}

1;
