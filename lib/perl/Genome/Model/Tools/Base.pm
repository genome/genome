package Genome::Model::Tools::Base;

use strict;
use warnings;

use Genome;
use Carp qw(confess);

class Genome::Model::Tools::Base {
    is => 'Command::V2',
    is_abstract => 1,
    subclass_description_preprocessor => 'Genome::Model::Tools::Base::_expand_version_param',
    attributes_have => [
        tool_arg_name => {
            is => 'String',
            is_optional => 1,
        },
        tool_bare_arg_position => {
            is => 'Integer',
            is_optional => 1,
        },
        tool_input_file => {
            is => 'Boolean',
            is_optional => 1,
        },
    ],
    has_optional_output => [
        output_files => {
            is_many => 1,
        },
        input_files => {
            is_many => 1,
        },
    ],
};

sub versions {
    confess('please define versions in sub class!');
    # example:
    # my $class = shift;
    # my %VERSIONS = (
    #     'v1' => 'path',
    # );
    # return %VERSIONS;
}

sub _expand_version_param {
    my ($class,$desc) = @_;

    #only preprocess if we're an immediate child
    my $is = exists $desc->{is} ? $desc->{is} : [];
    unless (grep { $_ eq __PACKAGE__ } @$is) {
        return $desc;
    }

    my $subclass_name = $desc->{class_name};

    $desc->{has}->{version} = {
        is => 'UR::Value::Text',
        is_param => 1,
        doc => 'The version of the tool to use.',
        valid_values => [available_versions($subclass_name)],
    };

    return $desc;
}

sub available_versions {
    my $class = shift;
    my %versions = $class->versions;
    return sort keys %versions;
}

sub path_for_version {
    my ($class, $version) = @_;

    unless ($version) {
        die('No version defined!');
    }

    my %versions = $class->versions;
    my $path = $versions{$version};

    return $path if defined $path;

    $class->fatal_message('No path found for version: '. $version);
}

# Override in sub classes that have sub-command structure
# return a list of the base command and the sub command,
# ex. return ($self->path_for_version($self->version),$self->_tool_subcommand_name);
sub tool_path {
    my $self = shift;
    return $self->path_for_version($self->version);
}

sub _resolve_output_files {
    # override in sub class

    # example:
    # my @output_files = qw(file1.bam file2.bam);
    # return @output_files;
}

# Input Files
sub _resolve_input_files {
    my $self = shift;

    my @metas = $self->_tool_input_file_metas;

    my @input_files = map { $self->_tool_input_files_from_meta($_) } @metas;
    $self->input_files(\@input_files);
    return @input_files;
}

# tool params should be 'params' and also have tool_arg_name defined
sub _tool_param_metas {
    my $class = shift;
    my @param_property_metas = $class->__meta__->properties(is_param => 1);
    return grep {$_->can("tool_arg_name") && $_->tool_arg_name} @param_property_metas;
}

# tool inputs should be 'inputs' and also have tool_arg_name defined
sub _tool_input_metas {
    my $class = shift;
    my @input_property_metas = $class->__meta__->properties(is_input => 1);
    return grep {$_->can("tool_arg_name") && $_->tool_arg_name} @input_property_metas;
}

# tool inputs should be 'inputs' and also have tool_input_bare_arg_position defined
sub _tool_bare_arg_input_metas {
    my $class = shift;
    my @input_property_metas = $class->__meta__->properties(is_input => 1);
    my @tool_input_metas = grep {$_->can("tool_bare_arg_position") && $_->tool_bare_arg_position} @input_property_metas;
    return sort { $a->tool_bare_arg_position <=> $b->tool_bare_arg_position } @tool_input_metas;
}

# tool input files should be 'inputs' and also have tool_input_file defined
sub _tool_input_file_metas {
    my $class = shift;
    my @property_metas = $class->__meta__->properties(is_input => 1);
    return grep {$_->can("tool_input_file") && $_->tool_input_file} @property_metas;
}

# Handle conversion from boolean => 'true'/'false' and general tool arg
# formatting (argname=argvalue).
sub _format_tool_arg {
    my ($type, $name, $value) = @_;

    if ($type eq 'Boolean') {
        return '--'. $name;
    } else {
        return '--'.$name, $value;
    }
}

# given a property meta object (defining property name, and ostensibly
# tool_arg_name), return the list of cmdline args that should be passed
# to the tool based on the values set on the current object.
#
# this is where is_many properties get handled.
sub _tool_args_from_meta {
    my ($self, $meta) = @_;

    my $ur_name = $meta->property_name;
    my $type = $meta->data_type // 'String';
    my $tool_arg_name = $meta->tool_arg_name;

    return unless defined $self->$ur_name;
    if ($type eq 'Boolean' && $self->$ur_name eq '0') {
        return;
    }

    # This works for things that are is_many or not.
    return map {_format_tool_arg($type, $tool_arg_name, $_)} $self->$ur_name;
}

sub _tool_input_files_from_meta {
    my ($self, $meta) = @_;

    my $ur_name = $meta->property_name;

    return unless defined $self->$ur_name;

    # This works for things that are is_many or not.
    return map {$_} $self->$ur_name;
}

# command line args that come after the executable
#
# The default implementation just translates all the inputs that have
# tool_arg_name defined. Boolean values get translated to true,
# false are not translated,
# and is many attributes get repeated, e.g.,
#
#  exome => {..., is => 'Boolean'}
#       when set to something that evalutes to false the argument is not returned
#       when set to true it simply returns the correct tool argument as defined by tool_arg_name
#  ex. ('--exome')
#
# and
#
#  bams => {..., is_many => 1}
#       when set to ['a', 'b', 'c'] yields:
#  ex. ('--bam "a"', '--bam "b"','--bam "c"')
#

sub _tool_param_args {
    my $self = shift;

    my @metas = $self->_tool_param_metas;
    my @args = map {$self->_tool_args_from_meta($_)} @metas;

    return @args;
}

sub _tool_input_args {
    my $self = shift;

    my @metas = $self->_tool_input_metas;
    my @args = map {$self->_tool_args_from_meta($_)} @metas;

    my @bare_arg_metas = $self->_tool_bare_arg_input_metas;

    for my $bare_arg_meta (@bare_arg_metas) {
        my $ur_name = $bare_arg_meta->property_name;
        push @args, $self->$ur_name;
    }

    return @args;
}

# this is your chance to make a loud sound and die if a caller asks for
# something ridiculous. params can also be edited here if that is somehow
# appropriate.
sub _validate_params {
    # example:
    # my $self = shift;
    # $self->enforce_minimum_version('1.85');
    # die unless -s $self->input_file;
    # $self->be_noisy(0) unless $self->log_file;

    #TODO: enforce one of bams, normal_bam or tumor_bam inputs are required
}

# want to pass extra stuff to Genome::Sys->shellcmd? return a hash (not ref)
# of those params here.
sub _shellcmd_extra_params {
    my $self = shift;
    my @output_files = $self->_resolve_output_files;
    my @input_files = $self->_resolve_input_files;
    return (
        skip_if_output_is_present => 0,
        keep_dbh_connection_open => 0,
        input_files => \@input_files,
        output_files => \@output_files,
    );
}

sub build_cmdline_array_ref {
    my $self = shift;

    return [
        $self->tool_path,
        $self->_tool_param_args,
        $self->_tool_input_args,
    ];
}

# override with care.
sub execute {
    my $self = shift;

    $self->_validate_params;

    my %params = (
        cmd => $self->build_cmdline_array_ref,
        $self->_shellcmd_extra_params
    );

    return Genome::Sys->shellcmd(%params);
}


1;
