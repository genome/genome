package Genome::Model::Tools::Speedseq::Base;

use strict;
use warnings;

use Genome;
use Carp qw(confess);
use Data::Dumper;

class Genome::Model::Tools::Speedseq::Base {
    is => 'Command::V2',
    is_abstract => 1,
    attributes_have => [
        tool_param_name => {
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
    has_param => [
        version => {
            is => 'Version',
            doc => 'SpeedSeq version to be used',
            is_optional => 1,
        },
        config_file => {
            is => 'Text',
            doc => 'path to speedseq.config file (default: same directory as speedseq)',
            is_optional => 1,
            tool_param_name => 'K',
            tool_input_file => 1,
        },
        verbose => {
            is => 'Boolean',
            doc => 'verbose',
            is_optional => 1,
            tool_param_name => 'v',
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

sub help_detail {
    "A flexible framework for rapid genome analysis and interpretation"
}

my %TOOL_VERSIONS = (
    'test'     => '/gscmnt/sata849/info/speedseq_freeze/v1/speedseq/bin/speedseq',
);

sub available_versions {
    my $self = shift;
    return keys(%TOOL_VERSIONS);
}

sub path_for_version {
    my ($class, $version) = @_;
    unless ($version) {
        die('No version defined!');
    }
    my $path = $TOOL_VERSIONS{$version};
    return $path if defined $path;
    die 'No path found for version: '. $version;
}

sub tool_path {
    my $self = shift;
    return $self->path_for_version($self->version);
}

# tool subcommand name as a string (e.g., 'realign')
sub _tool_subcommand_name {
    confess "sub _tool_subcommand_name must be overridden in child class";
}

sub _resolve_output_files {
    # example:
}

sub _resolve_input_files {
    my $self = shift;

    my @metas = $self->_tool_input_file_metas;

    my @input_files = map { $self->_tool_input_files_from_meta($_) } @metas;
    $self->input_files(\@input_files);
    return @input_files;
}

# tool params should be 'params' and also have tool_param_name defined
sub _tool_param_metas {
    my $class = shift;
    my @param_property_metas = $class->__meta__->properties(is_param => 1);
    return grep {$_->can("tool_param_name") && $_->tool_param_name} @param_property_metas;
}

# tool inputs should be 'inputs' and also have tool_input_bare_arg_position defined
sub _tool_input_metas {
    my $class = shift;
    my @input_property_metas = $class->__meta__->properties(is_input => 1);
    my @tool_input_metas = grep {$_->can("tool_bare_arg_position") && $_->tool_bare_arg_position} @input_property_metas;
    return sort { $a->tool_bare_arg_position <=> $b->tool_bare_arg_position } @tool_input_metas;
}

sub _tool_input_file_metas {
    my $class = shift;
    my @property_metas = $class->__meta__->properties();
    return grep {$_->can("tool_input_file") && $_->tool_input_file} @property_metas;
}

# Handle conversion from boolean => 'true'/'false' and general tool arg
# formatting (argname=argvalue).
sub _format_tool_arg {
    my ($type, $name, $value) = @_;
    return '-'. $name if ($type eq 'Boolean');
    $value = qq{"$value"} if $type eq 'Text';
    return sprintf '-%s %s', $name, $value;
}

# given a property meta object (defining property name, and ostensibly
# tool_param_name), return the list of cmdline args that should be passed
# to the tool based on the values set on the current object.
#
# this is where is_many properties get handled.
sub _tool_param_args_from_meta {
    my ($self, $meta) = @_;

    my $ur_name = $meta->property_name;
    # using // instead of ? : below made a lims perl test fail
    my $type = defined $meta->data_type ? $meta->data_type : 'String';
    my $tool_param_name = $meta->tool_param_name;

    return unless defined $self->$ur_name;

    # This works for things that are is_many or not.
    return map {_format_tool_arg($type, $tool_param_name, $_)} $self->$ur_name;
}

sub _tool_input_args_from_meta {
    my ($self, $meta) = @_;

    my $ur_name = $meta->property_name;
    
    return unless defined $self->$ur_name;

    # This works for things that are is_many or not.
    return map {$_} $self->$ur_name;
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
# speedseq_param_name defined. Boolean values get translated to true/false,
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

sub _tool_param_args {
    my $self = shift;

    my @metas = $self->_tool_param_metas;
    my @args = map {$self->_tool_param_args_from_meta($_)} @metas;

    return @args;
}

sub _tool_input_args {
    my $self = shift;

    my @metas = $self->_tool_input_metas;
    my @args = map {$self->_tool_input_args_from_meta($_)} @metas;

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
        set_pipefail => 0,
        input_files => \@input_files,
        output_files => \@output_files,
    );
}

sub build_cmdline_list {
    my $self = shift;

    return (
        $self->tool_path,
        $self->_tool_subcommand_name,
        $self->_tool_param_args,
        $self->_tool_input_args,
    );
}

sub build_cmdline_string {
    my $self = shift;

    my $cmd = join(' ', $self->build_cmdline_list);
    return $cmd;
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
