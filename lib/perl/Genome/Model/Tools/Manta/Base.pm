package Genome::Model::Tools::Manta::Base;

use strict;
use warnings;

use Genome;
use Carp qw(confess);

my %TOOL_VERSIONS = (
    '0.29.6'     => '/gscmnt/gc13001/info/model_data/jwalker_scratch/src/manta-0.29.6.centos5_x86_64/bin',
);

class Genome::Model::Tools::Manta::Base {
    is => 'Command',
    is_abstract => 1,
    attributes_have => [
        tool_arg_name => {
            is => 'String',
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
            doc => 'The version of Manta to use.',
            valid_values => [sort keys %TOOL_VERSIONS],
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
    "Manta calls structural variants (SVs) and indels from mapped paired-end sequencing reads."
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
    return File::Spec->join($self->path_for_version($self->version),$self->_tool_subcommand_name);
}

# tool subcommand name as a string (e.g., 'configManta.py')
sub _tool_subcommand_name {
    confess "sub _tool_subcommand_name must be overridden in child class";
}

sub _resolve_output_files {
    # example:
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
    return '--'. $name if ($type eq 'Boolean');
    $value = qq{"$value"} if $type eq 'Text';
    return sprintf('--%s %s', $name, $value);
}

# given a property meta object (defining property name, and ostensibly
# tool_arg_name), return the list of cmdline args that should be passed
# to the tool based on the values set on the current object.
#
# this is where is_many properties get handled.
sub _tool_args_from_meta {
    my ($self, $meta) = @_;

    my $ur_name = $meta->property_name;
    # using // instead of ? : below made a lims perl test fail
    my $type = defined $meta->data_type ? $meta->data_type : 'String';
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
    my @args = map {$self->_tool_args_from_meta($_)} @metas;

    return @args;
}

sub _tool_input_args {
    my $self = shift;

    my @metas = $self->_tool_input_metas;
    my @args = map {$self->_tool_args_from_meta($_)} @metas;

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
        set_pipefail => 0,
        input_files => \@input_files,
        output_files => \@output_files,
    );
}

sub build_cmdline_list {
    my $self = shift;

    return (
        $self->tool_path,
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
