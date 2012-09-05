package Genome::Model::Tools::Sx::ExternalCmdBase;

use strict;
use warnings;

use Genome;

require Cwd;

class Genome::Model::Tools::Sx::ExternalCmdBase {
    is => 'Genome::Model::Tools::Sx::Base',
    is_abstract => 1,
    has => [
        _cmds => {
            is => 'Array',
            is_optional => 1,
            is_transient => 1,
            default_value => [],
            doc => 'Commands that were run.',
        },
    ],
};

sub _tmpdir {
    my $self = shift;
    $self->{_tmpdir} = Genome::Sys->create_temp_directory if not $self->{_tmpdir};
    return $self->{_tmpdir};
}

sub _rm_tmpdir {
    my $self = shift;
    return 1 if not $self->{_tmpdir};
    return File::Path::rmtree( $self->{_tmpdir} );
}

sub executable_path {
    my $self = shift;

    if ( not $self->version ) {
        $self->error_message('No version to get executable path!');
        return;
    }

    my %versions = $self->_cmd_versions;
    if ( not $versions{ $self->version } ) {
        $self->error_message('Invalid version! '.$self->version);
        return;
    }

    return $versions{ $self->version };
}

sub cmd_property_names {
    my $self = shift;
    my %cmd_properties = $self->_cmd_properties;
    return sort keys %cmd_properties;
}

sub build_command {
    my $self = shift;

    my $cmd = $self->executable_path;
    return if not $cmd;

    my $meta = $self->__meta__;
    for my $key ( $self->cmd_property_names ) {
        my $property = $meta->property_meta_for_name($key);
        my $value = $self->$key;
        next if not defined $value;
        $key =~ s/\_/\-/g;
        $cmd .= sprintf(
            ' %s%s%s',
            ( length($key) == 1 ? '-' : '--'),                       # - or --
            $key,                                                    # param name
            ( $property->data_type eq 'Boolean' ? '' : ' '.$value ), # value or empty string for boolean
        );
    }

    return $cmd;
}

sub _run_command {
    my ($self, $cmd) = @_;

    my $cmds = $self->_cmds;
    push @$cmds, $cmd;
    $self->_cmds($cmds);

    my $cwd = Cwd::getcwd();
    chdir $self->_tmpdir;

    my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
    chdir $cwd;

    if ( not $rv ) {
        $self->error_message($@) if $@;
        $self->error_message("Failed to run: $cmd");
        return;
    }

    return 1;
}

1;

