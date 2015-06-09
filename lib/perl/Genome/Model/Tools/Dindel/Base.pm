package Genome::Model::Tools::Dindel::Base;

use strict;
use warnings;

use Genome;
use File::Spec;

class Genome::Model::Tools::Dindel::Base {
    is => 'Command',
    is_abstract => 1,
    has_input => [
        output_directory => {
            is => 'Path',
        }
    ],
};

sub help_synopsis {
    return <<EOS
EOS
}

sub help_detail {
    return <<EOS
EOS
}

use String::ShellQuote 'shell_quote';
sub shellcmd_arrayref {
    my ($self, %params) = @_;

    my $cmd = $params{cmd};
    unless (ref($cmd) eq 'ARRAY') {
        Carp::confess("cmd parameter expected to be an ARRAY reference not: " . ref($cmd));
    }
    $params{cmd} = join(' ', map {shell_quote $_} @$cmd);

    return Genome::Sys->shellcmd(%params);
}


sub create_output_directory {
    my $self = shift;
    if (! -d $self->output_directory) {
        Genome::Sys->create_directory($self->output_directory);
    }
    if (! -d $self->output_directory) {
        die "Couldn't create output_directory " . $self->output_directory;
    }
}

sub dindel_executable {
    my $self = shift;
    return '/usr/bin/dindel-tgi1.01-wu1';
}

sub python_script {
    my ($self, $name) = @_;
    return File::Spec->join($self->python_dir, $name . '.py');
}

sub python_dir {
    return '/usr/share/dindel-tgi1.01-wu1';
}

sub run_python_shellcmd {
    my $self = shift;

    # This is to get around Dindel's python scripts unconventional package
    # management strategy
    local $ENV{PYTHONPATH} = sprintf("%s:%s",
        File::Spec->join($self->python_dir, 'utils'),
        $ENV{PYTHONPATH} || '');

    my $rv = $self->shellcmd_arrayref(@_);
    return $rv;
}

1;
