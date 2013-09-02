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

sub dindel_install_dir {
    my $self = shift;
    return "/gscmnt/gc2146/info/medseq/dindel";
}

sub dindel_executable {
    my $self = shift;
    return File::Spec->join($self->dindel_install_dir, "binaries", "dindel-1.01-linux-64bit");
}

sub python_script {
    my ($self, $name) = @_;
    my $path = File::Spec->join($self->dindel_install_dir, 'dindel-1.01-python', $name . '.py');
    if (-s $path) {
        return $path;
    } else {
        die "Couldn't find Dindel python script named '$name' at: $path\n";
    }
}

1;
