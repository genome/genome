package Genome::Model::Tools::Joinx;

use strict;
use warnings;

use Genome;
use Carp qw/confess/;
use Sys::Hostname;

my $DEFAULT_VER = '1.8';
my $MINIMUM_VER_FOR_RLIB = 1.5;

class Genome::Model::Tools::Joinx {
    is  => 'Command',
    is_abstract => 1,
    has_input => [
        use_version => {
            is  => 'Version', 
            doc => "joinx version to be used.  default_value='$DEFAULT_VER'",
            is_optional   => 1, 
            default_value => $DEFAULT_VER,
        },
    ],
};


sub help_brief {
    "Tools to run joinx, a variant file manipulation tool.";
}

sub help_synopsis {
    "gmt joinx ...";
}

sub help_detail {                           
    "used to invoke joinx commands";
}

sub get_default_version { return $DEFAULT_VER }

sub joinx_path {
    my ($self, $ver) = @_;
    if(!defined $ver && ref $self) {
        $ver = $self->use_version || "";
    }
    my $path = "/usr/bin/joinx$ver";
    my $hostname = Sys::Hostname::hostname();
    if (! -x $path) {
        confess "Failed to find executable joinx version $ver at $hostname:$path!";
    }
    return $path;
}

sub rlib_path {
    my $self = shift;
    my $ver = $self->use_version;
    if($ver >= $MINIMUM_VER_FOR_RLIB) {
        my $path = "/usr/share/joinx$ver/";
        my $hostname = Sys::Hostname::hostname();
        if(! -d $path) {
            confess "Failed to find directory for joinx version $ver scripts at $hostname:$path!";
        }
        return $path;
    }
    else {
        confess "joinx version $ver does not install scripts.";
    }
}

1;

