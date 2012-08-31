package Genome::Model::Tools::Abyss;

use strict;
use warnings;

use Genome;

my %ABYSS_VERSIONS = (
    '1.2.7' => '/opt/abyss-1.2.7'
);

class Genome::Model::Tools::Abyss {
    is  => 'Command',
    is_abstract  => 1,
    has_input => [
        version => {
            is   => 'String',
            doc  => 'Abyss version.',
        },
    ],
};

sub help_brief {
    "Tools to run abyss, a parallel read assembler";
}

sub help_synopsis {
    my $self = shift;
    return <<"EOS"
gmt abyss ...
EOS
}

sub help_detail {
    return <<EOS
EOS
}

sub create {
    my $self = shift;
    my $obj = $self->SUPER::create(@_);
    if (!defined $obj->version or !exists $ABYSS_VERSIONS{$obj->version}) {
        die "Invalid abyss version: " . ($obj->version||"''") .
            ". Valid versions are: " . join(", ", keys %ABYSS_VERSIONS);
    }
    return $obj;
}

sub resolve_version {
    my $self = shift;

    my ($type) = ref($self) =~ /\:\:(\w+)$/;
    $type = 'velvet'.lc(substr $type, 0, 1);

    my $ver = $self->version;
    $ver = 'velvet_'.$ver unless $ver eq 'installed';
    
    my @uname = POSIX::uname();
    $ver .= '-64' if $uname[4] eq 'x86_64';
    
    my $exec = $ENV{GENOME_SW} . "/velvet/$ver/$type";
    unless (-x $exec) {
        $self->error_message("$exec is not excutable");
        return;
    }

    return $exec;
}

sub bindir {
    my $self = shift;
    my $basedir = $ABYSS_VERSIONS{$self->version};
    return join('/', $basedir, "bin");
}

1;
