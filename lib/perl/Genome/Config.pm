package Genome::Config;

use strict;
use warnings;

use Genome::Carp qw(croakf);
use Genome::ConfigSpec qw();
use Path::Class qw();
use YAML::Syck qw();

sub get {
    my $key = shift;
    my $spec = spec($key);

    my $value = _lookup_value($spec);
    if (!$spec->has_default_value && !defined($value)) {
        croakf('required config key missing: %s', $key);
    }

    if (!defined($value)) {
        $value = $spec->default_value;
    }

    for my $v (@{$spec->validators}) {
        if (my $error = $v->($value, $spec)) {
            croakf('%s must be %s', $key, $error);
        }
    }

    return $value;
}

sub spec {
    my $key = shift;
    my $subpath = Path::Class::File->new('genome', $key . '.yaml');
    my $file = _lookup_file($subpath);
    unless (defined $file) {
        croakf('unable to locate spec: %s', $key);
    }
    return Genome::ConfigSpec->new_from_file($file);
}

sub _lookup_value {
    my $spec = shift;

    if ($spec->has_env && exists $ENV{$spec->env}) {
        return $ENV{$spec->env};
    }

    my @files = _lookup_files(Path::Class::File->new('genome', 'config.yaml'));
    for my $f (@files) {
        my $data = YAML::Syck::LoadFile($f);
        if ($data->{$spec->key}) {
            return $data->{$spec->key};
        }
    }

    return;
}

sub _lookup_files {
    my $subpath = shift;
    my @dirs = (_home_dir(), _dirs());
    my @files = map { Path::Class::File->new($_, $subpath) } @dirs;
    my @matches = grep { -f $_ } @files;
    return @matches;
}

sub _lookup_file {
    my @matches = _lookup_files(@_);
    return $matches[0];
}

sub _home_dir {
    return $ENV{XGENOME_CONFIG_HOME} || '$HOME/.config';
}

sub _dirs {
    my $dirs = $ENV{XGENOME_CONFIG_DIRS} || '/etc';

    # TODO One way to support per-snapshot dir.  Maybe not needed.
    my @path = File::Spec->splitdir(__FILE__);
    my @chop = ('lib', 'perl', split('::', __PACKAGE__));
    my @chopped = splice(@path, -1 * @chop);
    my @base_dir = @path;

    return (File::Spec->join(@base_dir, 'etc'), split(/:/, $dirs));
}

1;
