package Genome::Config;

use strict;
use warnings;

use File::Find::Rule qw();
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

    my @errors = $spec->validate($value);
    if (@errors) {
        my $msg = $spec->validation_error(@errors);
        croakf($msg);
    }

    if ($spec->sticky) {
        unless ($spec->has_env) {
            croakf('sticky values must have `env` set');
        }
        $ENV{$spec->env} = $value;
    }

    return $value;
}

sub spec {
    my $key = shift;
    my $subpath = Path::Class::File->new('genome', $key . '.yaml');
    my $file = (_lookup_files($subpath, global_dirs()))[0];
    unless (defined $file) {
        croakf('unable to locate spec: %s', $key);
    }
    return Genome::ConfigSpec->new_from_file($file);
}

sub all_specs {
    my @genome_dirs = map { Path::Class::Dir->new($_, 'genome') } global_dirs();
    my @spec_files = File::Find::Rule->file()
                                     ->name('*.yaml')
                                     ->not(File::Find::Rule->new->name('config.yaml'))
                                     ->in(@genome_dirs);
    my @specs = map { Genome::ConfigSpec->new_from_file($_) } @spec_files;

    my %specs;
    for my $spec (@specs) {
        next if exists $specs{$spec->key};
        $specs{$spec->key} = $spec;
    }

    return values %specs;
}

sub config_subpath { Path::Class::File->new('genome', 'config.yaml') }

sub local_values_enabled {
    return ! ! $ENV{XGENOME_ENABLE_USER_CONFIG};
}

sub home_dir {
    my $path = ( $ENV{XGENOME_CONFIG_HOME} || File::Spec->join($ENV{HOME}, '.config') );
    return Path::Class::Dir->new($path);
}

sub snapshot_dir {
    my @path = File::Spec->splitdir(__FILE__);
    my @chop = ('lib', 'perl', split('::', __PACKAGE__));
    my @chopped = splice(@path, -1 * @chop);
    my @base_dir = @path;
    return File::Spec->join(@base_dir, 'etc');
}

sub global_dirs {
    my $dirs = $ENV{XGENOME_CONFIG_DIRS} || '/etc';
    return map { Path::Class::Dir->new($_) } (snapshot_dir(), split(/:/, $dirs));
}

sub _lookup_value {
    my $spec = shift;
    my $value = _lookup_local_value($spec);
    if (defined $value) {
        return $value;
    }
    return _lookup_global_value($spec);
}

sub _lookup_local_value {
    my $spec = shift;

    return unless local_values_enabled();

    my $config_subpath = config_subpath();
    if ($spec->has_env && exists $ENV{$spec->env}) {
        return $ENV{$spec->env};
    }
    my @files = _lookup_files($config_subpath, home_dir());
    return _lookup_value_from_files($spec, @files);
}

sub _lookup_global_value {
    my $spec = shift;
    my $config_subpath = config_subpath();
    my @files = _lookup_files($config_subpath, global_dirs());
    return _lookup_value_from_files($spec, @files);
}

sub _lookup_value_from_files {
    my ($spec, @files) = @_;
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
    my @dirs = @_;
    my @files = map { Path::Class::File->new($_, $subpath) } @dirs;
    my @matches = grep { -f $_ } @files;
    return @matches;
}

1;
