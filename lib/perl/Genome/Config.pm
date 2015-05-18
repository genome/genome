package Genome::Config;

use strict;
use warnings;

use File::Find::Rule qw();
use Genome::Carp qw(croakf);
use Genome::ConfigSpec qw();
use Path::Class qw();
use YAML::Syck qw();

require Scope::Guard;

=item get()

C<get($key)> retrieves the configuration value specified by the input key.  It
does a hierarchical search for the value in the following locations:

- environment variable
- user config file
- snapshot config file
- global config files

User's config files are discovered in C<XGENOME_CONFIG_HOME> which defaults to
C<$HOME/.config>.  Global config files are discovered in C<XGENOME_CONFIG_DIRS>
which defaults to C</etc> and can be a colon-delimited list of directories like
C<PATH>.

Values are specified in C<genome/config.yaml> in one of the source directories
or via the corresponding environment variable.

Configuration values are "registered" by creating a spec file, see
L<Genome::ConfigSpec>.  Spec files go in the snapshot config directory, e.g.
C<etc/genome/foo.yaml>.

If you want the configuration value to be optional you must specify a
C<default_value>. The C<default_value> can be zero or an empty string but it
must be defined.

If the spec file enables the C<sticky> attribute then the C<env> attribute must
also be set so that the value can be stored in an environment variable.  Since
an environment variable is the primary source this will have the effect of
ignoring values from any other location.  Since environment variables are
global it cannot be guaranteed that it will not be overwritten.

=cut

sub get {
    my $key = shift;
    my $spec = spec($key);

    my $value = _lookup_value($spec);
    my $error = $spec->validate($value);
    if (defined $error) {
        my $msg = $spec->validation_error($error);
        croakf($msg);
    }

    if ($spec->sticky) {
        $ENV{$spec->env} = $value;
    }

    return $value;
}

sub validate {
    my $key = shift;
    my $spec = spec($key);
    my $value = _lookup_value($spec);
    return $spec->validate($value);
}

sub spec {
    my $key = shift;
    my $subpath = Path::Class::File->new('genome', $key . '.yaml');
    my $file = (_lookup_files($subpath, snapshot_dir()))[0];
    unless (defined $file) {
        croakf('unable to locate spec: %s', $key);
    }
    return Genome::ConfigSpec->new_from_file($file);
}

sub set_env {
    my ($key, $value) = @_;
    my $spec = Genome::Config::spec($key);
    unless ($spec->has_env) {
        croakf('configuration does not specify an environment variable: %s', $key);
    }

    my $env_key = $spec->env;

    my $guard_closure;
    if (exists $ENV{$env_key}) {
        my $orig_value = $ENV{$env_key};
        $guard_closure = sub { $ENV{$env_key} = $orig_value };
    }
    else {
        $guard_closure = sub { delete $ENV{$env_key} };
    }

    $ENV{$env_key} = $value;

    if (defined wantarray) {
        my $guard = Scope::Guard->new($guard_closure);
        return $guard;
    }

    return 1;
}

sub all_specs {
    my $snapshot_dir = Path::Class::Dir->new(snapshot_dir(), 'genome');
    my @spec_files = File::Find::Rule->file()
                                     ->name('*.yaml')
                                     ->not(File::Find::Rule->new->name('config.yaml'))
                                     ->in($snapshot_dir);
    my @specs = map { Genome::ConfigSpec->new_from_file($_) } @spec_files;

    my %specs;
    for my $spec (@specs) {
        next if exists $specs{$spec->key};
        $specs{$spec->key} = $spec;
    }

    return values %specs;
}

sub all_keys {
    return map { $_->key } all_specs();
}

sub config_subpath { Path::Class::File->new('genome', 'config.yaml') }

sub home_dir {
    my $path = ( $ENV{XGENOME_CONFIG_HOME} || File::Spec->join($ENV{HOME}, '.config') );
    return Path::Class::Dir->new($path);
}

sub snapshot_dir {
    if ($ENV{XGENOME_CONFIG_SNAP}) {
        return $ENV{XGENOME_CONFIG_SNAP};
    }

    my @path = File::Spec->splitdir(__FILE__);
    my @chop = ('lib', 'perl', split('::', __PACKAGE__));
    my @chopped = splice(@path, -1 * @chop);
    my @base_dir = @path;
    return File::Spec->join(@base_dir, 'etc');
}

sub global_dirs {
    my $dirs = $ENV{XGENOME_CONFIG_DIRS} || '/etc';
    return map { Path::Class::Dir->new($_) } split(/:/, $dirs);
}

sub _lookup_value {
    my $spec = shift;

    my $config_subpath = config_subpath();
    if ($spec->has_env && exists $ENV{$spec->env}) {
        return $ENV{$spec->env};
    }
    my @files = _lookup_files($config_subpath, home_dir());
    my $value = _lookup_value_from_files($spec, @files);
    if (defined $value) {
        return $value;
    }

    @files = _lookup_files($config_subpath, snapshot_dir(), global_dirs());
    $value = _lookup_value_from_files($spec, @files);

    if (!defined($value) && $spec->has_default_value) {
        $value = $spec->default_value;
    }

    return $value;
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

sub reference_sequence_directory {
    return join('/', Genome::Config::get('model_root'), 'reference_sequences');
}

1;
