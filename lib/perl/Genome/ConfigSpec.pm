package Genome::ConfigSpec;

use strict;
use warnings;

use Mouse;

use File::Basename qw(fileparse);
use Genome::Carp qw(croakf);
use List::Util qw(first);
use Module::Runtime qw();
use Path::Class qw();
use Try::Tiny qw(try catch);
use YAML::Syck qw();

require Genome::ConfigValidator::Defined;

has 'key'        => (is => 'ro', isa => 'Str', required => 1);
has 'validators' => (is => 'ro', isa => 'ArrayRef', default => sub { [] });

has 'default_value' => (is => 'ro', isa => 'Str', predicate => 'has_default_value');
has 'default_from'  => (is => 'ro', isa => 'Str', predicate => 'has_default_from');
has 'env'           => (is => 'ro', isa => 'Str', predicate => 'has_env');
has 'sticky'        => (is => 'ro', isa => 'Bool');

sub BUILD {
    my $self = shift;

    if ($self->sticky && !$self->has_env) {
        croakf('`sticky` requires `env`');
    }

    if ($self->has_default_value && $self->has_default_from) {
        croakf('`default_value` and `default_from` are mutually exclusive');
    }

    unshift @{$self->validators}, Genome::ConfigValidator::Defined->new();
}

sub new_from_file {
    my $class = shift;
    my $file = Path::Class::File->new(shift);

    my $key = fileparse($file, qr/\.yaml/i);
    my $data = YAML::Syck::LoadFile($file->stringify);

    my @validators;
    for my $v (@{$data->{validators}}) {
        my $module_name = join('::', 'Genome', 'ConfigValidator', ucfirst($v));
        try { Module::Runtime::require_module($module_name) }
        catch { croakf('failed to load validator: %s', $v) };
        push @validators, $module_name->new();
    }

    my %params = (
        key => $key,
        validators => \@validators,
    );
    for my $k (qw(env default_from default_value sticky)) {
        if (defined $data->{$k}) {
            $params{$k} = $data->{$k};
        }
    }
    return $class->new(%params);
}

sub validate {
    my ($self, $value) = @_;
    return first { $_->validate($value) } @{$self->validators};
}

sub validation_error {
    my ($self, $error) = @_;

    unless (defined $error) {
        croakf('no errors to display');
    }

    return sprintf('%s must be %s', $self->key, $error->message);
}


1;

=pod

=head1 NAME

Genome::ConfigSpec

=head1 DESCRIPTION

Genome::ConfigSpec defines the properties a given configuration value should
have.  These properties are: C<key>, C<validators>, C<env>, C<default_value>,
and C<sticky>.

=head2 Required Properties

C<key> is the identifier for the configuration value.  When loaded from a file
the file's basename is used as the key.

=head2 Optional Properties

C<validators> are C<Genome::ConfigValidator> objects.  When loaded from a file
the basename of one or more C<Genome::ConfigValidator> is used.

C<env> is the name of the environment variable that the configuration value is
bound to.

C<default_value> is the default value a configuration value should take when
not otherwise specified.

C<default_from> can be used to delegate the default value to another
configuration variable.  If a configuration variable does not have a
C<default_value> then the Genome::Config API will use C<default_from> to try
looking up a default value for the specified configuration variable.

C<sticky> is an attribute that suggests the value should not be allowed to
change.

=head 1 EXAMPLE

For example, if you want to register a configuration variable for the key,
'foo', you could create C<foo.yaml>,

  ---
  default_value: bar

=cut
