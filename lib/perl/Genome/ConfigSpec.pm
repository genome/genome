package Genome::ConfigSpec;

use strict;
use warnings;

use Mouse;

use File::Basename qw(fileparse);
use Genome::Carp qw(croakf);
use Module::Runtime qw();
use Path::Class qw();
use Try::Tiny qw(try catch);
use YAML::Syck qw();

has 'key'        => (is => 'ro', isa => 'Str', required => 1);
has 'type'       => (is => 'ro', isa => 'Str', required => 1);
has 'validators' => (is => 'ro', isa => 'ArrayRef[CodeRef]', default => sub { [] });

has 'default_value' => (is => 'ro', isa => 'Str', predicate => 'has_default_value');
has 'env'           => (is => 'ro', isa => 'Str', predicate => 'has_env');
has 'sticky'        => (is => 'ro', isa => 'Bool');

sub BUILD {
    my $self = shift;

    if ($self->sticky && !$self->has_env) {
        croakf('`sticky` requires `env`');
    }
};

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

        my $coderef = $module_name->can('validate');
        unless ($coderef) {
            croakf('validator missing implementation: %s', $v);
        }
        push @validators, $coderef;
    }

    my %params = (
        key => $key,
        validators => \@validators,
    );
    for my $k (qw(env default_value type sticky)) {
        if (defined $data->{$k}) {
            $params{$k} = $data->{$k};
        }
    }
    return $class->new(%params);
}

1;

=pod

=head1 NAME

Genome::ConfigSpec

=head1 DESCRIPTION

Genome::ConfigSpec defines the properties a given configuration value should
have.  These properties are: C<key>, C<type>, C<validators>, C<env>,
C<default_value>, and C<sticky>.

=head2 Required Properties

C<key> is the identifier for the configuration value.  When loaded from a file
the file's basename is used as the key.

C<type> specifies the type of the value.  This may often be a UR::Value type.

=head2 Optional Properties

C<validators> are anonymous subs that validate values.  When loaded from a file
the basename of one or more C<Genome::ConfigValidator> is used.

C<env> is the name of the environment variable that the configuration value is
bound to.

C<default_value> is the default value a configuration value should take when
not otherwise specified.

C<sticky> is an attribute that suggests the value should not be allowed to
change.

=head 1 EXAMPLE

For example, if you want to register a configuration variable for the key,
'foo', you could create C<foo.yaml>,

  ---
  type: Text
  default_value: bar

=cut
