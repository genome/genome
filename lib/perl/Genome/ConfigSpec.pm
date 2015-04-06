package Genome::ConfigSpec;

use strict;
use warnings;

use Mouse;

use Genome::Carp qw(croakf);
use Module::Runtime qw();
use Path::Class qw();
use Try::Tiny qw(try catch);
use YAML::Syck qw();

has 'key'  => (is => 'ro', isa => 'Str', required => 1);
has 'type' => (is => 'ro', isa => 'Str', required => 1);
has 'validators'    => (is => 'ro', isa => 'ArrayRef[CodeRef]', default => sub { [] });

has 'default_value' => (is => 'ro', isa => 'Str', predicate => 'has_default_value');
has 'env'           => (is => 'ro', isa => 'Str', predicate => 'has_env');

sub new_from_file {
    my $class = shift;
    my $file = Path::Class::File->new(shift);

    (my $key = $file->basename) =~ s/\.yaml$//i;
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
    for my $k (qw(env default_value type)) {
        if (defined $data->{$k}) {
            $params{$k} = $data->{$k};
        }
    }
    return $class->new(%params);
}

1;
