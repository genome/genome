use strict;
use warnings;

use above 'Genome';
use Test::More tests => 8;

use Genome::Utility::NamedArgs qw(named_args);

my %args;
my %base_args = (
    args => [
        first_name => 'Nathan',
        last_name  => 'Nutter',
        gender     => 'male',
        age        => 30,
    ],
    required => [
        'first_name',
        'last_name',
    ],
    optional => [
        'gender',
        'age',
    ],
);

$@ = '';

%args = %base_args;
$args{invalid_key} = 'some_value';
eval { named_args(%args) };
like($@, qr/^unexpected named argument/,
    'errored when invalid key passed to named_args');
$@ = '';

%args = %base_args;
delete $args{args};
eval { named_args(%args) };
like($@, qr/^missing args/,
    'errored when missing args');
$@ = '';

%args = %base_args;
delete $args{required};
delete $args{optional};
eval { named_args(%args) };
like($@, qr/^missing argument specification/,
    'errored when missing both required and optional');
$@ = '';

%args = %base_args;
delete $args{required};
$args{args} = [
    gender      => 'male',
    age         => 30,
];
eval { named_args(%args) };
is($@, '', 'did not error when optional are specified');
$@ = '';

%args = %base_args;
delete $args{optional};
$args{args} = [
    first_name  => 'Nathan',
    last_name   => 'Nutter',
];
eval { named_args(%args) };
is($@, '', 'did not error when required are specified');
$@ = '';

%args = %base_args;
$args{args} = [
    last_name  => 'Nutter',
    gender     => 'male',
    age        => 30,
];
eval { named_args(%args) };
like($@, qr/^missing required argument.*first_name/,
    'errored when missing required argument');
$@ = '';

%args = %base_args;
$args{args} = [
    first_name => 'Nathan',
    last_name  => 'Nutter',
    gender     => 'male',
];
eval { named_args(%args) };
is($@, '', 'did not error when missing optional argument');
$@ = '';

%args = %base_args;
$args{args} = [
    first_name  => 'Nathan',
    last_name   => 'Nutter',
    gender      => 'male',
    age         => 30,
    invalid_key => 'some_value',
];
eval { named_args(%args) };
like($@, qr/^unexpected named argument.*invalid_key/,
    'errored when invalid argument passed in args');
$@ = '';
