use strict;
use warnings;

use Test::More;

use above "Genome";

use_ok('Genome::SoftwareResult'); #called with :: syntax below

class Genome::Bar {
    id => {
        is => 'Number',
    },
};

class Genome::Foo {
    is => 'Genome::SoftwareResult',
    has_input => {
        bar => {
            is => 'Genome::Bar',
            doc => 'to test lock name resolution with an object',
        },
        string => {
            is => 'Text',
            doc => 'to test lock name resolution with a string',
        },
    },
};

my ($foo_lock1, $bar_address1) = get_lock_name(1, "test1");
my ($foo_lock2, $bar_address2) = get_lock_name(1, "test1");

isnt($bar_address1, $bar_address2, "new object has different memory address");
is($foo_lock1, $foo_lock2, 'two locks for same inputs are the same');

my ($foo_lock3) = get_lock_name(2, "test1");
isnt($foo_lock1, $foo_lock3, "got different lock with different object input");

my ($foo_lock4) = get_lock_name(1, "test2");
isnt($foo_lock1, $foo_lock4, "got different lock with different string input");

done_testing();


sub get_lock_name {
    my $bar_id = shift;
    my $string = shift;

    my $bar = Genome::Bar->create(id => $bar_id);
    my $bar_memory_address = "".$bar;

    my $class = 'Genome::Foo';
    my $lookup_hash = $class->calculate_lookup_hash_from_arguments(bar => $bar, string => $string);
    my $foo_lock = $class->_resolve_lock_name($lookup_hash);
    ok($foo_lock, 'got a lock name');
    $bar->delete;

    return ($foo_lock, $bar_memory_address);
}
