use strict;
use warnings;

use above 'Genome';
use Test::More tests => 6;

use Sub::Install qw(reinstall_sub);

my $type = UR::Object::Type->define(
    class_name => 'TestResult',
    is => 'Genome::SoftwareResult',
);
ok($type, 'defined subclass of Genome::SoftwareResult: TestResult');

reinstall_sub({
    into => 'TestResult',
    as => 'create',
    code => sub {
        return UR::Object::create(@_);
    },
});

subtest 'lookup_hash set on create with test_name' => sub {
    plan tests => 1;
    my $tr = TestResult->create(test_name => 'test');
    is($tr->lookup_hash, $tr->calculate_lookup_hash, 'lookup_hash matches');
};

subtest 'lookup_hash set when test_name added' => sub {
    plan tests => 4;
    my $tr = TestResult->create();
    my $orig_lookup_hash = $tr->lookup_hash;
    ok(!$tr->test_name, 'no test_name initially');
    $tr->test_name('test');
    ok($tr->test_name, 'test_name added');
    isnt($tr->lookup_hash, $orig_lookup_hash, 'lookup_hash changed');
    is($tr->lookup_hash, $tr->calculate_lookup_hash, 'lookup_hash matches');
};

subtest 'lookup_hash set when test_name modified' => sub {
    plan tests => 3;
    my $tr = TestResult->create(test_name => 'test');
    my $orig_lookup_hash = $tr->lookup_hash;
    ok($tr->test_name('test2'), 'test_name modified');
    isnt($tr->lookup_hash, $orig_lookup_hash, 'lookup_hash changed');
    is($tr->lookup_hash, $tr->calculate_lookup_hash, 'lookup_hash matches');
};

subtest 'lookup_hash set when test_name set to undef' => sub {
    plan tests => 2;
    my $tr = TestResult->create(test_name => 'test');
    my $orig_lookup_hash = $tr->lookup_hash;
    $tr->test_name(undef);
    isnt($tr->lookup_hash, $orig_lookup_hash, 'lookup_hash changed');
    is($tr->lookup_hash, $tr->calculate_lookup_hash, 'lookup_hash matches');
};

subtest 'lookup_hash set when test_name removed' => sub {
    plan tests => 2;
    my $tr = TestResult->create(test_name => 'test');
    my $orig_lookup_hash = $tr->lookup_hash;
    $tr->remove_test_name();
    isnt($tr->lookup_hash, $orig_lookup_hash, 'lookup_hash changed');
    is($tr->lookup_hash, $tr->calculate_lookup_hash, 'lookup_hash matches');
};
