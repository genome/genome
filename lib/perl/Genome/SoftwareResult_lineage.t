use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;
use feature ':5.10';
use Sub::Install qw(reinstall_sub);

my @results;
subtest 'setup' => sub {
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

    for my $n (0..9) {
        my $tr = TestResult->create();
        ok($tr, "created \$results[$n]");
        push @results, $tr;
    }

    # Tree: 0 1
    #        2
    #      3 4 5
    #          6
    #        7
    #       8 9
    $results[0]->add_user(label => 'uses', user => $results[2]);
    $results[1]->add_user(label => 'uses', user => $results[2]);
    $results[2]->add_user(label => 'uses', user => $results[3]);
    $results[2]->add_user(label => 'uses', user => $results[4]);
    $results[2]->add_user(label => 'uses', user => $results[5]);
    $results[5]->add_user(label => 'uses', user => $results[6]);
    $results[3]->add_user(label => 'uses', user => $results[7]);
    $results[4]->add_user(label => 'uses', user => $results[7]);
    $results[6]->add_user(label => 'uses', user => $results[7]);
    $results[7]->add_user(label => 'uses', user => $results[8]);
    $results[7]->add_user(label => 'uses', user => $results[9]);
};

subtest 'parents' => sub {
    plan tests => (scalar @results);
    my %cases = (
        0 => [],
        1 => [],
        2 => [@results[0, 1]],
        3 => [$results[2]],
        4 => [$results[2]],
        5 => [$results[2]],
        6 => [$results[5]],
        7 => [@results[3, 4, 6]],
        8 => [$results[7]],
        9 => [$results[7]],
    );
    for my $n (sort keys %cases) {
        my $exp = [sort { $a->id cmp $b->id } @{$cases{$n}}];
        my $got = [sort { $a->id cmp $b->id } $results[$n]->parents];
        is_deeply($got, $exp, "got expected parents for \$results[$n]");
    }
};

subtest 'children' => sub {
    plan tests => (scalar @results);
    my %cases = (
        0 => [$results[2]],
        1 => [$results[2]],
        2 => [@results[3, 4, 5]],
        3 => [$results[7]],
        4 => [$results[7]],
        5 => [$results[6]],
        6 => [$results[7]],
        7 => [@results[8, 9]],
        8 => [],
        9 => [],
    );
    for my $n (sort keys %cases) {
        my $exp = [sort { $a->id cmp $b->id } @{$cases{$n}}];
        my $got = [sort { $a->id cmp $b->id } $results[$n]->children];
        is_deeply($got, $exp, "got expected children for \$results[$n]");
    }
};

subtest 'ancestors' => sub {
    plan tests => (scalar @results);
    my %cases = (
        0 => [],
        1 => [],
        2 => [@results[0..1]],
        3 => [@results[0..2]],
        4 => [@results[0..2]],
        5 => [@results[0..2]],
        6 => [@results[0..2, 5]],
        7 => [@results[0..6]],
        8 => [@results[0..7]],
        9 => [@results[0..7]],
    );
    for my $n (sort keys %cases) {
        # can't guarantee sort order due to loss of dimension (tree => list)
        my $exp = [sort { $a->id cmp $b->id } @{$cases{$n}}];
        my $got = [sort { $a->id cmp $b->id } $results[$n]->ancestors];
        is_deeply($got, $exp, "got expected ancestors for \$results[$n]");
    }
};

subtest 'descendents' => sub {
    plan tests => (scalar @results);
    my %cases = (
        0 => [@results[2..9]],
        1 => [@results[2..9]],
        2 => [@results[3..9]],
        3 => [@results[7..9]],
        4 => [@results[7..9]],
        5 => [@results[6..9]],
        6 => [@results[7..9]],
        7 => [@results[8..9]],
        8 => [],
        9 => [],
    );
    for my $n (sort keys %cases) {
        # can't guarantee sort order due to loss of dimension (tree => list)
        my $exp = [sort { $a->id cmp $b->id } @{$cases{$n}}];
        my $got = [sort { $a->id cmp $b->id } $results[$n]->descendents];
        is_deeply($got, $exp, "got expected descendents for \$results[$n]");
    }
};
