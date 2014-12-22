use strict;
use warnings;

use autodie qw(chmod);
use Fcntl ':mode';
use File::Temp qw();
use Test::More;

use Genome::Utility::File::Mode qw(mode);

my @bits = (
    S_IROTH,
    S_IWOTH,
    S_IXOTH,
    S_ISUID,
    S_ISGID,
    S_IRGRP,
    S_IWGRP,
    S_IXGRP,
    S_IRUSR,
    S_IWUSR,
    S_IXUSR,
);

my %bit = Genome::Utility::File::Mode::bit();
my %modes = (
    '00000' => 00000, '02000' => 02000, '04000' => 04000,
    '00100' => 00100, '00200' => 00200, '00400' => 00400,
    '00010' => 00010, '00020' => 00020, '00040' => 00040,
    '00001' => 00001, '00002' => 00002, '00004' => 00004,
);

my %multimodes = (
    (S_IRWXO) => {
        pass => [00007, 00057, 00077],
        fail => [00005, 00075],
    }
);

plan tests => scalar(keys %modes)
    + scalar(keys %bit)
    + scalar(map { @$_ } map { values %$_ } values %multimodes);

do {
    my $file = File::Temp->new();
    my $filename = $file->filename;

    for my $name (sort keys %modes) {
        my $expected = $modes{$name};
        subtest "mode=$name" => sub {
            plan tests => scalar(@bits);

            chmod $expected, $filename;

            my $got = mode($filename);
            for (@bits) {
                if ($expected == $_) {
                    ok($got->has_bit($_));
                } else {
                    ok(!$got->has_bit($_));
                }
            }
        };
    }
};

do {
    my $file = File::Temp->new();
    my $filename = $file->filename;

    for my $bitmask (sort keys %multimodes) {
        for my $mode (@{$multimodes{$bitmask}{pass}}) {
            chmod $mode, $filename;

            my $got = mode($filename);
            ok($got->has_bit($bitmask),
                sprintf("Bitmask '0%04o' passes for mode '0%04o'", $bitmask, $mode));
        }
        for my $mode (@{$multimodes{$bitmask}{fail}}) {
            chmod $mode, $filename;

            my $got = mode($filename);
            ok(!$got->has_bit($bitmask),
                sprintf("Bitmask '0%04o' fails for mode '0%04o'", $bitmask, $mode));
        }
    }
};

do {
    my $file = File::Temp->new();
    my $path = $file->filename;

    for my $name (keys %bit) {
        my $bit = $bit{$name};
        subtest "add/rm $name" => sub {
            plan tests => 4;

            my $initial_mode = 0;
            my $mode = mode($path);
            $mode->set_mode($initial_mode);
            is($mode->perm, $initial_mode, 'set initial permissions');

            my $add = 'add_' . $name;
            $mode->$add;
            ok($mode->has_bit($bit), $add);

            my $rm = 'rm_' . $name;
            $mode->$rm;
            ok(!$mode->has_bit($bit), $rm);

            is($mode->perm, $initial_mode, 'back to initial permissions');
        };
    };
};
