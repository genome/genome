use strict;
use warnings;

use autodie qw(chmod);
use Fcntl ':mode';
use File::Temp qw();
use Test::More;

use Genome::Utility::Test::Stat qw(has_bit hasnt_bit);

my %bits = (
#    S_IRWXO, 'o=rwx',
#    S_IRWXG, 'g=rwx',
#    S_IRWXU, 'u=rwx',
    S_IROTH, 'o=r',
    S_IWOTH, 'o=w',
    S_IXOTH, 'o=x',
    S_ISUID, 'setuid',
    S_ISGID, 'setgid',
    S_IRGRP, 'g=r',
    S_IWGRP, 'g=w',
    S_IXGRP, 'g=x',
    S_IRUSR, 'u=r',
    S_IWUSR, 'u=w',
    S_IXUSR, 'u=x',
);

my %modes = (
    '00000' => 00000, '02000' => 02000, '04000' => 04000,
    '00100' => 00100, '00200' => 00200, '00400' => 00400,
    '00010' => 00010, '00020' => 00020, '00040' => 00040,
    '00001' => 00001, '00002' => 00002, '00004' => 00004,
);
plan tests => scalar(keys %modes);

my $file = File::Temp->new();
my $filename = $file->filename;

for my $name (sort keys %modes) {
    my $expected = $modes{$name};
    subtest "mode=$name" => sub {
        my @bits = sort keys %bits;
        plan tests => scalar(@bits);

        chmod $expected, $filename;

        my $got = [stat($filename)]->[2];
        for (@bits) {
            if ($expected == $_) {
                has_bit($got, $_, q(has ) . $bits{$_});
            } else {
                hasnt_bit($got, $_, q(doesn't have ) . $bits{$_});
            }
        }
    };
}
