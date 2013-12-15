# Since we want to support OS X someday we should avoid taking advantage of
# case sensitive filesystems.

use strict;
use warnings;

use Git::Repository qw();
use Test::More;

my %count;
my $repo = Git::Repository->new();
my @files = $repo->run('ls-files');

plan tests => scalar(@files);

for my $path (@files) {
    my $key = lc($path);
    $count{$key}{count}++;
    push @{$count{$key}{paths}}, $path;
}

for my $path (@files) {
    my $key = lc($path);
    my $has_no_conflicts = ($count{$key}{count} <= 1);
    ok($has_no_conflicts, $path) or do {
        my @paths = @{$count{$key}{paths}};
        my @conflicts = grep { $_ ne $path } @paths;
        diag ' ' x 4, 'Conflicts with ', join(', ', @conflicts);
    };
}
