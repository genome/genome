# Since we want to support OS X someday we should avoid taking advantage of
# case sensitive filesystems.

use strict;
use warnings;

use Test::More;
use IPC::System::Simple qw(capture);

my %count;
my @files = git('ls-files');

plan tests => scalar(@files) + 1;

ok(scalar(@files) > 1, 'found more than one file');

for my $path (@files) {
    my $key = lc($path);
    $count{$key}{count}++;
    push @{$count{$key}{paths}}, $path;
}

for my $path (@files) {
    my $key = lc($path);
    my $has_no_conflicts = ($count{$key}{count} <= 1);
    ok($has_no_conflicts, 'File: ' . $path) or do {
        my @paths = @{$count{$key}{paths}};
        my @conflicts = grep { $_ ne $path } @paths;
        diag ' ' x 4, 'Conflicts with ', join(', ', @conflicts);
    };
}

sub git {
    my @args = shift;
    my @output = capture('git', @args);
    chomp @output;
    return @output;
}
