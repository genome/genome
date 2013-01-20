use strict;
use warnings;

use Test::More;
use above 'Genome';

BEGIN {
    use_ok 'Genome::Utility::Test', qw(sub_test compare_ok);
}

sub_test('compare_ok' => sub {
    my $a_fh = File::Temp->new(TMPDIR => 1);
    my $a_fn = $a_fh->filename;
    $a_fh->print("a\n");
    $a_fh->close();

    my $b_fh = File::Temp->new(TMPDIR => 1);
    my $b_fn = $b_fh->filename;
    $b_fh->print("b\n");
    $b_fh->close();

    my $aa_fh = File::Temp->new(TMPDIR => 1);
    my $aa_fn = $aa_fh->filename;
    $aa_fh->print("a\n"); # like a, not aa!
    $aa_fh->close();

    {
        my $compare_ok = compare_ok($a_fn, $b_fn, test => 0);
        my $diff    = (system(qq(diff -u "$a_fn" "$b_fn" > /dev/null)) == 0 ? 1 : 0);
        is($diff, 0, 'diff detected diff between different files');
        is($compare_ok, $diff, 'compare_ok detected diff between different files');
    }

    {
        my $compare_ok = compare_ok($a_fn, $aa_fn, test => 0);
        my $diff    = (system(qq(diff -u "$a_fn" "$aa_fn" > /dev/null)) == 0 ? 1 : 0);
        is($diff, 1, 'diff did not detect diff between similar files');
        is($compare_ok, $diff, 'compare_ok did not detect diff between similar files');
    }
});

done_testing();
