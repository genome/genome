use strict;
use warnings;

package Genome::Utility::Test::Stat;
use parent 'Test::Builder::Module';

use Exporter 'import';
use Fcntl ':mode';
use List::MoreUtils qw(any);

use Genome::Utility::File::Mode qw();

our @EXPORT_OK = qw(has_bit hasnt_bit);

my %names = (
    S_IRWXO, 'other %s read, write, execute',
    S_IROTH, 'other %s read',
    S_IWOTH, 'other %s write',
    S_IXOTH, 'other %s execute',
    S_ISUID, '%s setuid',
    S_ISGID, '%s setgid',
    S_IRWXG, 'group %s read, write, execute',
    S_IRGRP, 'group %s read',
    S_IWGRP, 'group %s write',
    S_IXGRP, 'group %s execute',
    S_IRWXU, 'user %s read, write, execute',
    S_IRUSR, 'user %s read',
    S_IWUSR, 'user %s write',
    S_IXUSR, 'user %s execute',
);

sub has_bit {
    my ($mode, $bit, $name) = @_;
    unless ($name) {
        $name = sprintf($names{$bit}, 'has');
    }
    my $tb = __PACKAGE__->builder;
    return $tb->ok(Genome::Utility::File::Mode::_has_bit($mode, $bit), $name);
}

sub hasnt_bit {
    my ($mode, $bit, $name) = @_;
    unless ($name) {
        $name = sprintf($names{$bit}, 'does not have');
    }
    my $tb = __PACKAGE__->builder;
    return $tb->ok(!Genome::Utility::File::Mode::_has_bit($mode, $bit), $name);
}

1;
