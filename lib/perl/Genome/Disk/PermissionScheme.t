use strict;
use warnings;

use Test::More tests => 3;

use Genome;
use Genome::Utility::File::Mode qw(mode);

use File::Path qw(make_path);
use File::stat qw(stat);
use File::Temp qw();
use POSIX qw(getgroups);
use Test::Fatal qw(exception);

my $gid = (getgroups())[1]; # 0 = primary
unless (defined $gid) {
    die 'must have supplementary groups to run this test';
}
my $owner = Genome::Sys::User->get(username => Genome::Sys->username);
my $permission_scheme = Genome::Disk::PermissionScheme->create(
    gid => $gid,
    min_mode => oct(2550),
    max_mode => oct(7770),
    file_min_mask => oct(7111),
);

ok(exception { Genome::Disk::PermissionScheme->create(uid => 0) },
    'an exception is thrown when creating a permission scheme with a uid');

my $staging_path = File::Temp->newdir();
my $path = sub { File::Spec->join($staging_path, @_) };
my $perm = sub { mode(File::Spec->join($staging_path, @_))->mode };
my $gid_of = sub { stat(File::Spec->join($staging_path, @_))->gid };
setup_staging_path($staging_path->dirname);

subtest 'before apply' => sub {
    plan tests => 6;
    isnt($perm->('file'), $permission_scheme->file_mode($path->('file')));
    isnt($perm->('dir'), $permission_scheme->dir_mode($path->('dir')));
    isnt($perm->('dir', 'file'), $permission_scheme->file_mode($path->('dir', 'file')));
    isnt($gid_of->('file'), $gid);
    isnt($gid_of->('dir'), $gid);
    isnt($gid_of->('dir', 'file'), $gid);
};

$permission_scheme->apply($staging_path->dirname);

subtest 'after apply' => sub {
    plan tests => 6;
    is($perm->('file'), $permission_scheme->file_mode($path->('file')));
    is($perm->('dir'), $permission_scheme->dir_mode($path->('dir')));
    is($perm->('dir', 'file'), $permission_scheme->file_mode($path->('dir', 'file')));
    is($gid_of->('file'), $gid);
    is($gid_of->('dir'), $gid);
    is($gid_of->('dir', 'file'), $gid);
};

sub setup_staging_path {
    my $staging_path = shift;
    local $ENV{GENOME_SYS_UMASK} = 0;
    my $old_umask = umask 0;
    chmod oct(755), $staging_path;
    Genome::Sys->touch($path->('file'));
    make_path($path->('dir'));
    Genome::Sys->touch($path->('dir', 'file'));
    umask $old_umask;
}
