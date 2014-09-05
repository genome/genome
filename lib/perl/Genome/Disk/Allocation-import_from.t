use strict;
use warnings;

use Test::More tests => 2;

use Genome;
use Genome::Test::Factory::DiskAllocation qw();

use File::Path qw(make_path);
use Path::Class::Dir qw();

my $owner = Genome::Sys::User->get(username => Genome::Sys->username);
my $allocation = Genome::Test::Factory::DiskAllocation->setup_object(owner => $owner);

my $staging_path = File::Temp->newdir();
my $staging_path_for = sub { File::Spec->join($staging_path, @_) };
my $allocation_path_for = sub { File::Spec->join($allocation->absolute_path, @_) };;
setup_staging_path($staging_path->dirname);
$DB::single = 1;
my $rv = $allocation->import_from($staging_path->dirname);
ok($rv, 'import_from returned successful');
subtest 'allocation paths exist' => sub {
    plan tests => 3;
    ok(-f $allocation_path_for->('file'));
    ok(-d $allocation_path_for->('dir'));
    ok(-f $allocation_path_for->('dir', 'file'));
};

sub setup_staging_path {
    my $staging_path = shift;
    local $ENV{GENOME_SYS_UMASK} = 0;
    my $old_umask = umask 0;
    chmod oct(755), $staging_path;
    Genome::Sys->touch($staging_path_for->('file'));
    make_path($staging_path_for->('dir'));
    Genome::Sys->touch($staging_path_for->('dir', 'file'));
    umask $old_umask;
}
