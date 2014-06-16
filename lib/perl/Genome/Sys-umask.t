use strict;
use warnings;

use above "Genome";

use Fcntl ':mode';
use File::Spec qw();

use Genome::Utility::Test qw(abort run_ok);

use Genome::Utility::File::Mode qw(mode);

use Test::More tests => 2;

subtest 'create_directory overrides umask' => sub {
    plan tests => 7;

    # use a string because we get a string from shell, not an octal
    local $ENV{GENOME_SYS_UMASK} = '0002';

    # setup
    my $td_path = File::Temp->newdir();
    ok(-d $td_path, 'made a temp directory to work in') or abort;
    run_ok(['chmod', 2775, $td_path], 'chmod that directory to have gid sticky') or abort;
    ok(umask oct(77), 'set umask so that no group permissions are allowed') or abort;

    # make sure create_directory overrides umask
    my $cd_path = File::Spec->join($td_path, 'cd');
    Genome::Sys->create_directory($cd_path);
    ok(-d $cd_path, 'made a subdirectory');
    my $cd_mode = mode($cd_path);
    ok($cd_mode->is_group_writable, 'subdirectory made with Genome::Sys->create_directory has group write permissions');

    # verify mkdir, without overrides create_directory has, does not
    my $mkdir_path = File::Spec->join($td_path, 'mkdir');
    mkdir $mkdir_path;
    ok(-d $mkdir_path, 'made a subdirectory');
    my $mkdir_mode = mode($mkdir_path);
    ok(!$mkdir_mode->is_group_writable, 'subdirectory made with mkdir does not have group write permissions');
};

subtest "GENOME_SYS_UMASK='0027'" => sub {
    # setting GENOME_SYS_UMASK='0027' revealed that by not casting
    # GENOME_SYS_UMASK as octal we were not using the expected umask
    plan tests => 7;

    # use a string because we get a string from shell, not an octal
    local $ENV{GENOME_SYS_UMASK} = '0027';

    # setup
    my $td_path = File::Temp->newdir();
    ok(-d $td_path, 'made a temp directory to work in') or abort;

    my $mode = oct(2775);
    ok((chmod $mode, $td_path), 'chmod that directory to have gid sticky') or abort;

    my $umask = oct(77);
    ok(umask $umask, 'set umask so that no group permissions are allowed') or abort;

    my $cd_path = File::Spec->join($td_path, 'cd');
    Genome::Sys->create_directory($cd_path);
    ok(-d $cd_path, 'made a subdirectory with create_directory');

    subtest 'create_directory permissions' => sub {
        plan tests => 11;
        # default dir perms = 0777
        # with setgid and GENOME_SYS_UMASK = 2750
        my $cd_mode = mode($cd_path);
        ok(!$cd_mode->is_setuid);
        ok( $cd_mode->is_setgid);
        ok( $cd_mode->is_user_readable);
        ok( $cd_mode->is_user_writable);
        ok( $cd_mode->is_user_executable);
        ok( $cd_mode->is_group_readable);
        ok(!$cd_mode->is_group_writable);
        ok( $cd_mode->is_group_executable);
        ok(!$cd_mode->is_other_readable);
        ok(!$cd_mode->is_other_writable);
        ok(!$cd_mode->is_other_executable);
    };

    # verify mkdir, without overrides create_directory has, does not
    my $mkdir_path = File::Spec->join($td_path, 'mkdir');
    mkdir $mkdir_path;
    ok(-d $mkdir_path, 'made a subdirectory with mkdir');

    subtest 'mkdir permissions' => sub {
        plan tests => 11;
        # default dir perms = 0777
        # with setgid and umask = 2700
        my $mkdir_mode = mode($mkdir_path);
        ok(!$mkdir_mode->is_setuid);
        ok( $mkdir_mode->is_setgid);
        ok( $mkdir_mode->is_user_readable);
        ok( $mkdir_mode->is_user_writable);
        ok( $mkdir_mode->is_user_executable);
        ok(!$mkdir_mode->is_group_readable);
        ok(!$mkdir_mode->is_group_writable);
        ok(!$mkdir_mode->is_group_executable);
        ok(!$mkdir_mode->is_other_readable);
        ok(!$mkdir_mode->is_other_writable);
        ok(!$mkdir_mode->is_other_executable);
    };
};
