use strict;
use warnings;

use above "Genome";
use Fcntl ':mode';
use File::Spec qw();
use Genome::Utility::Test qw(abort run_ok);
use Genome::Utility::Test::Stat qw(has_bit hasnt_bit);

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
    my $cd_mode = [stat($cd_path)]->[2];
    has_bit($cd_mode, S_IWGRP, 'subdirectory made with Genome::Sys->create_directory has group write permissions');

    # verify mkdir, without overrides create_directory has, does not
    my $mkdir_path = File::Spec->join($td_path, 'mkdir');
    mkdir $mkdir_path;
    ok(-d $mkdir_path, 'made a subdirectory');
    my $mkdir_mode = [stat($mkdir_path)]->[2];
    hasnt_bit($mkdir_mode, S_IWGRP, 'subdirectory made with mkdir does not have group write permissions');
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
        my $cd_mode = [stat($cd_path)]->[2];
        hasnt_bit($cd_mode, S_ISUID);
        has_bit  ($cd_mode, S_ISGID);
        has_bit  ($cd_mode, S_IRUSR);
        has_bit  ($cd_mode, S_IWUSR);
        has_bit  ($cd_mode, S_IXUSR);
        has_bit  ($cd_mode, S_IRGRP);
        hasnt_bit($cd_mode, S_IWGRP);
        has_bit  ($cd_mode, S_IXGRP);
        hasnt_bit($cd_mode, S_IROTH);
        hasnt_bit($cd_mode, S_IWOTH);
        hasnt_bit($cd_mode, S_IXOTH);
    };

    # verify mkdir, without overrides create_directory has, does not
    my $mkdir_path = File::Spec->join($td_path, 'mkdir');
    mkdir $mkdir_path;
    ok(-d $mkdir_path, 'made a subdirectory with mkdir');

    subtest 'mkdir permissions' => sub {
        plan tests => 11;
        # default dir perms = 0777
        # with setgid and umask = 2700
        my $mkdir_mode = [stat($mkdir_path)]->[2];
        hasnt_bit($mkdir_mode, S_ISUID);
        has_bit  ($mkdir_mode, S_ISGID);
        has_bit  ($mkdir_mode, S_IRUSR);
        has_bit  ($mkdir_mode, S_IWUSR);
        has_bit  ($mkdir_mode, S_IXUSR);
        hasnt_bit($mkdir_mode, S_IRGRP);
        hasnt_bit($mkdir_mode, S_IWGRP);
        hasnt_bit($mkdir_mode, S_IXGRP);
        hasnt_bit($mkdir_mode, S_IROTH);
        hasnt_bit($mkdir_mode, S_IWOTH);
        hasnt_bit($mkdir_mode, S_IXOTH);
    };
};
