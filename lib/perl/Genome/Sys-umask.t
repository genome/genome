use strict;
use warnings;

use above "Genome";
use Fcntl ':mode';
use File::Spec qw();
use Genome::Utility::Test qw(abort run_ok);
use Genome::Utility::Test::Stat qw(has_bit hasnt_bit);

use Test::More tests => 1;

subtest 'create_directory overrides umask' => sub {
    plan tests => 7;

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

