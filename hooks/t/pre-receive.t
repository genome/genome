#!/usr/bin/env genome-perl

use strict;
use warnings FATAL => qw(all);

use File::Copy      qw(copy);
use Git::Repository qw(Test);
use Path::Class     qw(file dir);

use Test::More;

my $hook_script = file('hooks/pre-receive');
my $proto_repo  = initialize_repo({
    'README' => 'contents of readme file',
    'ur/test' => 'test file inside UR dir'
});

my $setup = sub {
    my $remote = bare_clone_repo($proto_repo);
    $remote->install_hook($hook_script, 'pre-receive');
    my $local = clone_repo($remote);

    return ($local, $remote);
};

test_modify_ur_file( $setup->() );
test_modify_readme(  $setup->() );
done_testing();

sub test_modify_ur_file {
    my ($local, $remote) = @_;

    dir($local->work_tree)->subdir('ur')->file('test')->openw('Modified file contents');

    commit_all($local, 'Modified the test UR file');
    $local->run_exit_is(1, push => $remote->git_dir, 'master');
}

sub test_modify_readme {
    my ($local, $remote) = @_;

    dir($local->work_tree)->file('README')->openw('New readme documentation!');

    commit_all($local, 'Modified the README file');
    $local->run_exit_ok(push => $remote->git_dir, 'master');
}




sub initialize_repo {
    my ($dir_spec) = @_;

    my $git = Git::Repository->new_tmp_repo;
    while (my ($path, $content) = each %$dir_spec) {
        my $f = dir($git->work_tree)->file($path);
        $f->dir->mkpath;
        $f->touch;
        $f->openw->say($content);
    }

    commit_all($git, 'Initial commit');

    return $git;
}

sub bare_clone_repo {
    my ($origin) = @_;

    my $git = Git::Repository->new_tmp_repo('--bare');
    $git->run('fetch', $origin->git_dir, 'master:master');

    return $git;
}

sub clone_repo {
    my ($origin) = @_;

    # This involves a little hackery because fetching to the checked-out branch
    # of a non-bare repo does not work

    my $git = Git::Repository->new_tmp_repo;
    $git->run('checkout', '-b', 'fake');
    $git->run('fetch', $origin->git_dir, 'master:master');
    $git->run('checkout', 'master');

    return $git;
}

sub commit_all {
    my ($git, $message) = @_;
    $git->run('add', '-A');
    $git->run('commit', '-am', $message);
}

1;
