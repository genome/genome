#!/usr/bin/env genome-perl

use strict;
use warnings FATAL => qw(all);

use File::Copy      qw(copy);
use File::Temp      qw();
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
    install_hook($remote, $hook_script);
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
    $local->run_exit_is(1, push => 'origin', 'master');
}

sub test_modify_readme {
    my ($local, $remote) = @_;

    dir($local->work_tree)->file('README')->openw('New readme documentation!');

    commit_all($local, 'Modified the README file');
    $local->run_exit_ok(push => 'origin', 'master');
}




sub initialize_repo {
    my ($dir_spec) = @_;

    my $work_dir = dir( File::Temp->newdir() );
    while (my ($path, $content) = each %$dir_spec) {
        my $f = $work_dir->file($path);
        $f->dir->mkpath;
        $f->touch;
        $f->openw->say($content);
    }

    Git::Repository->run(init => $work_dir);
    my $git = Git::Repository->new( work_tree => $work_dir );
    commit_all($git, 'Initial commit');

    return $git;
}

sub bare_clone_repo {
    my ($origin) = @_;

    my $new_git_dir = dir( File::Temp->newdir() );
    Git::Repository->run(clone => '--bare', $origin->git_dir, $new_git_dir);

    return Git::Repository->new( git_dir => $new_git_dir );
}

sub clone_repo {
    my ($origin) = @_;

    my $new_work_tree = dir( File::Temp->newdir() );
    Git::Repository->run(clone => $origin->git_dir, $new_work_tree);

    return Git::Repository->new( work_tree => $new_work_tree );
}

sub commit_all {
    my ($git, $message) = @_;
    $git->run('add', '-A');
    $git->run('commit', '-am', $message);
}

sub install_hook {
    my ($git, $hook_script) = @_;
    my $install_path = dir($git->git_dir)->subdir('hooks')->file($hook_script->basename);

    copy($hook_script, $install_path) or die "file copy failed $!";
    chmod 0755, $install_path;

    return $install_path;
}
