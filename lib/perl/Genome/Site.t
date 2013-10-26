#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 5;

use Carp qw(croak);
use File::Path qw(make_path);
use File::Temp qw(); # tempdir is wrapped

BEGIN {
    use_ok('Genome::Site');
}

subtest 'shallow match' => sub {
    plan tests => 3;

    my $tmp_dir_path = tempdir();

    local %INC = %INC;
    local @INC = @INC;
    unshift @INC, $tmp_dir_path;

    my @site_dirs = qw(edu wustl genome);
    no warnings 'redefine';
    local *Genome::Site::site_dirs = sub { return @site_dirs };
    use warnings 'redefine';

    my $foo = Genome::Site::site_pkg($site_dirs[0]);
    my $foo_filename = Genome::Site::module_to_filename($foo);
    my $foo_path = create_module($tmp_dir_path, $foo,
        'package ' . $foo . ';',
        '1;',
    );
    ok(-s $foo_path, "$foo_filename has size");

    ok(!defined($INC{$foo_filename}), "$foo has not been loaded yet");
    eval 'use Genome::Site';
    ok(defined($INC{$foo_filename}), "$foo has been loaded");
};

subtest 'deep match' => sub {
    plan tests => 3;

    my $tmp_dir_path = tempdir();

    local %INC = %INC;
    local @INC = @INC;
    unshift @INC, $tmp_dir_path;

    my @site_dirs = qw(edu wustl genome);
    no warnings 'redefine';
    local *Genome::Site::site_dirs = sub { return @site_dirs };
    use warnings 'redefine';

    my $foo = Genome::Site::site_pkg(@site_dirs);
    my $foo_filename = Genome::Site::module_to_filename($foo);
    my $foo_path = create_module($tmp_dir_path, $foo,
        'package ' . $foo . ';',
        '1;',
    );
    ok(-s $foo_path, "$foo_filename has size");

    ok(!defined($INC{$foo_filename}), "$foo has not been loaded yet");
    eval 'use Genome::Site';
    ok(defined($INC{$foo_filename}), "$foo has been loaded");
};

subtest 'catch site compile errors' => sub {
    plan tests => 4;

    my $tmp_dir_path = tempdir();

    local %INC = %INC;
    local @INC = @INC;
    unshift @INC, $tmp_dir_path;

    my @site_dirs = qw(edu wustl genome);
    no warnings 'redefine';
    local *Genome::Site::site_dirs = sub { return @site_dirs };
    use warnings 'redefine';

    my $foo = Genome::Site::site_pkg($site_dirs[0]);
    my $foo_filename = Genome::Site::module_to_filename($foo);
    my $foo_path = create_module($tmp_dir_path, $foo,
        'package ' . $foo . ';',
        '0;',
    );
    ok(-s $foo_path, "$foo_filename has size");

    ok(!defined($INC{$foo_filename}), "$foo has not been loaded yet");
    eval 'use Genome::Site';
    my $eval_error = $@;
    ok(!defined($INC{$foo_filename}), "$foo should fail to load");
    like($eval_error, qr/$foo_filename/, "eval error should mention site module: $foo_filename");
};

subtest 'catch missing dependencies' => sub {
    plan tests => 4;

    my $tmp_dir_path = tempdir();

    local %INC = %INC;
    local @INC = @INC;
    unshift @INC, $tmp_dir_path;

    my @site_dirs = qw(edu wustl genome);
    no warnings 'redefine';
    local *Genome::Site::site_dirs = sub { return @site_dirs };
    use warnings 'redefine';

    my $dep = 'Chautauquan';

    my $foo = Genome::Site::site_pkg($site_dirs[0]);
    my $foo_filename = Genome::Site::module_to_filename($foo);
    my $foo_path = create_module($tmp_dir_path, $foo,
        'package ' . $foo . ';',
        'use ' . $dep . ';',
        '1;',
    );
    ok(-s $foo_path, "$foo_filename has size");

    ok(!defined($INC{$foo_filename}), "$foo has not been loaded yet");
    eval 'use Genome::Site';
    my $eval_error = $@;
    ok(!defined($INC{$foo_filename}), "$foo should fail to load");
    like($eval_error, qr/$dep/, "eval error should mention missing dependency: $dep");
};

sub tempdir {
    # File::Temp->newdir() would DESTROY at eval 'use ...' so I switch to
    # tempdir and extracted method for documentation purposes.
    return File::Temp::tempdir();
}

sub create_module {
    my $dir_path = shift;
    my $module = shift;
    my @content = @_;

    unless ($dir_path && $module && @content) {
        croak 'invalid inputs';
    }

    my $filename = File::Spec->join(
        $dir_path,
        Genome::Site::module_to_filename($module),
    );

    my $dirname = (File::Spec->splitpath($filename))[1];
    make_path($dirname);

    my $file = IO::File->new($filename, 'w');
    $file->print(join("\n", @content), "\n");
    $file->close();

    return $filename;
}
