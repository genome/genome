use strict;
use warnings;

use above "Genome";
use Genome::Utility::Test qw(run_ok);
use Test::Builder;
use Test::More;

use File::Spec qw();
use IPC::Cmd qw(can_run);
use IPC::System::Simple qw(capture);

my $start_point = resolve_start_point();

my @files_to_compile = files_to_compile($start_point);;
plan tests => scalar(@files_to_compile) + 1;
ok(1, 'need at least one test to pass'); # need a test for when @files_to_compile == 0

for my $file (@files_to_compile) {
    my $pid = fork();
    if (!defined $pid) {
        die 'failed to fork';
    } elsif ($pid) {
        waitpid $pid, 0;
        my $rel_path = File::Spec->abs2rel($file);
        ok($? == 0, qq(compiled '$rel_path'));
    } else {
        my $exit = compile_file($file);
        exit($exit);
    }
}


sub compile_file {
    my $file = shift;
    my @output = qx(perl -c "$file" 2>&1);
    my $exit = $? >> 8;
    if ($exit != 0) {
        diag @output;
    }
    return $exit;
}


sub is_perl_file {
    my $file = shift;
    return ($file =~ /\.(pm|t|pl)$/ && -f $file);
}


sub resolve_start_point {
    my $start_point_arg = (shift @ARGV) || '';

    if ($start_point_arg) {
        return $start_point_arg;
    }

    if ($ENV{JENKINS_STABLE_REVISION}) {
        return $ENV{JENKINS_STABLE_REVISION};
    }

    return '@{u}';
}


sub is_blacklisted {
    my $file = shift;
    my @blacklist = (
        qr(Db/Ensembl/Command/Vep.d/),
        qr(Site/CLIA.pm$),
        qr(Site/CLIA.t$),
    );
    return grep { $file =~ /$_/ } @blacklist;
}


sub files_to_compile {
    my $start_point = shift;
    my @files;

    my $git_dir = capture('git', 'rev-parse', '--show-toplevel');
    chomp $git_dir;

    my @cmd = ('git', 'diff', '--name-only', $start_point);
    @files = capture(@cmd);
    chomp @files;

    @files = map { File::Spec->join($git_dir, $_) } @files;
    @files = grep { is_perl_file $_ } @files;
    @files = grep { ! is_blacklisted $_ } @files;

    return @files;
}
