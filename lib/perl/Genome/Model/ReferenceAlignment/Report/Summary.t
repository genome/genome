#!/usr/bin/env genome-perl

use strict;
use warnings;

use File::Temp;


BEGIN {
    #use locks local to this test--must be set before Genome is loaded
    my $lock_tempdir = File::Temp::tempdir(CLEANUP => 1);
    $ENV{XGENOME_SITE_LOCK_DIR} = $lock_tempdir;

    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use above "Genome";
use Test::More skip_all => 'Test relies on live data that is no longer available.';
#use Test::More tests => 8;
use Genome::Sys;

my $tmp = Genome::Sys->create_temp_directory();

sub rmtree {
    system "/bin/rm -r $tmp";
    if (-e $tmp) {
        die "failed to remove directory $tmp: $!";
    }
}

die "temp in odd location! $tmp" unless $tmp =~ /\/tmp/;

if (-e $tmp) {
    warn "removing previously left-behind directory $tmp...";
    rmtree();
}

mkdir($tmp) or die "Failed to create directory $tmp: $!";

my $model = Genome::Model->get('id' => '2862809551');
unless ($model) {
    die "Can't find a model to work with";
}
my $build = $model->last_complete_build;
my $build_id = $build->id;
ok($build, "build found with id $build_id");

my $r = Genome::Model::ReferenceAlignment::Report::Summary->create(
    build_id => $build_id,
);
ok($r, "created a new report");

my @t = $r->report_templates;
is(scalar(@t),2, "got 2 templates") or diag(@t);

my $v = $r->generate_report;
ok($v, "generation worked");

my $result = $v->save($tmp);
ok($result, "saved to $tmp");

my $name = $r->name;
$name =~ s/ /_/g;

ok(-d "$tmp/$name", "report directory $tmp/$name is present");
ok(-e "$tmp/$name/report.txt", 'text report is present');
ok(-e "$tmp/$name/report.html", 'html report is present');

rmtree();

1;
