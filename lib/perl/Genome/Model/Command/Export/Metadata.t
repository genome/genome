#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Test::More tests => 4;

my $id = 2891454740;
my $model = Genome::Model->get($id);
ok($model, "got test model");

my $outfile = Genome::Sys->create_temp_file_path(); 
my $result = Genome::Model::Command::Export::Metadata->execute(models => [$model], output_path => $outfile);
ok($result, "ran");
ok(-e $outfile, "outfile $outfile exists");

my $expected = __FILE__ . '.expected-output';
my @diff = grep { $_ !~ /total_kb/ } `sdiff -s -w 500 $expected $outfile`;
is(scalar(@diff), 0, "no differences")
    or diag(@diff);



