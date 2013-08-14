#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Genome::Utility::Test;
use Test::More skip_all => 'bypass until expected data is updated';

my $id = 2891454740;
my $model = Genome::Model->get($id);
ok($model, "got test model");

my $outfile = Genome::Sys->create_temp_file_path(); 
my $result = Genome::Model::Command::Export::Metadata->execute(models => [$model], output_path => $outfile);
ok($result, "ran");
ok(-e $outfile, "outfile $outfile exists");

my $expected = __FILE__ . '.expected-output';
Genome::Utility::Test::compare_ok(
    $outfile, 
    $expected,
    'output matches',
    filters => [
        sub{ 
            my $line = shift;
            return if $line =~ /total_kb/i; # ignore disk volumes
            return $line;
        },
    ],
);

done_testing(4);
