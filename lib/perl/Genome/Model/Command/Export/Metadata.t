#!/usr/bin/env genome-perl
use strict;
use warnings;
use above "Genome";
use Genome::Utility::Test;
use Test::More;

my $id = 2891454740;
my $model = Genome::Model->get($id);
ok($model, "got test model");

my $expected = __FILE__ . '.expected-output';

my $intermediate_outfile = Genome::Sys->create_temp_file_path();

my $scrubbed_outfile;
if ($ARGV[0] && $ARGV[0] eq 'REBUILD') {
    print "\n\nRebuilding test result\n";
    $scrubbed_outfile = $expected;    
    unlink $expected;
}
else {
    $scrubbed_outfile = Genome::Sys->create_temp_file_path(); 
}

my $result = Genome::Model::Command::Export::Metadata->execute(models => [$model], output_path => $intermediate_outfile, verbose => 1);
ok($result, "ran");
ok(-e $intermediate_outfile, "intermediate_outfile $intermediate_outfile exists");

Genome::Sys->shellcmd(cmd => "grep -v Genome::Disk::Allocation <$intermediate_outfile | grep -v Genome::Disk::Volume >  $scrubbed_outfile");

Genome::Utility::Test::compare_ok(
    $scrubbed_outfile, 
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
