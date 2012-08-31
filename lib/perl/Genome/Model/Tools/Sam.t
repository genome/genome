#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 5;

use_ok('Genome::Model::Tools::Sam');

my $sam = Genome::Model::Tools::Sam->create();

isa_ok($sam, 'Genome::Model::Tools::Sam');

# should get default versions since we do not specify
my $sam_version    = $sam->use_version();  

ok(-e $sam->path_for_samtools_version($sam_version), "samtools version $sam_version exists");

my $time = Genome::Model::Tools::Sam->time();
my $time_format_is_good = ($time =~ /\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}Z/) ? 1 : 0;
ok($time_format_is_good, "time returns something like YYYY-MM-DDThh:mm:ssTZD");

my $date = Genome::Model::Tools::Sam->date();
my $date_format_is_good = ($date =~ /\d{4}-\d{2}-\d{2}/) ? 1 : 0;
ok($date_format_is_good, "date returns something like YYYY-MM-DD");
