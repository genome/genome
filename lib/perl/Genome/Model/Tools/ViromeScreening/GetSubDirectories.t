#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More tests => 4;

print "user:  " . Genome::Sys->username . "\n";

BEGIN {use_ok('Genome::Model::Tools::ViromeScreening::GetSubDirectories');}

my ($dir) = ('/gscmnt/sata835/info/medseq/virome/test17');

#create
my $gsd = Genome::Model::Tools::ViromeScreening::GetSubDirectories->create(
                dir => $dir,
);
isa_ok($gsd, 'Genome::Model::Tools::ViromeScreening::GetSubDirectories');

ok($gsd->execute, 'getting sub directories');

my $sub_dirs = $gsd->sub_directories;
my $exists = 1;

foreach my $sub_dir (@$sub_dirs)
{
   $exists = 0 unless(-e $sub_dir); 
}

ok($exists, "subdirectories successfully retrieved");
