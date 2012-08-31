#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;
use Cwd;

use Test::More tests => 5;

BEGIN {
        use_ok('Genome::Model::Tools::Hgmi::Merge');
        use_ok('Genome::Model::Tools::Hgmi::DirBuilder');
}
#my $testpath = '/tmp/disk/analysis/HGMI/B_catenulatum/Bifidobacterium_catenulatum_BIFCATDFT_1.0_newb/Version_1.0/BAP/Version_1.0';
#my $tmpdir = tempdir("HGMI_XXXXXX", DIR => '/tmp/disk/analysis', CLEANUP => 1);
my $tmpdir = Genome::Sys->create_temp_directory();

my $d = Genome::Model::Tools::Hgmi::DirBuilder->create(
                    #'path' => "/tmp/disk/analysis/HGMI",
                    'path' => $tmpdir,
                    'org_dirname' => "B_catenulatum",
                    'assembly_version_name' => "Bifidobacterium_catenulatum_BIFCATDFT_1.0_newb",
                    'assembly_version' => "Version_1.0",
                    'pipe_version' => "Version_1.0",
                    'cell_type' => "BACTERIA");
isa_ok($d,'Genome::Model::Tools::Hgmi::DirBuilder');
ok($d->execute());

my $testpath = $tmpdir ."/B_catenulatum/Bifidobacterium_catenulatum_BIFCATDFT_1.0_newb/Version_1.0/BAP/Version_1.0";
my $orig_dir = getcwd;
diag($orig_dir);
&chdir($testpath);

my $m = Genome::Model::Tools::Hgmi::Merge->create(
  'organism_name' => "Bifidobacterium_catenulatum",
  'locus_tag' => "BIFCATDFT",
  'project_type' => "HGMI",
  'work_directory' => $testpath,
  'dev' => 1
 );

isa_ok($m,'Genome::Model::Tools::Hgmi::Merge');
&chdir($orig_dir);
remove \1, qw{ $tmpdir } ;

#my $cmd = $m->gather_details();
#diag( join(" ", @{$cmd[0]})."\n");

exit 0;
