#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use File::Remove qw/ remove /;
use File::Temp qw/ tempdir/;
use English;

use Test::More tests => 3;

BEGIN {
        use_ok('Genome::Model::Tools::Hgmi::DirBuilder');
}

my $dir = tempdir("HGMI_XXXXXX", CLEANUP => 1, DIR => File::Spec->tmpdir() );

my $tool_db = Genome::Model::Tools::Hgmi::DirBuilder->create(
                    path => $dir,
                    'org_dirname' => "B_catenulatum",
                    'assembly_version_name' => "Bifidobacterium_catenulatum_BIFCATDFT_1.0_newb",
                    'assembly_version' => "Version_1.0",
                    'pipe_version' => "Version_1.0",
                    'cell_type' => "BACTERIA");
isa_ok($tool_db,'Genome::Model::Tools::Hgmi::DirBuilder');
ok($tool_db->execute,'execute dir builder');

