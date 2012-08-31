#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 6;

BEGIN {
        use_ok('Genome::Model::Tools::BioDbFasta::Build');
}

my $seq_dir1 = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-BioDbFasta/Build1";
my $seq_dir2 = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-BioDbFasta/Build-rebuild";

my $b = Genome::Model::Tools::BioDbFasta::Build->create( 'dir' => $seq_dir1 );
isa_ok($b,'Genome::Model::Tools::BioDbFasta::Build');
ok($b);


my $b1 = Genome::Model::Tools::BioDbFasta::Build->create( 'dir' => "/This/Directory/Dont/Exist" );

is($b1->execute(),0,'directory non-existent');

ok($b->execute(),'build index');

my $b_rebuild = Genome::Model::Tools::BioDbFasta::Build->create( 'dir' => $seq_dir2,
                                                                 'rebuild' => 1);

ok($b_rebuild->execute(),'rebuild index');
