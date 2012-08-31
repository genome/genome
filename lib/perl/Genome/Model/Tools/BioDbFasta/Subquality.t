#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";
use Test::More tests => 5;

BEGIN {
        use_ok('Genome::Model::Tools::BioDbFasta::Subquality');
}

my $start = 1;
my $stop = 5;

my $dir_noexist = "/idontexist";
my $dir_noindex = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-BioDbFasta/Subquality/not-indexed";

my $seq_name0 = "noname";

my $s_noex = Genome::Model::Tools::BioDbFasta::Subquality->create(
                                                                  dir => $dir_noexist,
                                                                  name => $seq_name0,
                                                                  start => $start,
                                                                  stop => $stop
                                                                  );

is($s_noex->execute(),0,'dir doesn\'t exist');

my $s_index = Genome::Model::Tools::BioDbFasta::Subquality->create(
                                                                  dir => $dir_noindex,
                                                                  name => $seq_name0,
                                                                  start => $start,
                                                                  stop => $stop
                                                                  );

is($s_index->execute(),0,'no index in directory');

my $bad_index =  $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-BioDbFasta/Subquality/bad-index";
my $s_bi = Genome::Model::Tools::BioDbFasta::Subquality->create(
                                                             dir => $bad_index,
                                                             name => "seq1",
                                                             start => $start,
                                                             stop => $stop
                                                             );
my $ret ;
eval {$ret = $s_bi->execute();};
ok(!$ret,'should confess on the bad directory index');

my $qual_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-Tools-BioDbFasta/Subquality/sequence";


my $s = Genome::Model::Tools::BioDbFasta::Subquality->create(
                                                             dir => $qual_dir,
                                                             name => "seq1",
                                                             start => $start,
                                                             stop => $stop
                                                             );


ok($s->execute());

