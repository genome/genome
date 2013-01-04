#!/usr/bin/env genome-perl
use strict;
use warnings;
use above 'Genome';
use Test::More tests => 3;

my $s = Genome::Sample->create(id => -1, name => 'TEST-patient1-sampleX');
ok($s, 'created a test sample');

class Genome::Model::Foo { is => 'Genome::Model' };
sub Genome::Model::Foo::publication_description { "My methods!" }

my $p = Genome::ProcessingProfile::Foo->create(id => -1, name => "test profile X");

my $m1 = Genome::Model::Foo->create(subject => $s, processing_profile => $p);
ok($m1, "created a model");

my $out;
open TESTOUT, '>', \$out or die "failed to open output fh to variable";
select TESTOUT;
Genome::Model::Command::Export::PublicationMethods->execute(models => [$m1]);
select STDOUT;

my $expected_out = "\nMy methods!\n\n";
is($out, $expected_out, "got expected outputs");
