#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above "Genome";

use Test::More tests => 12;

use_ok('Genome::FeatureList::Command::Import');

my $reference = Genome::Model::Build::ReferenceSequence->get(name => 'NCBI-human-build36');

my @test_beds = (undef); #fill 0th position
push @test_beds, map { sprintf('%s.d/%s.bed', __FILE__, $_) } 1..4;


#1 has unknown track name
#3 has mismatched chromosomes
my @should_fail = (1,3);

#2 is normal
#4 has browser lines and alternate track names
my @should_pass = (2,4);

my @new_feature_lists;
for my $passing_bed (@test_beds[@should_pass]) {
    my $cmd = Genome::FeatureList::Command::Import->create(
        file_path => $passing_bed,
        name => 'test ' . substr($passing_bed, -5),
        reference => $reference,
    );
    isa_ok($cmd, 'Genome::FeatureList::Command::Import', 'created command');
    ok($cmd->execute, 'executed command');
    my $new_feature_list = $cmd->new_feature_list;
    isa_ok($new_feature_list, 'Genome::FeatureList', 'created new feature-list');
    push @new_feature_lists, $new_feature_list;
}
is($new_feature_lists[0]->file_content_hash, $new_feature_lists[1]->file_content_hash, 'sanitization of browser lines led to same BED file');

for my $failing_bed(@test_beds[@should_fail]) {
    my $cmd = Genome::FeatureList::Command::Import->create(
        file_path => $failing_bed,
        name => 'test ' . substr($failing_bed, -5),
        reference => $reference,
    );
    isa_ok($cmd, 'Genome::FeatureList::Command::Import', 'created command');
    my $rv = eval { $cmd->execute; };
    my $error = $@;
    ok(!$rv && $error, 'failed to execute');
}
