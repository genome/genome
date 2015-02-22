#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above "Genome";

use Test::More tests => 52;

use_ok('Genome::FeatureList::Command::Merge');

my $test_bed_dir = join('.', __FILE__, 'd');

my $test_refseq = Genome::Model::Build::ReferenceSequence->get(name => 'NCBI-human-build36');

my @input_feature_lists = (undef); #fill 0th position
for (1..6) {
    my $file = join("/", $test_bed_dir,  "${_}.bed");
    my $next_feature_list = Genome::FeatureList->create(
        file_path => $file,
        file_content_hash => Genome::Sys->md5sum($file),
        format => ($_ == 1? 'true-BED' : 'multi-tracked'),
        reference => $test_refseq,
        content_type => 'targeted',
        source => 'G:FL:C:Merge test',
        name => "G:FL:C:Merge test FeatureList #$_",
    );
    isa_ok($next_feature_list, 'Genome::FeatureList', "created test feature list $_");
    push @input_feature_lists, $next_feature_list;
}

#Should Pass:
#2,3 should combine both tracks
#2,4 should combine only probes
#2,5 should combine only targets

#Should Fail:
#1,2 has mismatch in track count
#2,6 has unknown track name
#4,5 has no non-empty tracks in common
my @should_pass = ([2,3, ['probes','targets']], [2,4,['probes']], [2,5,['targets']]);
my @should_fail = ([1,2,qr(all be single tracked or all be multi-tracked)], [2,6,qr(mysterious track name)], [4,5,qr(Could not find any complete track)]);

for my $test (@should_pass) {
    my $name = 'test merging of ' . join(' and ', @$test[0,1]);
    my $merge_cmd = Genome::FeatureList::Command::Merge->create(
        source_lists => [@input_feature_lists[@$test[0,1]]],
        name => $name,
    );
    isa_ok($merge_cmd, 'Genome::FeatureList::Command::Merge', "created command to $name");
    ok($merge_cmd->execute, "executed command to $name");

    my $new_feature_list = Genome::FeatureList->get(name => $name);
    ok($new_feature_list, 'created new feature list');
    is($new_feature_list->format, 'multi-tracked', 'format set properly on merged list');
    is($new_feature_list->reference, $test_refseq, 'reference set properly on merged list');

    my $merged_bed_content = Genome::Sys->read_file($new_feature_list->file_path);
    ok($merged_bed_content, "new BED file was produced for $name");
    my $hashed_content = Genome::FeatureList::Command::Merge->hash_bed_content_by_tracks($merged_bed_content);
    my @tracks = keys %$hashed_content;
    is_deeply([sort @tracks], $test->[2], "expected tracks were added to BED file for $name");

    my $merge_cmd_again = Genome::FeatureList::Command::Merge->create(
        source_lists => [@input_feature_lists[@$test[0,1]]],
        name => $name . ' again',
    );
    isa_ok($merge_cmd_again, 'Genome::FeatureList::Command::Merge', "created second command to $name");
    ok($merge_cmd_again->execute, "executed second command to $name");

    my $repeated_feature_list = Genome::FeatureList->get(name => $name);
    is($repeated_feature_list, $new_feature_list, 'returned same FeatureList when merging same two inputs again');
}

for my $test (@should_fail) {
    my $name = 'test merging of ' . join(' and ', @$test[0,1]);

    my $merge_cmd = Genome::FeatureList::Command::Merge->create(
        source_lists => [@input_feature_lists[@$test[0,1]]],
        name => $name,
    );
    isa_ok($merge_cmd, 'Genome::FeatureList::Command::Merge', "created command to $name");

    my $result = eval { $merge_cmd->execute; };
    my $error = $@;
    ok(!$result, 'command failed to execute');
    ok($error, 'execute crashed due to errors');
    like($error, $test->[2], 'crashed with expected error');
    my $new_feature_list = Genome::FeatureList->get(name => $name);
    ok(!$new_feature_list, 'no new feature-list created');
}
