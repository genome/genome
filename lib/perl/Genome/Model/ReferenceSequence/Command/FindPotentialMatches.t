#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use Test::More tests => 3;
use Test::Deep qw(cmp_bag);

use above 'Genome';

use Genome::Test::Factory::Build;
use Genome::Test::Factory::Model::ReferenceSequence;

my $pkg = 'Genome::Model::ReferenceSequence::Command::FindPotentialMatches';
use_ok($pkg);

my $taxon = Genome::Taxon->create(
    name => 'test taxon for find-potential-matches',
);

my @refseqs;

for (1..5) {
    my $model = Genome::Test::Factory::Model::ReferenceSequence->setup_object(
        subject => $taxon,
    );
    my $data_dir = File::Spec->join(__FILE__.'.d', $_);
    my $build = Genome::Test::Factory::Build->setup_object(
        model_id => $model->id,
        data_directory => $data_dir,
        status => 'Succeeded',
    );

    push @refseqs, $build;
}

subtest 'gather query chromosomes from file' => sub {
    my $file_path = Genome::Sys->create_temp_file_path;
    my @chromosome_list = (1..5, 'X', 'something|else');
    Genome::Sys->write_file($file_path, map { "$_\tignored\tcolumns\there\n" } @chromosome_list);

    my $cmd = $pkg->create(
        query_file => $file_path,
    );
    isa_ok($cmd, $pkg);

    my $result = $cmd->_query_chromosomes;
    cmp_bag($result, \@chromosome_list, 'found expected chromosomes');
};

subtest 'command produces results' => sub {

    my $file_path = Genome::Sys->create_temp_file_path;
    my @chromosome_list = (1..4);
    Genome::Sys->write_file($file_path, map { "$_\n" } @chromosome_list);

    subtest 'command with default matching options' => sub {
        my $cmd = $pkg->create(
            taxon => $taxon,
            query_file => $file_path,
        );
        isa_ok($cmd, $pkg);

        ok($cmd->execute, 'executed command');

        my $exact_matches = $cmd->exact_matches;
        is(scalar(@$exact_matches), 1, 'found an exact match');
        my $supersets = $cmd->supersets;
        is(scalar(@$supersets), 1, 'found a superset');
        my $subsets = $cmd->subsets;
        is($subsets, undef, 'subsets not set when not requested');
        my $partial_matches = $cmd->partial_matches;
        is($partial_matches, undef, 'partial matches not set when not requested');
    };

    subtest 'command with all matching options' => sub {
        my $cmd = $pkg->create(
            taxon => $taxon,
            query_file => $file_path,
            include_subsets => 1,
            include_supersets => 1,
            include_partial_matches => 1,
        );
        isa_ok($cmd, $pkg);

        ok($cmd->execute, 'executed command');

        my $exact_matches = $cmd->exact_matches;
        is(scalar(@$exact_matches), 1, 'found an exact match');
        my $supersets = $cmd->supersets;
        is(scalar(@$supersets), 1, 'found a superset');
        my $subsets = $cmd->subsets;
        is(scalar(@$subsets), 1, 'found a subset');
        my $partial_matches = $cmd->partial_matches;
        is(scalar(@$partial_matches), 1, 'found a partial match');
    };
};
