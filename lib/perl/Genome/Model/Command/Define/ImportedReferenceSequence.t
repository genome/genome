#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN { delete $ENV{UR_DBI_NO_COMMIT} }

use above 'Genome';
use Test::More;
use Genome::Utility::Test qw(assert_using_test_db);

assert_using_test_db();

Genome::Report::Email->silent();

if (Genome::Sys->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

my $version = 1;
my $cmd_class = 'Genome::Model::Command::Define::ImportedReferenceSequence';
use_ok($cmd_class);

my $data_dir = __FILE__.'.d';
my $pp = Genome::ProcessingProfile::ImportedReferenceSequence->create(name => 'test_ref_pp');
my $taxon = Genome::Taxon->get(name => 'human');
my $patient = Genome::Individual->create(name => "test-patient", common_name => 'testpat', taxon => $taxon);
my $sample = Genome::Sample->create(name => "test-patient", common_name => 'tumor', source => $patient);
ok($sample, 'created sample');

my $sequence_uri = "http://genome.wustl.edu/foo/bar/test.fa.gz";

my $fasta_file1 = "$data_dir/data.fa";
my $fasta_file2 = "$data_dir/data2.fa";

my @params = (
    "--fasta-file=$fasta_file1",
    "--model-name=test-ref-seq-1",
    "--processing-profile=".$pp->id,
    "--species-name=human",
    "--subject=".$sample->id,
    "--version=42",
    "--sequence-uri=".$sequence_uri,
    );

print Data::Dumper::Dumper(\@params) . "\n";

my $rv = $cmd_class->_execute_with_shell_params_and_return_exit_code(@params);
is($rv, 0, 'executed command');
my $model = Genome::Model::ImportedReferenceSequence->get(name => 'test-ref-seq-1');
ok($model, 'Found newly created model');
my $build = $model->last_complete_build;
ok($build, 'Found a completed build');
is($build->version, 42, 'Build has correct version');
is($build->sequence_uri, $sequence_uri, "sequence uri matches");

# specify derived_from
@params = (
    "--derived-from=".$build->name,
    "--fasta-file=$fasta_file1",
    "--model-name=test-ref-seq-2",
    "--processing-profile=".$pp->id,
    "--species-name=human",
    "--subject=".$sample->id,
    "--version=26",
    "--sequence-uri=".$sequence_uri,
    );
$rv = $cmd_class->_execute_with_shell_params_and_return_exit_code(@params);
is($rv, 0, 'executed command');
my $coords_model = Genome::Model::ImportedReferenceSequence->get(name => 'test-ref-seq-2');
ok($coords_model, 'Found newly created model');
my $d1_build = $coords_model->last_complete_build;
ok($d1_build, 'Found a completed build');
is($d1_build->version, 26, 'Build has correct version');
is($d1_build->derived_from->id, $build->id, 'derived_from property is correct');
is($d1_build->coordinates_from->id, $build->id, 'coordinates_from property is correct');
ok($d1_build->is_compatible_with($build), 'coordinates_from build is_compatible_with parent build');
is($d1_build->sequence_uri, $sequence_uri, "sequence uri matches");
ok($build->is_compatible_with($d1_build), 'parent build is_compatible_with coordinates_from build');

# derive from d1_build
@params = (
    "--derived-from=".$d1_build->id,
    "--fasta-file=$fasta_file1",
    "--model-name=test-ref-seq-3",
    "--processing-profile=".$pp->id,
    "--species-name=human",
    "--subject=".$sample->id,
    "--version=96",
    "--sequence-uri=".$sequence_uri,
    );
$rv = $cmd_class->_execute_with_shell_params_and_return_exit_code(@params);
is($rv, 0, 'executed command');
my $derived_model = Genome::Model::ImportedReferenceSequence->get(name => 'test-ref-seq-3');
ok($derived_model, 'Found newly created model');
my $d2_build = $derived_model->last_complete_build;
ok($d2_build, 'Found a completed build');
is($d2_build->version, 96, 'Build has correct version');
is($d2_build->derived_from->id, $d1_build->id, 'derived_from property is correct');
is($d2_build->coordinates_from->id, $build->id, 'coordinates_from property is correct');
ok($d2_build->is_compatible_with($d1_build), 'derived build is_compatible_with parent build');
ok($d1_build->is_compatible_with($d2_build), 'derived build is_compatible_with parent build');
ok($d2_build->is_compatible_with($build), 'derived build is_compatible_with parent build');
is($d2_build->sequence_uri, $sequence_uri, "sequence uri matches");
ok($build->is_compatible_with($d2_build), 'parent build is_compatible_with derived build');

@params = (
    "--append-to=".$d1_build->id,
    "--fasta-file=$fasta_file2",
    "--model-name=test-ref-seq-4",
    "--processing-profile=".$pp->id,
    "--species-name=human",
    "--subject=".$sample->id,
    "--version=append",
    "--sequence-uri=".$sequence_uri,
    );
$rv = $cmd_class->_execute_with_shell_params_and_return_exit_code(@params);
is($rv, 0, 'executed command');
my $append_model = Genome::Model::ImportedReferenceSequence->get(name => 'test-ref-seq-4');
ok($append_model, 'Found newly created model');
my $a_build = $append_model->last_complete_build;
ok($a_build, 'Found a completed build');
is($a_build->version, 'append', 'Build has correct version');
is($a_build->derived_from->id, $d1_build->id, 'derived_from property is correct');
is($a_build->append_to->id, $d1_build->id, 'append_to property is correct');
is($a_build->coordinates_from->id, $build->id, 'coordinates_from property is correct');
ok($a_build->is_compatible_with($d1_build), 'derived build is_compatible_with parent build');
ok($d1_build->is_compatible_with($a_build), 'derived build is_compatible_with parent build');
ok($a_build->is_compatible_with($build), 'derived build is_compatible_with parent build');
is($a_build->sequence_uri, $sequence_uri, "sequence uri matches");
ok($build->is_compatible_with($a_build), 'parent build is_compatible_with derived build');
ok($a_build->primary_consensus_path('fa') ne $a_build->full_consensus_path('fa'), 'primary sequence != full sequence in appended build');
my $primary_sz = -s $a_build->primary_consensus_path('fa');
my $full_sz = -s $a_build->full_consensus_path('fa');
is($primary_sz, -s $fasta_file2, 'primary sequence file has correct size');
is($full_sz, (-s $fasta_file1) + (-s $fasta_file2), 'full sequence is correct size');


done_testing();
