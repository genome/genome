#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above "Genome";
use Genome::Test::Factory::Model::ReferenceSequence;

use Test::More tests => 12;
use Test::Exception;

#test refseq checking
my ($incompatible_reference, $combined_reference, @compatible_references) = map { Genome::Test::Factory::Model::ReferenceSequence->setup_reference_sequence_build } (1..5);

write_seqdict_file($compatible_references[0], 1,2,3);
write_seqdict_file($compatible_references[1], 1,2);
write_seqdict_file($compatible_references[2], 1);
write_seqdict_file($combined_reference, 1,2,3,4);
write_seqdict_file($incompatible_reference, 4);

$compatible_references[0]->derived_from($compatible_references[1]);
$compatible_references[0]->coordinates_from($compatible_references[2]);
$compatible_references[1]->derived_from($compatible_references[2]);
$compatible_references[1]->coordinates_from($compatible_references[2]);
ok($compatible_references[1]->is_compatible_with($compatible_references[0]), 'test references are compatible');
ok($compatible_references[2]->is_compatible_with($compatible_references[0]), 'test references are compatible');

my $common_reference;
lives_and(sub {
    $common_reference = Genome::FeatureList::Command::Merge->find_common_reference(@compatible_references);
    is($common_reference, $compatible_references[0]);
}, 'found the common reference between three compatible references');

dies_ok(sub {
    Genome::FeatureList::Command::Merge->find_common_reference(@compatible_references, $incompatible_reference);
}, 'cannot find common reference between incompatible references with no converters');

lives_and(sub {
    $common_reference = Genome::FeatureList::Command::Merge->find_common_reference(($incompatible_reference) x 5);
    is($common_reference, $incompatible_reference);
}, 'returns the same reference when all have the same reference');

my @refs_to_combine = ($compatible_references[0], $compatible_references[1], $incompatible_reference);
map { $combined_reference->add_input(
    name => 'combines',
    value_id => $_->id,
    value_class_name => $_->class,
) } @refs_to_combine;
lives_and(sub {
    $common_reference = Genome::FeatureList::Command::Merge->find_common_reference(@refs_to_combine);
    is($common_reference, $combined_reference);
}, 'selected combined reference');

dies_ok(sub {
    $common_reference = Genome::FeatureList::Command::Merge->find_common_reference(@refs_to_combine, $compatible_references[2]);
}, 'does not find combined reference when one reference was not a part of it');

my $converter = Genome::Model::Build::ReferenceSequence::Converter->create(
    source_reference_build => $compatible_references[2],
    destination_reference_build => $incompatible_reference,
    algorithm => 'prepend_chr',
);
isa_ok($converter, 'Genome::Model::Build::ReferenceSequence::Converter', 'defined a converter between references');
lives_and(sub {
    $common_reference = Genome::FeatureList::Command::Merge->find_common_reference($compatible_references[2], $incompatible_reference);
    is($common_reference, $converter->destination_reference_build);
}, 'found convertible reference');

my ($source_reference, $intermediate_reference, $other_reference, $combined_reference) = map { Genome::Test::Factory::Model::ReferenceSequence->setup_reference_sequence_build } (1..4);
write_seqdict_file($combined_reference, 1,2);
write_seqdict_file($other_reference, 2);
write_seqdict_file($source_reference, 'chr1');
write_seqdict_file($intermediate_reference, 1);

for my $ref ($intermediate_reference, $other_reference) {
    $combined_reference->add_input(name => 'combines', value => $ref);
}

dies_ok(sub {
    my $found_reference = Genome::FeatureList::Command::Merge->find_common_reference($source_reference, $combined_reference);
}, 'no common reference known');

my $intermediate_converter = Genome::Model::Build::ReferenceSequence::Converter->create(
    source_reference_build => $source_reference,
    destination_reference_build => $intermediate_reference,
    algorithm => 'chop_chr',
);
isa_ok($converter, 'Genome::Model::Build::ReferenceSequence::Converter', 'defined another converter between references');

lives_and(sub {
    my $found_reference = Genome::FeatureList::Command::Merge->find_common_reference($source_reference, $combined_reference);
    is($found_reference, $combined_reference, 'found combination through conversion');
}, 'found combined reference');

sub write_seqdict_file {
    my $reference = shift;
    my @chrs_to_include = @_;

    my $seqdict_dir = File::Spec->join($reference->data_directory, 'seqdict');
    Genome::Sys->create_directory($seqdict_dir);
    my $file = File::Spec->join($seqdict_dir, 'seqdict.sam');
    warn $file;
    Genome::Sys->write_file(
        $file,
        map( sprintf('@SQ	SN:%s	LN:5	UR:file:///dev/null	AS:test	SP:Unknown' . "\n", $_), @chrs_to_include),
    );
}
