#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Genome::InstrumentData::InstrumentDataTestObjGenerator;
use Test::More;

my $bam_path = $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-AlignmentResult-Bwa/input.bam';

use_ok('Genome::InstrumentData::AlignmentResult');
# Test Intermediate AR Class
class Genome::InstrumentData::IntermediateAlignmentResult::Tester {
    is=>['Genome::InstrumentData::IntermediateAlignmentResult'],
};
# Test AR Class
class Genome::InstrumentData::AlignmentResult::Tester {
    is => 'Genome::InstrumentData::AlignmentResult',
};
my $iar_id;
sub Genome::InstrumentData::AlignmentResult::Tester::_run_aligner { 
    my $self = shift;

    # Create IAR to test if they are deleted
    my $iar = Genome::InstrumentData::IntermediateAlignmentResult::Tester->__define__();
    ok($iar, 'create tester iar');
    $iar_id = $iar->id;
    $iar->add_user(user => $self, label => 'intermediate result');
    my @users = $iar->users;
    ok(@users, 'alignment result uses the iar');

    # Copy bam
    my $copy_ok = Genome::Sys->copy_file($bam_path, $self->temp_staging_directory.'/all_sequences.bam');
    ok($copy_ok, 'copied bam');

    return 1;
};
sub Genome::InstrumentData::AlignmentResult::Tester::aligner_params_for_sam_header { 'align me bro!' };
sub Genome::InstrumentData::AlignmentResult::Tester::estimated_kb_usage { 0 };
sub Genome::InstrumentData::AlignmentResult::Tester::fillmd_for_sam { 0 };
sub Genome::InstrumentData::AlignmentResult::Tester::requires_fastqs_to_align { 0 };

my $inst_data = Genome::InstrumentData::InstrumentDataTestObjGenerator::create_solexa_instrument_data($bam_path);
ok($inst_data, 'create inst data');
my $reference_model = Genome::Model::ImportedReferenceSequence->get(name => 'TEST-human');
ok($reference_model, "got reference model");
my $reference_build = $reference_model->build_by_version('1');
ok($reference_build, "got reference build");

my $alignment_result = Genome::InstrumentData::AlignmentResult::Tester->create(
    id => -1337,
    instrument_data_id => $inst_data->id,
    reference_build => $reference_build,
    aligner_name => 'tester',
    aligner_version => '1',
    aligner_params => '',
);
ok($alignment_result, 'defined alignment result');
isa_ok($alignment_result, 'Genome::InstrumentData::AlignmentResult::Tester');

# Check that the iar is deleted
my @users = Genome::SoftwareResult::User->get(user => $alignment_result, label => 'intermediate result');
ok(!@users, 'alignment result is not using any intermediate results');
ok(!Genome::InstrumentData::IntermediateAlignmentResult::Tester->get($iar_id), 'did not get the intemediate result');

#  define qc result
my $qc_result = Genome::InstrumentData::AlignmentResult::Merged::BamQc->__define__(alignment_result_id => $alignment_result->id);
ok($qc_result, 'define qc result for alienment result');
ok($alignment_result->add_user(user => $qc_result, label => 'uses'), 'add qc result as user of alignment result');

# delete
ok($alignment_result->delete, 'delete');
ok(ref($qc_result) eq 'UR::DeletedRef', 'deleted qc result');

done_testing();
