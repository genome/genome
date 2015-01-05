#!/usr/bin/env genome-perl

use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_DUMP_STATUS_MESSAGES} = 1;
    $ENV{UR_DUMP_DEBUG_MESSAGES} = 1;
}

use above 'Genome';
use Test::More;
use Genome::InstrumentData::InstrumentDataTestObjGenerator;

my $TEST_BWA_VERSION = '0.7.10';
my $TEST_SAMTOOLS_VERSION = '0.1.19';
my $samtools_path = Genome::Model::Tools::Sam->path_for_samtools_version($TEST_SAMTOOLS_VERSION);

my $bam_path = $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-AlignmentResult-Bwa/input.bam';
my $inst_data = Genome::InstrumentData::InstrumentDataTestObjGenerator::create_solexa_instrument_data(
    $bam_path
    );

ok($inst_data, 'create inst data');
my $reference_model = Genome::Model::ImportedReferenceSequence->get(name => 'TEST-human');
ok($reference_model, "got reference model");
my $reference_build = $reference_model->build_by_version('1');
ok($reference_build, "got reference build");

my $alnidx = Genome::Model::Build::ReferenceSequence::AlignerIndex->get_with_lock(
    aligner_name => "bwa",
    reference_build_id => $reference_build->id,
    aligner_version => '0.7.10'
    );

ok($alnidx, "got aligner index");

subtest "bwamem" => make_test_suite("Genome::InstrumentData::AlignmentResult::Bwamem", "bwamem");
subtest "bwamem-stream" => make_test_suite("Genome::InstrumentData::AlignmentResult::BwamemStream", "bwamem-stream");

done_testing();

sub make_test_suite {
    my ($pkg, $aligner_name) = @_;

    use_ok($pkg);

    return sub {
        subtest "parse params" => sub {
            my @invalid_strings = (
                '-o -m -g',
                'banana',
                '-c should_be_int',
                '-t should_be_int',
                '-r should_be_float',
                '-c should_be_int',
                '-D should_be_float',
                '-W should_be int',
                '-m should_be_int',
                '-S should_be_null',
                '-P should_be_null',
                '-e should_be_null',
                );

            for my $invalid (@invalid_strings) {
                eval {
                    $pkg->_param_string_to_hash($invalid);
                };

                ok($@, "Param string $invalid is invalid");
            }

            my $params = $pkg->_param_string_to_hash("-t 4 -M");
            ok(exists $params->{M}, "-M flag exists");
            ok(exists $params->{t}, "-t flag exists");
            ok(!$params->{M}, "-M doesn't have a value");
            is($params->{t}, 4, "parsed 4 threads");
        };

        subtest "required rusage" => sub {
            for my $i (4..8) {
                my $rusage = $pkg->required_rusage(
                    aligner_params => "-t $i"
                    );
                ok($rusage =~ /cpus >= $i .* -n $i/, "threads set correctly (-t $i)");
            }
        };

        subtest "align" => sub {
            my $ar = make_alignment_result($aligner_name, "-t 2 -M");
            ok($ar, "Created alignment result");

            is($ar->_aligner_index_fasta, $alnidx->full_consensus_path("fa"),
                "aligner index fasta points to the expected location");

            my $dir = $ar->disk_allocations->absolute_path;
            ok(-d $dir, "Output directory exists");
            my @expected_files = qw(
                all_sequences.bam
                all_sequences.bam.bai
                all_sequences.bam.flagstat
                all_sequences.bam.md5
                );

            for my $f (@expected_files) {
                my $path = File::Spec->catfile($dir, $f);
                ok(-s $path, "File $f exists");
            }

            print qx{cat $dir/all_sequences.bam.md5};

            subtest "validate header" => sub {
                my $bam = File::Spec->catfile($dir, "all_sequences.bam");
                my @header = qx{$samtools_path view -H $bam};
                chomp @header;
                my @pg_line = grep {/^\@PG\t/} @header;
                is(scalar @pg_line, 1, "Found exactly one PG line");
                ok($pg_line[0] =~ /bwa mem/, "PG line contains 'bwa mem'");
                ok($pg_line[0] !~ /bwa mem.* -t/, "PG line does not contain '-t' parameter'");
                ok($pg_line[0] =~ /bwa mem.* -M/, "PG line does contain '-M' parameter'");

                my @rg_line = grep {/^\@RG\t/} @header;
                is(scalar @rg_line, 1, "Found exactly one RG line");
                my @rg_fields = split("\t", $rg_line[0]);
                shift @rg_fields;
                my %rg_data = map {split(":", $_, 2)} @rg_fields;
                my $expected_rg = $ar->read_and_platform_group_tag_id;
                my $expected_lb = $inst_data->library_name;
                my $expected_sm = $inst_data->sample_name;
                is($rg_data{ID}, $expected_rg, "Read group is correct");
                is($rg_data{LB}, $expected_lb, "Library is correct");
                is($rg_data{SM}, $expected_sm, "Sample is correct");
            };

        };
    };
}

my $next_id = -1;
sub make_alignment_result {
    my ($aligner_name, $bwa_params) = @_;

    my $alignment_result = Genome::InstrumentData::AlignmentResult->create(
        id => $next_id,
        instrument_data_id => $inst_data->id,
        reference_build => $reference_build,
        aligner_name => $aligner_name,
        aligner_version => $TEST_BWA_VERSION,
        samtools_version => $TEST_SAMTOOLS_VERSION,
        aligner_params => $bwa_params,
    );

    ok($alignment_result, 'defined alignment result');
    --$next_id;

    return $alignment_result;
}
