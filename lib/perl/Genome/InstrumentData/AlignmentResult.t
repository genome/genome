#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
}

use strict;
use warnings;

use above 'Genome';

use Genome::InstrumentData::InstrumentDataTestObjGenerator;
use Path::Class;
use File::Slurp qw(read_file write_file);
use Test::More;

my $bam_path = $ENV{GENOME_TEST_INPUTS} . '/Genome-InstrumentData-AlignmentResult-Bwa/input.bam';

use_ok('Genome::InstrumentData::AlignmentResult');

# Test AR Class
class Genome::InstrumentData::AlignmentResult::Tester {
    is => 'Genome::InstrumentData::AlignmentResult',
};
my $iar_id;
sub Genome::InstrumentData::AlignmentResult::Tester::_run_aligner { 
    my $self = shift;

    # 'Align' [copy] bam for down stream ops
    my $copy_ok = Genome::Sys->copy_file($bam_path, $self->temp_staging_directory.'/all_sequences.bam');
    ok($copy_ok, 'copied bam');

    return 1;
};
sub Genome::InstrumentData::AlignmentResult::Tester::aligner_params_for_sam_header { 'aln cmd -x' };
sub Genome::InstrumentData::AlignmentResult::Tester::estimated_kb_usage { 0 };
sub Genome::InstrumentData::AlignmentResult::Tester::fillmd_for_sam { 0 };
sub Genome::InstrumentData::AlignmentResult::Tester::requires_fastqs_to_align { 0 };

my $next_id = -1337;
sub make_alignment_result {

    my $inst_data = Genome::InstrumentData::InstrumentDataTestObjGenerator::create_solexa_instrument_data($bam_path);
    ok($inst_data, 'create inst data');
    my $reference_model = Genome::Model::ImportedReferenceSequence->get(name => 'TEST-human');
    ok($reference_model, "got reference model");
    my $reference_build = $reference_model->build_by_version('1');
    ok($reference_build, "got reference build");

    my $alignment_result = Genome::InstrumentData::AlignmentResult::Tester->create(
        id => $next_id,
        instrument_data_id => $inst_data->id,
        reference_build => $reference_build,
        aligner_name => 'tester',
        aligner_version => '1',
        aligner_params => '',
    );

    ok($alignment_result, 'defined alignment result');
    isa_ok($alignment_result, 'Genome::InstrumentData::AlignmentResult::Tester');
    --$next_id;

    return $alignment_result;
}

subtest 'Temporary Input Files Queue Usage' => sub {
    my $ar = make_alignment_result();
    my $tmpdir = Path::Class::Dir->new($ar->temp_scratch_directory);

    diag("Creating test file system in : $tmpdir");
    my ($root, @files) = create_test_file_system($tmpdir);

    is(-d "$root", 1, "'$root' exists on file system");

    my @items = $ar->temporary_input_files_queue();
    ok(@items == 0, "queue is empty");

    ok($ar->add_to_temporary_input_files_queue(@files, $root), "adding temp files to queue");
    @items = $ar->temporary_input_files_queue();
    ok(@items == 3, "there are 3 items in the queue");

    ok($ar->show_temporary_input_files_queue(), "showing temp files in queue");

    ok($ar->clear_temporary_input_files_queue(), "clearing out temp files in queue");
    @items = $ar->temporary_input_files_queue();
    ok(@items == 0, "the queue is again empty");

    is(-d "$root", undef, "'$root' no longer exists on file system");
    done_testing();
};

subtest "construct groups file" => sub {
    my $alignment_result = make_alignment_result();

    my $header_extra = $alignment_result->_sam_header_extra();
    my @keys = sort keys %$header_extra;
    is_deeply(\@keys, ["PG", "RG"], "got RG and PG tags for sam header")
        or "header extra: " . diag($header_extra);

    my $inst_data = $alignment_result->instrument_data;
    my $pu_tag = sprintf("%s.%s", $inst_data->run_identifier, $inst_data->subset_name);

    my %expected_rg = (
        ID => sprintf("%d", $inst_data->id),
        PL => "illumina",
        PU => $pu_tag,
        PI => "0",
        DS => "paired end",
        LB => $inst_data->library->name,
        SM => $inst_data->sample->name,
        CN => "WUGSC",
        DT => $inst_data->run_start_date_formatted,
        );

    my %expected_pg = (
        ID => sprintf("%d", $inst_data->id),
        VN => "1",
        CL => "aln cmd -x",
        );

    my @rg_fields = split("\t", $header_extra->{RG});
    my @pg_fields = split("\t", $header_extra->{PG});
    my $tmp = shift @rg_fields;
    is($tmp, "\@RG", "SAM read group line begins with \@RG");
    $tmp = shift @pg_fields;
    is($tmp, "\@PG", "SAM program line begins with \@PG");

    my %actual_rg = map {split ":", $_, 2} @rg_fields;
    my %actual_pg = map {split ":", $_, 2} @pg_fields;

    is_deeply(\%actual_rg, \%expected_rg, "RG tag data is correct");
    is_deeply(\%actual_pg, \%expected_pg, "PG tag data is correct");

    my $tmpdir = Genome::Sys->create_temp_directory;
    my $path = sprintf("%s/groups.sam", $tmpdir);

    # Write some data to test that the file is appended to
    my $data = "append stuff to me";
    write_file($path, "$data\n");
    ok($alignment_result->construct_groups_file($path), "Created groups.sam file");
    my @lines = read_file($path);
    chomp @lines;
    is(scalar @lines, 3, "Got 3 lines as expected");
    is($lines[0], $data, "File was appended to, not overwritten");
    my @actual_tags = sort @lines[1..2];
    my @expected_tags = sort values %{$alignment_result->_sam_header_extra()};
    is_deeply(\@actual_tags, \@expected_tags, "groups file looks good");
};

subtest "construct qc" => sub {
    my $alignment_result = make_alignment_result();

    my $qc_result = Genome::InstrumentData::AlignmentResult::Merged::BamQc->__define__(
        alignment_result_id => $alignment_result->id
    );
    ok($qc_result, 'define qc result for alienment result');
    ok($alignment_result->add_user(user => $qc_result, label => 'uses'),
        'add qc result as user of alignment result');

    # delete
    ok($alignment_result->delete, 'delete');
    ok(ref($qc_result) eq 'UR::DeletedRef', 'deleted qc result');
};

done_testing();

sub create_test_file_system {
    my $tmp = shift;

    my $root = $tmp->subdir('anonymous0');
    ok($root->mkpath, "Created root path: $root");

    my @files = ();
    for my $i (1..2) {
        my $name = join('_', 's', 'unknown', $i, 'sequence') . '.txt';
        my $f = $root->file($name);
        ok($f->touch, "Creating file: $f");
        push(@files, $f);
    }

    return $root, @files;
}
