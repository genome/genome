package Genome::Model::RnaSeq::DetectFusionsResult::Chimerascan::TestHelpers;

use strict;
use warnings;

use above 'Genome';
use Test::More;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(
    test_create
    setup
);

sub test_create {
    my %parameters = @_;
    my $chimerascan_result_class = delete $parameters{chimerascan_result_class};
    my $chimerascan_version = delete $parameters{chimerascan_version};
    my $picard_version = delete $parameters{picard_version};
    my $alignment_result = delete $parameters{alignment_result};
    my $annotation_build = delete $parameters{annotation_build};

    my $class = $chimerascan_result_class;
    use_ok($class, "Can use result class");

    my %params = (
        alignment_result => $alignment_result,
        version => $chimerascan_version,
        detector_params => "--reuse-bam 0 --bowtie-version=",
        annotation_build => $annotation_build,
        picard_version => $picard_version,
        %parameters,
    );

    test_for_error($class, \%params, "You must supply a bowtie version");

    # needs --bowtie-version
    $params{'detector_params'} = "--reuse-bam 0";
    test_for_error($class, \%params, "Could not find parameter");

    # invalid bowtie version for chimerascan-vrl
    $params{'detector_params'} = "--bowtie-version 2.0.0 --reuse-bam 0", # --bowtie-version=2.0.0
              # space or = are both valid syntax  ^ here    or here ^              or here ^
    test_for_error($class, \%params, "Currently chimerascan only supports");

    # invalid value for --reuse-bam
    $params{'detector_params'} = "--bowtie-version 0.12.7 --reuse-bam bad";
    test_for_error($class, \%params, "You must specify either 1 (true) or 0 (false) for parameter");

    # should fail since -n must be an integer
    #   chimerascan_run.py: error: option -n: invalid integer value: 'a'
    $params{'detector_params'} = "--bowtie-version 0.12.7 --reuse-bam 0 -n a";
    test_for_error($class, \%params, "ERROR RUNNING COMMAND");

    # should fail because trimmed reads are shorter than segment length
    #   min_read_length after trimming: -74
    #   max_read_length after trimming: 1
    #   seed length (25) cannot be longer than read length (-74)
    #   Checking for 'bowtie-build' binary... found
    #   Checking for 'bowtie' binary... found
    #   Checking for chimerascan index directory... found
    #   Checking for bowtie index file... found
    #   Invalid run configuration, aborting.
    $params{'detector_params'} = "--bowtie-version 0.12.7 --reuse-bam 0 --trim5 100";
    test_for_error($class, \%params, "ERROR RUNNING COMMAND");
}

sub test_for_error {
    my ($class, $params, $expected_error) = @_;

    eval {
        my $result = $class->get_or_create(%{$params});
        die "failed test";
    };
    if ($@) {
        my $error_str = $@;
        chomp $error_str;
        diag "Got: \"$error_str\"";
        ok($error_str =~ m/\Q$expected_error\E/, "Crashed as expected with \"$expected_error\"");
    }
}

sub setup {
    my %parameters = @_;
    my $test_data_version = $parameters{test_data_version};
    my $chimerascan_version = $parameters{chimerascan_version};
    my $chimerascan_result_class = $parameters{chimerascan_result_class};
    my $picard_version = $parameters{picard_version};

    my $data_dir = $ENV{GENOME_TEST_INPUTS} . "/Genome-Model-RnaSeq-DetectFusionsResult-ChimerascanResult/$chimerascan_version";

    my $tophat_dir = $data_dir . "/tophat_data$test_data_version";
    die "Couldn't find tophat_dir at '$tophat_dir'" unless -d $tophat_dir;
    my @bam_files = glob(File::Spec->join($tophat_dir, "merged*.bam"));
    diag "Setting up with \n\ttest_data_version:$test_data_version\n\t" .
            "chimerascan_version:$chimerascan_version\n\t" .
            "picard_version:$picard_version\n\t" .
            "tophat_dir:$tophat_dir\n\t";

    my $t = Genome::Taxon->__define__(name => 'human');
    my $p = Genome::Individual->create(name => "test-human-patient", common_name => 'testpatient', taxon => $t);
    my $s = Genome::Sample->create(name => "test-human-patient", common_name => 'tumor', source => $p);
    my ($species_names, $versions) = @_;

    my $ref_pp = Genome::ProcessingProfile::ImportedReferenceSequence->create(name => 'test_ref_pp');

    my $ref_model = Genome::Model::ImportedReferenceSequence->create(
        name                => "test_ref_sequence_human",
        processing_profile  => $ref_pp,
        subject_class_name  => ref($s),
        subject_id          => $s->id,
    );

    my $reference_build = Genome::Model::Build::ImportedReferenceSequence->create(
        name            => "ref_sequence_$s-37",
        model           => $ref_model,
        fasta_file      => 'turkey_sammich',
        data_directory  => "$data_dir/ref",
        version         => '37',
    );

    my $annotation_model = Genome::Model::ImportedAnnotation->create(
        name => '1 chr test annotation',
        subject => $ref_model->subject,
        processing_profile => Genome::ProcessingProfile->get(name => 'imported-annotation.ensembl'),
        reference_sequence => $reference_build,
    );

    my $annotation_build = Genome::Model::Build::ImportedAnnotation->__define__(
        version => 'v1',
        model => $annotation_model,
        reference_sequence => $reference_build,
        data_directory => '/gscmnt/gc8002/info/model_data/2772828715/build125092315', #65_37j_v6
    );

    my $alignment_result = Genome::InstrumentData::AlignmentResult::Tophat->__define__(
        aligner_name => 'tophat',
        output_dir => $tophat_dir,
        reference_build_id => $reference_build->id,
        bowtie_version => '0.12.7'
    );
    $alignment_result->lookup_hash($alignment_result->calculate_lookup_hash());

    my $index_dir = File::Spec->join($data_dir, 'IndexResult');
    (my $index_class = $chimerascan_result_class) =~ s/::Result$/::Index/;
    my $index = $index_class->__define__(
        version => $chimerascan_version,
        bowtie_version => "0.12.7",
        reference_build => $reference_build,
        output_dir => $index_dir,
        annotation_build => $annotation_build,
        picard_version => $picard_version,
    );
    $index->lookup_hash($index->calculate_lookup_hash());

    Sub::Install::reinstall_sub({
        into => $chimerascan_result_class,
        as => '_staging_disk_usage',
        code => sub { return 40 * 1024 },
    });

    Sub::Install::reinstall_sub({
        into => $chimerascan_result_class,
        as => '_resolve_index_dir',
        code => sub { return $index_dir },
    });

    return $alignment_result, $annotation_build, @bam_files;
}

1;
