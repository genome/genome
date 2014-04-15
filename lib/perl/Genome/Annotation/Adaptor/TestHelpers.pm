package Genome::Annotation::Adaptor::TestHelpers;

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Exporter 'import';
use Test::Deep;
use Sub::Install qw(reinstall_sub);

our @EXPORT_OK = qw(
    test_accessors_without_vcf_results
    test_accessors_with_vcf_results
    setup_results
);

sub test_accessor {
    my ($pkg, $build, $add_vcf_results, $bam_result1, $bam_result2, $snv_vcf_result, $indel_vcf_result) = @_;

    # This will run in two modes, with vcf results and without vcf results added
    if ($add_vcf_results) {
        add_vcf_results($snv_vcf_result, $indel_vcf_result);
    } else {
        $snv_vcf_result = undef;
        $indel_vcf_result = undef;
    }

    my $cmd = $pkg->create(build => $build);
    ok($cmd->isa("$pkg"), "Command created correctly");
    ok($cmd->execute, "Command executed successfully");
    cmp_bag([$cmd->bam_results], [$bam_result1, $bam_result2], "Bam results set as expected");
    is_deeply($cmd->annotation_build, $build->annotation_build, "Annotation build set as expected");

    is_deeply($cmd->snv_vcf_result, $snv_vcf_result, "Snvs vcf result is as expected");
    is_deeply($cmd->indel_vcf_result, $indel_vcf_result, "Indel vcf result is as expected");
}

sub test_accessors_without_vcf_results {
    my ($pkg, $build, $bam_result1, $bam_result2, $snv_vcf_result, $indel_vcf_result) = @_;
    test_accessor($pkg, $build, 0, $bam_result1, $bam_result2, $snv_vcf_result, $indel_vcf_result);
}

sub test_accessors_with_vcf_results {
    my ($pkg, $build, $bam_result1, $bam_result2, $snv_vcf_result, $indel_vcf_result) = @_;
    test_accessor($pkg, $build, 1, $bam_result1, $bam_result2, $snv_vcf_result, $indel_vcf_result);
}

sub setup_results {
    my $bam_result1 = Genome::InstrumentData::AlignmentResult::Merged->__define__();
    my $bam_result2 = Genome::InstrumentData::AlignmentResult::Merged->__define__();

    my $snv_vcf_result = Genome::Model::Tools::DetectVariants2::Result::Vcf::Combine->__define__;
    my $indel_vcf_result = Genome::Model::Tools::DetectVariants2::Result::Vcf::Combine->__define__;

    return ($bam_result1, $bam_result2, $snv_vcf_result, $indel_vcf_result);
}

sub add_vcf_results {
    my ($snv_vcf_result, $indel_vcf_result) = @_;
    reinstall_sub({
        into => "Genome::Model::Build::RunsDV2",
        as => "get_detailed_snvs_vcf_result",
        code => sub { my $self = shift;
                      return $snv_vcf_result;
        },
    });
    reinstall_sub({
        into => "Genome::Model::Build::RunsDV2",
        as => "get_detailed_indels_vcf_result",
        code => sub { my $self = shift;
                      return $indel_vcf_result;
        },
    });
}
