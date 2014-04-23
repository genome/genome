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

{
    package Genome::Annotation::TestAdaptor;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::Annotation::TestAdaptor {
        is => 'Genome::Annotation::AdaptorBase',
    };

    sub resolve_plan_attributes {
        return;
    }

    1;
}

sub test_accessor {
    my ($build, $bam_result1, $bam_result2, $snv_vcf_result, $indel_vcf_result, $add_vcf_results) = @_;

    # This will run in two modes, with vcf results and without vcf results added
    if ($add_vcf_results) {
        add_vcf_results($snv_vcf_result, $indel_vcf_result);
    } else {
        $snv_vcf_result = undef;
        $indel_vcf_result = undef;
    }

    my $cmd = Genome::Annotation::TestAdaptor->create(build_id => $build->id, variant_type => 'snvs');
    ok($cmd, "Command created correctly");
    ok($cmd->execute, "Command executed successfully");
    cmp_bag([$cmd->bam_results], [$bam_result1, $bam_result2], "Bam results set as expected");
    is_deeply($cmd->output_result, $snv_vcf_result, "Snvs vcf result is as expected");

    # now again for indels
    $cmd = Genome::Annotation::TestAdaptor->create(build_id => $build->id, variant_type => 'indels');
    $cmd->execute;
    is_deeply($cmd->output_result, $indel_vcf_result, "Indels vcf result is as expected");
}

sub test_accessors_without_vcf_results {
    test_accessor(@_, 0);
}

sub test_accessors_with_vcf_results {
    test_accessor(@_, 1);
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
