package Genome::VariantReporting::Component::Adaptor::TestHelpers;

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Exporter 'import';
use Test::Deep;
use Sub::Install qw(reinstall_sub);

our @EXPORT_OK = qw(
    test_accessor
    setup_results
);

{
    package Genome::VariantReporting::TestAdaptor;

    use strict;
    use warnings FATAL => 'all';
    use Genome;

    class Genome::VariantReporting::TestAdaptor {
        is => 'Genome::VariantReporting::Component::Adaptor',
    };

    sub name {
        'test';
    }

    sub resolve_plan_attributes {
        return;
    }

    1;
}

sub test_accessor {
    my ($build, $bam_result1, $bam_result2) = @_;

    my $cmd = Genome::VariantReporting::TestAdaptor->create(
        build_id => $build->id,
        variant_type => 'snvs',
        plan_json => 'unused',
    );
    ok($cmd, "Command created correctly");
    ok($cmd->execute, "Command executed successfully");
    cmp_bag([$cmd->bam_results], [$bam_result1, $bam_result2], "Bam results set as expected");
}

sub setup_results {
    my $bam_result1 = Genome::InstrumentData::AlignmentResult::Merged->__define__();
    my $bam_result2 = Genome::InstrumentData::AlignmentResult::Merged->__define__();

    return ($bam_result1, $bam_result2);
}
