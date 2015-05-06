#!/usr/bin/env genome-perl

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
    $ENV{UR_USE_DUMMY_AUTOGENERATED_IDS} = 1;
};

use above 'Genome';

use Test::More;

if (Genome::Sys->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}
else {
    plan tests => 5;
}

use_ok('Genome::Model::Tools::Gatk::BaseRecalibrator');

# Inputs
my $test_data_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-Gatk-BaseRecalibrator/v1';
my $input_bam    = "$test_data_dir/test.bam";
my $input_ref_mt = "$test_data_dir/all_sequences.MT.fa";
my $input_snv_mt = "$test_data_dir/snvs.MT.hq.vcf";

# Outputs
my $output_dir = File::Temp::tempdir('GMTGatkBaseRecalibratorXXXXX', CLEANUP => 1, TMPDIR => 1);
my $output_grp = "$output_dir/test.grp";

# Expected
my $expected_grp = "$test_data_dir/expected.grp";
my $gatk_cmd = Genome::Model::Tools::Gatk::BaseRecalibrator->create(
        max_memory                 => "2",
        version                    => 2.4,
        number_of_cpu_threads      => 1,
        input_bam                  => $input_bam,
        reference_fasta            => $input_ref_mt,
        known_sites                => [$input_snv_mt],
        output_recalibration_table => $output_grp,
);

isa_ok($gatk_cmd, 'Genome::Model::Tools::Gatk::BaseRecalibrator', "Made the command");
is($gatk_cmd->base_recalibrator_command, $gatk_cmd->base_java_command . " -T BaseRecalibrator -I $test_data_dir/test.bam -R $test_data_dir/all_sequences.MT.fa -knownSites $test_data_dir/snvs.MT.hq.vcf -o $output_grp -nct 1", 'base recalibrator command');
ok($gatk_cmd->execute, "Executed the command");
ok(system("diff $output_grp $expected_grp") == 0, "Output and expected are not different.");
done_testing();
