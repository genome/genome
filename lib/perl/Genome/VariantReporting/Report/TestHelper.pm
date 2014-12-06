package Genome::VariantReporting::Report::TestHelper;

use strict;
use warnings;
use Genome;
use Genome::VariantReporting::Command::Wrappers::TestHelpers qw(
    compare_directories_and_files
);
use Sub::Install qw(reinstall_sub);
use Params::Validate qw(:types);
use Test::More;
use File::Copy::Recursive qw(dircopy);

use Exporter 'import';

our @EXPORT_OK = qw(
    test_report_result
);

my $FACTORY = Genome::VariantReporting::Framework::Factory->create();
my $global_data_dir = __FILE__ . ".d";

sub test_report_result {
    my %p = Params::Validate::validate(@_, {
        data_dir => {type => SCALAR},
        pkg => {type => SCALAR},
        variant_type => {type => SCALAR, default => 'snvs'},
        interpretations => {type => HASHREF},
    });

    my $plan = plan($p{data_dir});

    my $process = Genome::Process::Test::Process->create(statement => 'test');
    my $cmd = Genome::VariantReporting::Framework::GenerateReport->create(
        report_name => $p{pkg}->name,
        input_vcf => input_vcf(),
        plan_json => $plan->as_json,
        variant_type => $p{variant_type},
        process_id => $process->id,
    );

    mock_interpret_entry($p{interpretations});
    ok($cmd->execute, "Report command successfully executed");

    my $result = $cmd->output_result;
    if ($ENV{GENERATE_TEST_DATA}) {
        local $File::Copy::Recursive::RMTrgDir = 1;
        dircopy($result->output_dir, expected_directory($p{data_dir}));
    }
    compare_directories_and_files($result->output_dir,
        expected_directory($p{data_dir}));
}


sub mock_interpret_entry {
    my $all_interpretations = shift;

    while (my ($interpreter_name, $interpretations) = each %$all_interpretations) {
        reinstall_sub({
            into => $FACTORY->get_class("interpreters", $interpreter_name),
            as => "interpret_entry",
            code => sub {
                return %$interpretations;
            },
        });
    }
}

sub expected_directory {
    my $data_dir = shift;
    return File::Spec->join($data_dir, 'expected');
}

sub input_vcf {
    return File::Spec->join($global_data_dir, "input.vcf.gz");
}

sub plan_file {
    my $data_dir = shift;
}

sub plan {
    my $data_dir = shift;
    my $plan_file = File::Spec->join($data_dir, 'plan.yaml');
    return Genome::VariantReporting::Framework::Plan::MasterPlan->
        create_from_file($plan_file);
}


1;
