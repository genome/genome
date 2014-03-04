use strict;
use warnings;

BEGIN {
    $ENV{UR_DBI_NO_COMMIT} = 1;
};

use above 'Genome';
use Genome::Utility::Test;
use Test::More tests => 2;

my $base_dir = File::Basename::dirname(__FILE__);

subtest 'export' => sub {
    plan tests => 4;

    my $id = 2891454740;
    my $model = Genome::Model->get($id);
    ok($model, "got test model");

    my $expected = $base_dir . '/Export/Metadata.t.expected-output';

    my $intermediate_outfile = Genome::Sys->create_temp_file_path();

    my $scrubbed_outfile;
    if ($ARGV[0] && $ARGV[0] eq 'REBUILD') {
        print "\n\nRebuilding test result\n";
        $scrubbed_outfile = $expected;
        unlink $expected;
    }
    else {
        $scrubbed_outfile = Genome::Sys->create_temp_file_path();
    }

    my $result = Genome::Model::Command::Export::Metadata->execute(models => [$model], output_path => $intermediate_outfile, verbose => 1);
    ok($result, "ran");
    ok(-e $intermediate_outfile, "intermediate_outfile $intermediate_outfile exists");

    Genome::Sys->shellcmd(cmd => "grep -v 'Genome::Disk::Allocation\\|Genome::Disk::Assignment' <$intermediate_outfile | grep -v Genome::Disk::Volume >  $scrubbed_outfile");

    Genome::Utility::Test::compare_ok(
        $scrubbed_outfile,
        $expected,
        'output matches',
        filters => [
        sub{
            my $line = shift;
            return if $line =~ /total_kb/i; # ignore disk volumes
            return $line;
        },
        ],
    );
};

subtest 'import' => sub {
    plan tests => 5;

    my $input_path = $base_dir . '/Export/Metadata.t.expected-output';
    ok(-e $input_path, "input path exists: $input_path");

    my $expected_log_path = $base_dir . '/Import/Metadata.t.expected-output';
    ok(-e $expected_log_path, "expected log output file $expected_log_path exists");

    my $actual_log_path = Genome::Sys->create_temp_file_path();

    if ($ARGV[0] && $ARGV[0] eq 'REBUILD') {
        unlink $expected_log_path;
        $actual_log_path = $expected_log_path;
        warn "regenerating $expected_log_path...";
    }

    my $tx = UR::Context::Transaction->begin();
    my $result = Genome::Model::Command::Import::Metadata->execute(input_path => $input_path, log_path => $actual_log_path, verbose => 1, skip_file_db_install => 1);
    ok($result, "ran");
    $tx->rollback();
    ok(-e $actual_log_path, "actual_log_path $actual_log_path exists");

    Genome::Utility::Test::compare_ok(
        $actual_log_path,
        $expected_log_path,
        'log matches',
        filters => [
        sub{
            my $line = shift;
            return if $line =~ /total_kb/i; # ignore disk volumes
            return $line;
        },
        ],
    );
};
