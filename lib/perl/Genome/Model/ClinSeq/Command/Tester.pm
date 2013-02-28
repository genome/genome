package Genome::Model::ClinSeq::Command::Tester;
use strict;
use warnings;
use Test::More;
use Cwd;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(run_and_diff);

sub run_and_diff {
    my %params = @_;
    my $command = delete $params{command};
    my $results_version = delete $params{results_version};

    my $expected_output_dir = $ENV{GENOME_TEST_INPUTS} || '';
    
    my $class;
    if ($command =~ /^\S+\.pl\s/) {
        $class = delete $params{eventual_class};
        if (not $class) {
            fail 'when not using genome/gmt, expected "eventual_class" to be supplied so we can find the expected output directory', 
        }
    }
    else {
        my @argv = split(/\s+/,$command);
        print Data::Dumper::Dumper(\@argv);
        my $params;
        my $errors;
        ($class, $params, $errors) = $class->resolve_class_and_params_for_argv(@argv);
    }

    if (keys %params) {
        die "Unexpected args.  Expected command, results_version, and optional eventual_class: " . Data::Dumper::Dumper(%params);
    }

    my $class_as_dirname = $class;
    $class_as_dirname =~ s/::/-/g;

    $expected_output_dir .= "/$class_as_dirname/$results_version/expected-output";

    my $output_dir;
    if (@ARGV and $ARGV[0] eq 'REBUILD') {
        # the REBUILD flag on the command-line will generate new comparison data instead of using it
        # when this set the test will always pass because the output_dir and expected_output_dir are the same
        if (-d $expected_output_dir) {
            die "directory exists ...update to a new directory first!: $expected_output_dir";
        }
        Genome::Sys->create_directory($expected_output_dir);
        $output_dir = $expected_output_dir;
    }
    else {
        $output_dir = Genome::Sys->create_temp_directory();
    }
    
    ok($expected_output_dir, "resolved expected output dir: $expected_output_dir") or die "quitting..."; 
    
    unless(-e $expected_output_dir) {
        die "*\n* >>> cannot find expected output directory: run with REBUILD to generate it\n*\n";
    };

    my $prefix = UR::Util::used_libs_perl5lib_prefix();
    ok($prefix, "using local directory $prefix") or die "quitting..";
    my $script_dir = Cwd::abs_path(File::Basename::dirname(__FILE__)) . '/../original-scripts/';
    my $annotation_dir = '/gscmnt/sata132/techd/mgriffit/reference_annotations/';
    my $command_translated = eval qq{"$command"};
    if ($@) {
        fail "failed to translate $command.  Use \$script_dir, \$output_dir, \$annotation_dir: $@";
    }

    my $cmd = "PERL5LIB=$prefix:\$PERL5LIB perl " . $command_translated;

    Genome::Sys->shellcmd(cmd => $cmd);
    ok(-d $output_dir, "output dir exists: $output_dir") or die "quitting...";

    my @diff = `diff -r --brief $expected_output_dir $output_dir`;
    is(scalar(@diff),0, "no differences between expected and actual")
        or do {
            my $tmp = $ENV{TMP} || $ENV{TMPDIR} || '/tmp';
            my $dir = $tmp . '/last-failed-' . lc($class_as_dirname);
            my $cmd = "mv $output_dir $dir";
            print STDOUT $cmd,"\n";
            system $cmd;
        };
}

1;
