package Genome::Command::Tester;
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

    
    my $class;
    if (not defined $command or $command =~ /^[\w\_\:]+$/) {
        if (not defined $command) {
            my $base = $INC{"Genome.pm"} or die "Can't find Genome.pm for relative file lookup!";
            my $dir = File::Basename::dirname($base);
            my ($pkg, $file, $line) = caller();
            $class = Cwd::abs_path($file);
            $class =~ s/.t$// or die "Caller is not a .t file?: $file";
            $class =~ s/^$dir\/// or die "$file is not in directory $dir?";
            $class =~ s|/|::|g;
            UR::Object::Type->get($class) or die "cannot find class $class to go with file $file?";
        }
        else {
            $class = $command;
            UR::Object::Type->get($class) or die "cannot find class $class!";
        }
        $command = $class->help_synopsis;
        $command or die "No help_command() on $class!  Cannot autogenerate test.";
        print "SYN: $command\n";
    }
    elsif ($command =~ /^\S+\.pl\s/) {
        $class = delete $params{eventual_class};
        if (not $class) {
            fail 'when not using genome/gmt, expected "eventual_class" to be supplied so we can find the expected output directory', 
        }
    }
    else {
        my @argv = split(/\s+/,$command);
        my $params;
        my $errors;
        my $entry_class;
        
        my $cmd = shift @argv;
        if ($cmd eq 'genome') {
            $entry_class = 'Genome::Command';
        }
        elsif ($cmd eq 'gmt') {
            $entry_class = 'Genome::Model::Tools';
        }
        else {
            die "unexpected command $cmd! Expected 'genome' or 'gmt'!";
        }
        ($class, $params, $errors) = $entry_class->resolve_class_and_params_for_argv(@argv);
    }

    if (keys %params) {
        die "Unexpected args.  Expected command, results_version, and optional eventual_class: " . Data::Dumper::Dumper(%params);
    }

    my $class_as_dirname = $class;
    $class_as_dirname =~ s/::/-/g;

    my $input_dir = $ENV{GENOME_TEST_INPUTS} . "/$class_as_dirname/$results_version/inputs"; 
    my $expected_output_dir = $ENV{GENOME_TEST_INPUTS} . "/$class_as_dirname/$results_version/expected-output";

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

    # TODO: these values are clinseq-specific, though the module in overall general-purpose
    # Remove these after converting the last scripts and files into regular modules and file-based-databases.
    my $script_dir = Cwd::abs_path(File::Basename::dirname(__FILE__)) . '/../Model/ClinSeq/original-scripts/';
    my $annotation_dir = '/gscmnt/sata132/techd/mgriffit/reference_annotations/';
    
    # If the synopsis contains a /tmp/output_dir, swap it for the one we have dynamically created
    $command =~ s|/tmp/output_dir|$output_dir|g;
    $command =~ s|/tmp/input_dir|$input_dir|g;

    my $command_translated = eval qq{"$command"};
    if ($@) {
        fail "failed to translate $command.  Use \$script_dir, \$output_dir, \$annotation_dir: $@";
    }

    my $cmd = "PERL5LIB=$prefix:\$PERL5LIB genome-perl -S " . $command_translated . " 2>/dev/null 1>/dev/null";

    Genome::Sys->shellcmd(cmd => $cmd);
    ok(-d $output_dir, "output dir exists: $output_dir") or die "quitting...";

    my @diff = `diff -r --brief $expected_output_dir $output_dir`;
    is(scalar(@diff),0, "no differences between expected and actual")
        or do {
            my $tmp = $ENV{TMP} || $ENV{TMPDIR} || '/tmp';
            my $dir = $tmp . '/last-failed-' . lc($class_as_dirname);
            my $cmd = "mv $output_dir $dir";
            note($cmd);
            system $cmd;
            note("diff $expected_output_dir $dir");
        };
}

1;

