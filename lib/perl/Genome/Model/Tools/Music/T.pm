package Genome::Model::Tools::Music::T;
use Exporter 'import';

use strict;
use warnings;

use above 'Genome';
use Genome::Model::Tools::Music;
use Test::More;

our @EXPORT_OK = qw(run_test_case $input_dir $actual_output_dir);

# figure out where the test inputs are and expected outputs
# the package with this data is a dependency so this should work when deployed externally
my $test_data_dir = Genome::Sys->dbpath( 'genome-music-testdata', $Genome::Model::Tools::Music::VERSION) ;
my $expected_output_dir = $test_data_dir . '/expected_outputs/';
our $input_dir = $test_data_dir . '/inputs';
our $actual_output_dir;

# Decide where output goes depending on how this test was invoked
if( @ARGV ) {
    # Override output dir
    if( $ARGV[0] eq '--regenerate' ) {
        # Regenerate all output files as the new "correct" answer
        $actual_output_dir = $expected_output_dir;
    }
    else {
        # Use the dir the user specifies (for testing since tempdirs get destroyed)
        $actual_output_dir = shift @ARGV;
        mkdir $actual_output_dir unless -d $actual_output_dir;
        unless( -d $actual_output_dir ) {
            die "failed to create directory $actual_output_dir: $!";
        }
    }
}
else {
    # Use a temp dir if none were specified as arguments
    $actual_output_dir = Genome::Sys->create_temp_directory( "music" );
};

sub run_test_case {
    my %case = @_;

    if ($case{skip}) {
        plan skip_all => $case{skip};
        return;
    }

    # Pre-determine how many tests will run so the test harness knows if we exit early
    my $tests = 3; # These are the 3 tests looking for the test-data
    my $expect = $case{expect};

# TODO huh?
#    unless( $expect ) {
#        warn( "No expected output defined for: " . $case{run} );
#        next;
#    }

    $tests += 2 + ( scalar( @$expect ) * 2 );
    plan tests => $tests;

    ok( -d $test_data_dir, 'Directory with genome-music-testdata exists') or die;
    ok( -d $input_dir, 'Directory with genome-music-testdata inputs exists') or die;
    ok( -d $expected_output_dir, 'Directory with genome-music-testdata expected outputs exists') or die;

    my $cmd = $case{run};

    note( "use case: $cmd" );
    if( my $msg = $case{skip} ) {
        note "SKIPPING: $case{skip}\n";
        next;
    }

    # Make subdirs for the output if needed
    for my $expect_file ( @$expect ) {
        my $actual_full_path = $actual_output_dir . '/' . $expect_file;
        my $dir = $actual_full_path;
        use File::Basename;
        $dir = File::Basename::dirname( $dir );
        Genome::Sys->create_directory( $dir );
    }

    # Execute
    my @args = split( ' ', $cmd );
    my $exit_code = eval {
        Genome::Model::Tools->_execute_with_shell_params_and_return_exit_code( @args );
    };

    ok( !$@, "ran without crashing" ) or diag $@;
    is( $exit_code, 0, "returned a zero (good) exit code" ) or next;

    # Compare results
    for my $expect_file ( @$expect ) {
        my $expect_full_path = $expected_output_dir . '/'. $expect_file;
        my $actual_full_path = $actual_output_dir . '/' . $expect_file;

        ok( -e $actual_full_path, "has expected output file $expect_file" ) or next;

        my @diff = `diff -u $expect_full_path $actual_full_path`;
        my $diff_output;
        if( @diff > 20 ) {
            $diff_output = join( "\n", @diff[0..19] ) . "\ndiff output truncated.";
        }
        else {
            $diff_output = join( "\n", @diff );
        }
        is( scalar( @diff ), 0, "matches expectations for file $expect_file" )
            or diag( "\$ diff $expect_full_path $actual_full_path\n" . $diff_output . "\n" );
    }
}

1;
