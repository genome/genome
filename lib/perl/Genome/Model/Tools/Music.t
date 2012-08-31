#!/usr/bin/env genome-perl
use strict;
use warnings;
use above 'Genome';
use Genome::Model::Tools::Music;
use Test::More;

# figure out where the test inputs are and expected outputs
# the package with this data is a dependency so this should work when deployed externally
my $test_data_dir = Genome::Sys->dbpath( 'genome-music-testdata', $Genome::Model::Tools::Music::VERSION) ;
my $input_dir = $test_data_dir . '/inputs';
my $expected_output_dir = $test_data_dir . '/expected_outputs/';
my $actual_output_dir;

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

# Use-cases and expected outputs
my @cases = (
    {
        skip => 1, # Skip this tool, R arbitrarily changes precision causing expected output to differ
        run => "music clinical-correlation\n"
             . " --numeric-clinical-data-file $input_dir/numeric_clinical_data\n"
             . " --categorical-clinical-data-file $input_dir/categorical_clinical_data\n"
             . " --bam-list $input_dir/bam_list\n"
             . " --maf-file $input_dir/ucec_test.maf\n"
             . " --output-file $actual_output_dir/clin_correlations\n"
             . " --numerical-data-test-method wilcox\n"
             . " --genetic-data-type gene\n"
             . " --noskip-non-coding",
        expect => [
            'clin_correlations.categorical.csv',
            'clin_correlations.numeric.csv'
        ],
    },
    {
        skip => 1, # Skip this tool, R arbitrarily changes precision causing expected output to differ
        run => "music survival\n"
             . " --numeric-clinical-data-file $input_dir/numeric_clinical_data\n"
             . " --bam-list $input_dir/bam_list\n"
             . " --maf-file $input_dir/ucec_test.maf\n"
             . " --output-dir $actual_output_dir/survival\n"
             . " --phenotypes-to-include TP53,PIK3CA",
        expect => [
            'survival/survival_analysis_data_matrix.csv',
            'survival/survival_analysis_test_results.csv'
        ],
    },
    {
        skip => 1, # Skip this tool, because it's crap and needs to be rewritten
        run => "music cosmic-omim \n"
             . " --maf-file $input_dir/ucec_test.maf\n"
             . " --output-file $actual_output_dir/cosmic_omim.maf\n"
             . " --no-verbose",
        expect => [
            'cosmic_omim.maf'
        ],
    },
    {
        run => "music pfam\n"
             . " --maf-file $input_dir/ucec_test.maf\n"
             . " --output-file $actual_output_dir/pfam.maf",
        expect => [
            'pfam.maf'
        ],
    },
    {
        run => "music proximity\n"
             . " --maf-file $input_dir/ucec_test.maf\n"
             . " --output-dir $actual_output_dir\n"
             . " --noskip-non-coding",
        expect => [
            'proximity_report'
        ],
    },
    {
        run => "music path-scan\n"
             . " --bam-list $input_dir/bam_list\n"
             . " --maf-file $input_dir/ucec_test.maf\n"
             . " --pathway-file $input_dir/kegg_db_120910\n"
             . " --gene-covg-dir $input_dir/gene_covgs\n"
             . " --output-file $actual_output_dir/sm_pathways\n"
             . " --bmr 1.8E-06",
        expect => [
            'sm_pathways',
            'sm_pathways_detailed'
        ],
    },
    {
        run => "music smg\n"
             . " --gene-mr-file $input_dir/gene_mrs\n"
             . " --output-file $actual_output_dir/smgs\n"
             . " --max-fdr 0.55",
        expect => [
            'smgs',
            'smgs_detailed'
        ],
    },
    {
        skip => 1, # Skip this tool, R arbitrarily changes precision causing expected output to differ
        run => "music mutation-relation\n"
             . " --bam-list $input_dir/bam_list\n"
             . " --maf-file $input_dir/ucec_test.maf\n"
             . " --gene-list $input_dir/smgs\n"
             . " --output-file $actual_output_dir/mutation_relations\n"
             . " --mutation-matrix-file $actual_output_dir/mutation_matrix\n"
             . " --permutations 500\n"
             . " --noskip-non-coding",
        expect => [
            'mutation_relations',
            'mutation_matrix'
        ],
    },
);

# Pre-determine how many tests will run so the test harness knows if we exit early
my $tests = 3; # These are the 3 tests looking for the test-data
for my $case ( @cases ) {
    next if $case->{skip};
    my $expect = $case->{expect};
    unless( $expect ) {
        warn( "No expected output defined for: " . $case->{run} );
        next;
    }
    $tests += 2 + ( scalar( @$expect ) * 2 );
}
plan tests => $tests;

ok( -d $test_data_dir, 'Directory with genome-music-testdata exists') or die;
ok( -d $input_dir, 'Directory with genome-music-testdata inputs exists') or die;
ok( -d $expected_output_dir, 'Directory with genome-music-testdata expected outputs exists') or die;

# Run each case
my $n = 0;
for my $case ( @cases ) {
    my $cmd = $case->{run};
    my $expect = $case->{expect};

    ++$n;
    note( "use case $n: $cmd" );
    if( my $msg = $case->{skip} ) {
        note "SKIPPING: $case->{skip}\n";
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

    ok( !$@, " case $n ran without crashing" ) or diag $@;
    is( $exit_code, 0, " case $n ran returned a zero (good) exit code" ) or next;

    # Compare results
    for my $expect_file ( @$expect ) {
        my $expect_full_path = $expected_output_dir . '/'. $expect_file;
        my $actual_full_path = $actual_output_dir . '/' . $expect_file;

        ok( -e $actual_full_path, " case $n has expected output file $expect_file" ) or next;

        my @diff = `diff -u $expect_full_path $actual_full_path`;
        my $diff_output;
        if( @diff > 20 ) {
            $diff_output = join( "\n", @diff[0..19] ) . "\ndiff output truncated.";
        }
        else {
            $diff_output = join( "\n", @diff );
        }
        is( scalar( @diff ), 0, " case $n matches expectations for file $expect_file" )
            or diag( "\$ diff $expect_full_path $actual_full_path\n" . $diff_output . "\n" );
    }
}
