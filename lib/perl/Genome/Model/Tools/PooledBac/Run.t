#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';
use Test::More;
use Storable qw();
use Genome::Utility::Test qw(capture_ok compare_ok);

if (Genome::Sys->arch_os ne 'x86_64') {
    plan skip_all => 'requires 64-bit machine';
}

use_ok( 'Genome::Model::Tools::PooledBac::Run' ) or die;

my $version = 2;
my $test_dir = Genome::Config::get('test_inputs') . '/Genome-Model-Tools-PooledBac/Run_v'.$version;
ok( -d $test_dir, 'Test suite dir exists' ) or die;

my $run_dir = Genome::Sys->create_temp_directory();
ok( -d $run_dir, 'Test run dir created' );

# Update Storage file, Contig1.ci, which has a frozen path in it.
# Setting $Storable::canonical may help prevent side effects from sort changes
# due to rebuilding Contig1.ci.
my $temp_pooled_bac_dir = Genome::Sys->create_temp_directory();
capture_ok(['rsync', '-av', "$test_dir/pooled_bac/", "$temp_pooled_bac_dir/pooled_bac/"]) or die;
my $dft_path = "$temp_pooled_bac_dir/pooled_bac/consed/edit_dir/ace.dft";
my $contig_path = "$dft_path.idx/Contig1.ci";
my $contig = Storable::retrieve($contig_path);
$contig->{file_name} = $dft_path;
$Storable::canonical = 1;
Storable::store($contig, $contig_path);

my %params = (
    ace_file_name => 'ace.dft',
    pooled_bac_dir =>  "$temp_pooled_bac_dir/pooled_bac/",
    ref_seq_file => $test_dir.'/ref_seq/ref_small.txt',
    project_dir => $run_dir,
);

my $tool = Genome::Model::Tools::PooledBac::Run->create( %params );
ok( $tool, 'Created tool' ) or die;
ok( $tool->execute, 'Successfully executed tool' ) or die;

#compare output files: ace
my $project_name = 'H_GD-332K02';
my $project_ace = 'H_GD-332K02.fasta.screen.ace';
my $old_ace = $test_dir."/project/$project_name/edit_dir/$project_ace";
ok( -s $old_ace, 'Old ace exists' );
my $new_ace = $run_dir."/$project_name/edit_dir/$project_ace";
ok( -s $new_ace, 'New ace created' );
compare_ok($new_ace, $old_ace, 'Ace files match');

#compare output files: report
my @report_files = qw/
orphan_contigs
contigs_with_multiple_hits
contigs_that_link_to_matching_contigs
ambiguous_matching_contigs
contigs_only_consensus
assembly_size_report
complete_contig_list_with_orphan_contigs
complete_contig_list
contig_size_report
matching_contigs
matching_contigs
/;

for my $file ( @report_files ) {
    ok( -e $test_dir."/project/reports/$file", "Test dir $file exists" );
    ok( -e $run_dir."/reports/$file", "New $file created" );
    compare_ok("$run_dir/reports/$file", "$test_dir/project/reports/$file", "$file files match");
}

#<STDIN>;

done_testing();
