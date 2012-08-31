#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More;
use above "Genome";
require File::Compare;

use_ok( 'Genome::Model::Tools::Newbler::StandardOutputs' ) or die;

#test suite
my $version = 'v5';
my $test_suite = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Newbler/StandardOutputs-'.$version;
ok( -d $test_suite, "Test suite dir exists" );

#temp test dir
my $temp_dir = Genome::Sys->create_temp_directory();
ok( -d $temp_dir, "Created temp directory" );
Genome::Sys->create_directory( $temp_dir.'/consed' );
Genome::Sys->create_directory( $temp_dir.'/consed/edit_dir' );
ok( -d $temp_dir.'/consed/edit_dir', "Created consed/edit_dir in temp_dir" );

#link files for test
my @files_to_link = (
'2869511846-input.fastq',
'454AlignmentInfo.tsv',
'454AllContigs.fna',
'454AllContigs.qual',
'454ContigGraph.txt',
'454LargeContigs.fna',
'454LargeContigs.qual',
'454NewblerMetrics.txt',
'454NewblerProgress.txt',
'454ReadStatus.txt',
'454TrimStatus.txt',
'/consed/edit_dir/454Contigs.ace.1'
);

for my $file ( @files_to_link ) {
    ok( -s "$test_suite/$file", "Test suite $file file exists" );
    symlink( "$test_suite/$file", "$temp_dir/$file" );
    ok( -l "$temp_dir/$file", "Linked $file file to temp test dir" );
}

#create/execute tool
my $create = Genome::Model::Tools::Newbler::StandardOutputs->create(
    assembly_directory => $temp_dir,
    min_contig_length => 20,
);
ok( $create, "Created tool" );
ok( $create->execute, "Successfully executed tool" );

#compare output files
my @files_to_compare = qw/
Pcap.454Contigs.ace
gap.txt
supercontigs.fasta
supercontigs.agp
contigs.quals
contigs.bases
readinfo.txt
reads.placed
reads.unplaced
reads.unplaced.fasta
/;

for my $file ( @files_to_compare ) {
    ok( -s "$test_suite/consed/edit_dir/$file", "Test suite $file file exists" );
    ok( -s "$temp_dir/consed/edit_dir/$file", "Test created $file file" );
    ok( File::Compare::compare("$test_suite/consed/edit_dir/$file","$temp_dir/consed/edit_dir/$file")==0, "$file files match" );
}

#make sure no unexpected file created
for my $file ( glob( "$temp_dir/consed/edit_dir/*" ) ) {
    my $base_name = File::Basename::basename( $file );
    next if $base_name =~ /^454Contigs.ace.1$/; #an input file
    ok ( grep (/^$base_name$/, @files_to_compare), "Got file $base_name as expected" );
}

#<STDIN>;

done_testing();

exit;
