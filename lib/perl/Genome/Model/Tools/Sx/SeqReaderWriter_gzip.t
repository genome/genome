#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Genome';

use File::Temp;
use Test::More;

use_ok('Genome::Model::Tools::Sx::SeqReader') or die;
use_ok('Genome::Model::Tools::Sx::SeqWriter') or die;

my $dir = $ENV{GENOME_TEST_INPUTS} . '/Genome-Model-Tools-Sx/';
my $fasta = $dir.'/reader_writer.sam.fasta';
my $example_gzfasta = $fasta.'.gz';
ok(-s $example_gzfasta, 'example gzfasta exists') or die;
my $qual = $fasta.'.qual';
ok(-s $qual, 'qual exists') or die;
my $example_gzqual = $qual.'.gz';
ok(-s $example_gzqual, 'example gzqual exists') or die;

my $tmpdir = File::Temp::tempdir(CLEANUP => 1);
my $gzfasta = $tmpdir.'/out.fasta.gz';
my $gzqual = $tmpdir.'/out.qual.gz';

my $cmd = "gmt sx -input file=$example_gzfasta:qual_file=$example_gzqual -output type=phred:file=$gzfasta:qual_file=$gzqual"; 
my $rv = eval{ Genome::Sys->shellcmd(cmd => $cmd); };
ok($rv, 'execute sx cmd');

ok(-e $gzfasta, 'gzipped fasta created');
my $fasta_is_gzipped = eval{ Genome::Sys->shellcmd(cmd => "gunzip -t $gzfasta"); };
ok($fasta_is_gzipped, 'fasta is zipped');
my $zdiff_fasta = `zdiff $gzfasta $example_gzfasta`;
ok(!$zdiff_fasta, 'gzipped fasta match');

ok(-e $gzqual, 'gzipped qual created');
my $qual_is_gzipped = eval{ Genome::Sys->shellcmd(cmd => "gunzip -t $gzqual"); };
ok($qual_is_gzipped, 'qual is zipped');
my $zdiff_qual = `zdiff $gzqual $example_gzqual`;
ok(!$zdiff_qual, 'gzipped qual match');

#print "$tmpdir\n"; <STDIN>;
done_testing();
exit;

