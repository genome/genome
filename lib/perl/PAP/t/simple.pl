#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Workflow';

use Workflow::Simple;


#print Data::Dumper->new([\%ENV])->Dump;
print Data::Dumper->new([\@INC])->Dump;

my $xml_file   = $ARGV[0] || 'data/pap_outer.xml';
my $fasta_file = $ARGV[1] || 'data/B_coprocola.fasta';

my $output = run_workflow_lsf(
                              $xml_file,
                              'fasta file'       => $fasta_file,
                              'chunk size'       => 10,
                              'dev flag'         => 1,
                              'biosql namespace' => 'MGAP',
                              'gram stain'       => 'negative',
                              'report save dir'  => '/gscmnt/temp212/info/annotation/PAP_testing/blast_reports',
                          );


print Data::Dumper->new([$output,\@Workflow::Simple::ERROR])->Dump;

