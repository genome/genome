#!/usr/bin/env genome-perl

use strict;
use warnings;

use above 'Workflow';

use Workflow::Simple;


my $xml_file   = $ARGV[0] || 'data/repeatmasker_inner.xml';

my $output = run_workflow_lsf(
                              $xml_file,
                              'fasta file'     => 'data/BACSTEFNL_Contig694.fasta',
                              'repeat library' => '/gsc/var/lib/repeat/Trichinella_pseudospiralis_1.0_080103.rep', 
                              'species'        => 'elegans', 
                          );


print Data::Dumper->new([$output,\@Workflow::Simple::ERROR])->Dump;

