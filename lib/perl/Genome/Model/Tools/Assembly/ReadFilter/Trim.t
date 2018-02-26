#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "Genome";

use Bio::SeqIO;
use Test::More;

# use
use_ok('Genome::Model::Tools::Assembly::ReadFilter::Trim') or die;

# create failures
ok(!Genome::Model::Tools::Assembly::ReadFilter::Trim->create(), 'Create w/o trim length');
ok(!Genome::Model::Tools::Assembly::ReadFilter::Trim->create(trim_length => 'all'), 'Create w/ trim length => all');
ok(!Genome::Model::Tools::Assembly::ReadFilter::Trim->create(trim_length => 0), 'Create w/ trim length => 0');

# valid create and execution
my $trimmer = Genome::Model::Tools::Assembly::ReadFilter::Trim->create(trim_length => 10);
my $io = Bio::SeqIO->new(-file => Genome::Config::get('test_inputs') . "/Genome-Model/DeNovoAssembly/velvet_solexa_build_v0.3/collated.fastq", -format => 'fastq');
my $ok = 1;
my $count = 0;
while ((my $fq = $io->next_seq) && ($count++<500)) {
    my $length = length($fq->seq);
    my $qlength = scalar @{$fq->qual};
    
    $fq = $trimmer->trim($fq);
    my $new_length = length($fq->seq);
    my $new_qlength = scalar @{$fq->qual};
    if($length != ($new_length+10)||$qlength != ($new_qlength+10))
    {
        $ok=0;
        last;
    }
}

ok($ok, "Reads trimmed successfully");

done_testing();
