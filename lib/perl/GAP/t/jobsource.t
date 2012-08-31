#!/usr/bin/env genome-perl

use strict;
use warnings;
use lib '/gscmnt/temp212/info/annotation/bioperl-cvs/bioperl-live';
use lib '/gscmnt/temp212/info/annotation/bioperl-cvs/bioperl-run';

use File::Basename;
use Test::More tests => 23;

use Bio::SeqIO;

BEGIN {
    use_ok('GAP::JobSource::Composite');
    use_ok('GAP::JobSource::tRNAscan');
    use_ok('GAP::JobSource::RfamScan');
    use_ok('GAP::JobSource::RNAmmer');
    use_ok('GAP::JobSource::MetaRna');
   
}

#my $blast_db = '/gscmnt/temp110/analysis/blast_db/gsc_bacterial/bacterial_nr'; 
my $blast_db = '/gscmnt/gpfstest2/analysis/blast_db/gsc_bacterial/bacterial_nr/bacterial_nr';

my @job_sources = ( );


my $trnascan_job_source = GAP::JobSource::tRNAscan->new(
                                                        Bio::SeqIO->new(
                                                                        -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.fasta',
                                                                        -format => 'Fasta',
                                                                    ),
                                                        'bacteria',
                                                    );

my $rfamscan_job_source =  GAP::JobSource::RfamScan->new(
                                                         Bio::SeqIO->new(
                                                                         -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.fasta',
                                                                         -format => 'Fasta',
                                                                     ),
                                                     );

my $rnammer_job_source  = GAP::JobSource::RNAmmer->new(
                                                       Bio::SeqIO->new(
                                                                       -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.fasta',
                                                                       -format => 'Fasta',
                                                                   ),
                                                       'bacteria',
                                                   );

my $metarna_job_source  = GAP::JobSource::MetaRna->new(
                                                       Bio::SeqIO->new(
                                                                       -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.fasta',
                                                                       -format => 'Fasta',
                                                                   ),
                                                       'bacteria',
                                                   );

isa_ok($trnascan_job_source, 'GAP::JobSource::tRNAscan');
isa_ok($rfamscan_job_source, 'GAP::JobSource::RfamScan');
isa_ok($rnammer_job_source,  'GAP::JobSource::RNAmmer');
isa_ok($metarna_job_source,  'GAP::JobSource::MetaRna');

my $job_source = GAP::JobSource::Composite->new(
                                                $trnascan_job_source,
                                                $rfamscan_job_source,
                                                $rnammer_job_source,
                                                $metarna_job_source,
                                            );



isa_ok($job_source, 'GAP::JobSource');
isa_ok($job_source, 'GAP::JobSource::Composite');
can_ok($job_source, qw(get_job finish_job fail_job));

my @jobs = ( );

while (my $job = $job_source->get_job()) {
    
    isa_ok($job, 'GAP::Job');
    can_ok($job, qw(execute execution_host));
    
    push @jobs, $job;
    
}

my $fail_job = shift @jobs;

foreach my $job (@jobs) {
    $job->execute();
    $job_source->finish_job($job, 'fake_test_result');
}

$fail_job->execute();
$job_source->fail_job($fail_job, 'fake_test_result', 'fake_test_failure');

$fail_job = $job_source->get_job();

ok(defined($fail_job));
    
$job_source->finish_job($fail_job, 'fake_test_result');

my $null_job = $job_source->get_job();

ok(!defined($null_job));

my $ref = $job_source->feature_ref();

isa_ok($ref, 'ARRAY');

unlink('error.log');

