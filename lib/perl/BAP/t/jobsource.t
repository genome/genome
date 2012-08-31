#!/usr/bin/env genome-perl

use strict;
use warnings;

use Test::More qw(no_plan); # or 49, or 25, or 37?
use File::Basename;

use Bio::SeqIO;
use above "BAP";

BEGIN {
    use_ok('GAP::JobSource::Composite');
    use_ok('BAP::JobSource::Genemark');
    use_ok('BAP::JobSource::Glimmer2');
    use_ok('BAP::JobSource::Glimmer3');
    use_ok('BAP::JobSource::InterGenicBlastX');
    use_ok('BAP::JobSource::Phase2BlastP');
}


#my $blast_db = '/gscmnt/temp110/analysis/blast_db/gsc_bacterial/bacterial_nr/bacterial_nr'; 
my $blast_db = '/gscmnt/gpfstest2/analysis/blast_db/gsc_bacterial/bacterial_nr/bacterial_nr';

my @job_sources = ( );


my $genemark_job_source = BAP::JobSource::Genemark->new(
                                                        Bio::SeqIO->new(
                                                                        -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.fasta',
                                                                        -format => 'Fasta',
                                                                    ),
                                                        File::Basename::dirname(__FILE__).'/data/heu_11_46.mod',
                                                    );

my $glimmer2_job_source = BAP::JobSource::Glimmer2->new(
                                                        Bio::SeqIO->new(
                                                                        -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.fasta',
                                                                        -format => 'Fasta',
                                                                    ),
                                                        File::Basename::dirname(__FILE__).'/data/glimmer2.icm',
                                                    );

my $glimmer3_job_source = BAP::JobSource::Glimmer3->new(
                                                        Bio::SeqIO->new(
                                                                        -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.fasta',
                                                                        -format => 'Fasta',
                                                                    ),
                                                        File::Basename::dirname(__FILE__).'/data/glimmer3.icm',
                                                        File::Basename::dirname(__FILE__).'/data/glimmer3.pwm',
                                                    );

isa_ok($genemark_job_source, 'BAP::JobSource::Genemark');
isa_ok($glimmer2_job_source, 'BAP::JobSource::Glimmer2');
isa_ok($glimmer3_job_source, 'BAP::JobSource::Glimmer3');

push @job_sources, GAP::JobSource::Composite->new(
                                                  $genemark_job_source,
                                                  $glimmer2_job_source,
                                                  $glimmer3_job_source,
                                              );

my $intergenic_blastx_job_source =
    BAP::JobSource::InterGenicBlastX->new(
                                          Bio::SeqIO->new(
                                                          -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.fasta',
                                                          -format => 'Fasta',
                                                      ),
                                          Bio::Tools::GFF->new(
                                                               -file => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.gff',
                                                           ),
                                          $blast_db,
                                          1,
                                      );

my $phase2_blastp_job_source = BAP::JobSource::Phase2BlastP->new(
                                                                 $blast_db,
                                                                 File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig694.pep.fasta',
                                                                1, 
                                                             );

isa_ok($intergenic_blastx_job_source, 'BAP::JobSource::InterGenicBlastX');
isa_ok($phase2_blastp_job_source,     'BAP::JobSource::Phase2BlastP');

push @job_sources, $intergenic_blastx_job_source, $phase2_blastp_job_source;

foreach my $job_source (@job_sources) {
    
    isa_ok($job_source, 'GAP::JobSource');
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

    my $job_source_class = ref($job_source);

    if ($job_source_class eq 'GAP::Job::JobSource') {
        
        my $ref = $job_source->feature_ref();

        isa_ok($ref, 'ARRAY');
        
    }
    elsif ($job_source_class eq 'BAP::JobSource::InterGenicBlastX') {

        my $collection = $job_source->feature_collection();

        isa_ok($collection, 'Bio::SeqFeature::Collection');
        
    }
    elsif ($job_source_class eq 'BAP::JobSource::Phase2BlastP') {

        my $ref = $job_source->evidence();
        
        isa_ok($ref, 'HASH');
    }
    
}

#done_testing();
unlink('error.log');

