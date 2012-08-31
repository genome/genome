#!/usr/bin/env genome-perl

use strict;
use warnings;
use lib '/gscmnt/temp212/info/annotation/bioperl-cvs/bioperl-live';
use lib '/gscmnt/temp212/info/annotation/bioperl-cvs/bioperl-run';

use File::Basename;
#use Test::More tests => 40;
use Test::More tests => 42;

use Bio::SeqIO;

BEGIN {
    
    use_ok('GAP::Job::tRNAscan');
    use_ok('GAP::Job::RfamScan');
    use_ok('GAP::Job::RNAmmer');
    use_ok('GAP::Job::MetaRna');
   
}

my @jobs = ( );

{
    
    my $fasta = Bio::SeqIO->new(
                                -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig694.fasta',
                                -format => 'Fasta',
                            );
    
    my $seq = $fasta->next_seq();
       
    push @jobs, GAP::Job::tRNAscan->new(
                                        $seq,
                                        'bacteria',
                                        2112,
                                    );
    
    push @jobs, GAP::Job::RfamScan->new(
                                        $seq,
                                        2112,
                                    );
    
    push @jobs, GAP::Job::RNAmmer->new(
                                       $seq,
                                       'bacteria',
                                       2112,
                                   );

    push @jobs, GAP::Job::MetaRna->new(
                                       $seq,
                                       'bacteria',
                                       2112,
                                   );
    
}

foreach my $job (@jobs) {
    
    isa_ok($job, 'GAP::Job');
    can_ok($job, qw(execute execution_host));

    $job->execute();

    isnt($job->execution_host(), 'unknown');

    my $job_class = ref($job);

    my @genes = ( );

    @genes = $job->seq->get_SeqFeatures(); 
        
    foreach my $gene (@genes) {
	    isa_ok($gene, 'Bio::SeqFeature::Generic');
    }

    if (ref($job) eq 'GAP::Job::RfamScan') {
	
        foreach my $gene (@genes) {
	    
	        if ( $gene->source_tag('Infernal') and $gene->has_tag('id') ) {

		        ok($gene->has_tag('id'));
		        ok($gene->has_tag('rfam_prod'));
		    
		        my ($id)           = $gene->each_tag_value('id');
		        my ($rfam_product) = $gene->each_tag_value('rfam_prod');
		    
		        if ($id eq 'tRNA') { 
			        is($rfam_product, 'tRNA'); 
		        }
	    
	        }
	
        }

    }
    
}
