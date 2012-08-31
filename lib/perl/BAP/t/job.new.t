#!/usr/bin/env genome-perl

use strict;
use warnings;

use above "BAP";
use File::Basename;
use Test::More qw(no_plan);  # odd, seems to bounce between 192 and 196

use Bio::SeqIO;


BEGIN {
    use_ok('BAP::Job::Genemark');
    use_ok('BAP::Job::Glimmer');
    use_ok('BAP::Job::InterGenicBlastX');
    use_ok('BAP::Job::Phase2BlastP');
    use_ok('BAP::Job::RrnaBlast');
}

my @jobs = ( );

{
    
    my $fasta = Bio::SeqIO->new(
                                -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig694.fasta',
                                -format => 'Fasta',
                            );
    
    my $seq = $fasta->next_seq();
 
    push @jobs, BAP::Job::Genemark->new(
                                        $seq,
                                        File::Basename::dirname(__FILE__).'/data/heu_11_46.mod',
                                        2112,
                                    );
    
    push @jobs, BAP::Job::Glimmer->new(
                                       'glimmer2',
                                       $seq,
                                       File::Basename::dirname(__FILE__).'/data/glimmer2.icm',
                                       undef,
                                       0,
                                       2112,
                                   );
    
    push @jobs, BAP::Job::Glimmer->new(
                                       'glimmer3',
                                       $seq,
                                       File::Basename::dirname(__FILE__).'/data/glimmer3.icm',
                                       File::Basename::dirname(__FILE__).'/data/glimmer3.pwm',
                                       0,
                                       2112,
                                   );
    
}

{

    my $fasta = Bio::SeqIO->new(
                                -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.fasta',
                                -format => 'Fasta',
                            );
    
    my $seq = $fasta->next_seq();
    my $fasta2 = Bio::SeqIO->new(
                                 -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig694.fasta',
                                 -format => 'fasta',
                                );
    my $seq2 = $fasta2->next_seq();
 
    my @features = ( );
    
    my $gff = Bio::Tools::GFF->new(
                                   -file => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig26.1.gff',
                               );
    
    while (my $feature = $gff->next_feature()) {
        push @features, $feature;
    }
    
    #my $blast_db = '/gscmnt/temp110/analysis/blast_nr/nr';
    my $blast_db = '/gscmnt/gpfstest2/analysis/blast_db/gsc_bacterial/bacterial_nr/bacterial_nr';
    my $rrna_db = '/gscmnt/278/analysis/HGMI/rRNA_testing/16s_23srnadb';
    
    push @jobs, BAP::Job::InterGenicBlastX->new(
                                                2112,
                                                $seq,
                                                \@features,
                                                $blast_db,
                                                4,
                                            );

    push @jobs, BAP::Job::RrnaBlast->new(
                                         2112,
                                         $seq2,
                                         \@features,
                                         $rrna_db,
                                         4,
                                        );
   
    
    my $pep_fasta = Bio::SeqIO->new(
                                    -file   => File::Basename::dirname(__FILE__).'/data/BACSTEFNL_Contig694.pep.fasta',
                                    -format => 'Fasta',
                                );
    
    while (my $seq = $pep_fasta->next_seq()) {
        
        push @jobs, BAP::Job::Phase2BlastP->new(
                                                $seq,
                                                $blast_db,
                                                2112,
                                                4,
                                            );
        
    }
    
}

foreach my $job (@jobs) {
    
    isa_ok($job, 'GAP::Job');
    can_ok($job, qw(execute execution_host));

    $job->execute();

    isnt($job->execution_host(), 'unknown');

    my $job_class = ref($job);

    if ($job_class eq 'BAP::Job::Phase2BlastP') {
        
        my $evidence_ref = $job->evidence();
        
        isa_ok($evidence_ref, 'HASH');
        
        ok(exists($evidence_ref->{$job->seq->display_id()}));
        
    }
    elsif ($job_class eq 'BAP::Job::RrnaBlast')
    {
        my @rnas = @{$job->genes()};
        foreach my $rna (@rnas)
        {
            isa_ok($rna, 'Bio::SeqFeature::Generic');
            print "RNA ",$rna->primary_tag()," ", $rna->source_tag(),
                  " start ", $rna->start(), "-", $rna->end(),"\n";
        }
    }
    else {

        my @genes = ( );

        if ($job_class eq 'BAP::Job::InterGenicBlastX') {
            @genes = @{$job->genes()};
        }
        else {
            @genes = $job->seq->get_SeqFeatures(); 
        }
        
        foreach my $gene (@genes) {
            
            if (ref($job) eq 'BAP::Job::Genemark') {
                isa_ok($gene,'Bio::Tools::Prediction::Gene');
            }
            else {
                isa_ok($gene, 'Bio::SeqFeature::Generic');
            }
        }
        
    }

}

done_testing();
unlink('error.log');
