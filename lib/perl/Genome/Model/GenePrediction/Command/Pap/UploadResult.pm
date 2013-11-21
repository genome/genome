#$Id$

package Genome::Model::GenePrediction::Command::Pap::UploadResult;

use strict;
use warnings;


use Bio::DB::BioDB;
use Bio::DB::Query::BioQuery;


class Genome::Model::GenePrediction::Command::Pap::UploadResult {
    is  => ['Command::V1'],
    has => [
        dev_flag         => {
                             is  => 'SCALAR',
                             doc => 'if true, connect to dev biosql',
                             is_input => 1,
                            },
        biosql_namespace => { 
                             is  => 'SCALAR', 
                             doc => 'biosql namespace',
                             is_input => 1,
                            },
        bio_seq_features => { 
                              is  => 'ARRAY',
                              doc => 'array of Bio::Seq::Feature' ,
                             is_input => 1,
                            },
        lsf_queue => { is_param => 1, default_value => $ENV{GENOME_LSF_QUEUE_BUILD_WORKER},},
        lsf_resource => {is_param => 1, default_value => 'rusage[tmp=100]',},
    ],
};


sub sub_command_sort_position { 10 }

sub help_brief {
    "Store input gene predictions in the BioSQL schema using the specified namespace";
}

sub help_synopsis {
    return <<"EOS"
EOS
}

sub help_detail {
    return <<"EOS"
Need documenation here.
EOS
}

sub execute {
    
    my $self = shift;


    my $biosql_namespace = $self->biosql_namespace();
   
    my $dbadp;

    if ($self->dev_flag()) {
    
        $dbadp = Bio::DB::BioDB->new(
                                     -database => 'biosql',
                                     -user     => 'sg_user',
                                     -pass     => 'sgus3r',
                                     -dbname   => 'DWDEV',
                                     -driver   => 'Oracle',
                                    );
    
    }
    else {
        
        $dbadp = Bio::DB::BioDB->new(
                                     -database => 'biosql',
                                     -user     => 'sg_user',
                                     -pass     => 'sg_us3r',
                                     -dbname   => 'DWRAC',
                                     -driver   => 'Oracle',
                                    );
        
    }
    
    my $feature_adp = $dbadp->get_object_adaptor('Bio::SeqFeatureI');

    my $feature_query = Bio::DB::Query::BioQuery->new();
    
    $feature_query->datacollections(['Bio::SeqFeatureI f']);
    $feature_query->where(['f.display_name = ?']);

    my $seq_adp = $dbadp->get_object_adaptor('Bio::SeqI');

    my $seq_query = Bio::DB::Query::BioQuery->new();

    $seq_query->datacollections(['Bio::SeqI s']);
    $seq_query->where(["s.display_id = ?"]);

    my $interpro_count = 0;

    my @features          = ( );
    my @interpro_features = ( );
    my %interpro_features = ( );
    
    foreach my $ref (@{$self->bio_seq_features()}) {

        my @fixup = ( );
        
        if (ref($ref) eq 'ARRAY') {
            @fixup = @{$ref};
        }
        else {
            push @fixup, $ref;
            
        }
        
        FIXUP: foreach my $feature (@fixup) {

            ## Sanity Check!
            ## Due to some workflow funkiness related to 
            ## the conjoined cat seq feature operations
            ## we might see some undefined values...
            ## ...which we don't want to try and call
            ## methods on and such.  
            unless (defined($feature)) { next FIXUP; }
            
            my $source_tag = $feature->source_tag();

            ## Sort into two piles so as to be able to
            ## optimize the database activity.
            if (defined($source_tag) && ($source_tag eq 'InterPro')) {
                push @interpro_features, $feature;
            }
            else {
                push @features, $feature;
            }
            
        }

    }
    
  INTERPRO: foreach my $feature (@interpro_features) {
        
        my $display_name = $feature->display_name();
        
        my ($seq_id, $source, $number);
       
        ##FIXME: This is better than the previous lame, pathetic fragile
        ##       hack that was here, but this module really should not be
        ##       parsing stuff out of the feature display name.  However,
        ##       it's the only way (at the moment) of batching the darn
        ##       InterPro features together to get good performance.  
        ##       The proper way to deal with this is probably to refactor
        ##       the pipeline to take a list of genes as input instead of a
        ##       protein fasta file.  The darn things have to be in BioSQL 
        ##       anway.  That way, the seq_id could be passed along with the
        ##       feature.  
        if (
            ($display_name =~ /^(\w+\d+)\.(\w+)\.(\d+)$/) ||
            ($display_name =~ /^(\w+\d+\.\d+)\.(\w+)\.(\d+)$/)
           ) {
            ($seq_id, $source, $number) = ($1, $2, $3);
        }
        else {
            die "failed to parse '$display_name':\n - does not match expected format (seqid.predictor.sequence)";
        }
        
        push @{$interpro_features{$seq_id}}, $feature;
        
    }
    
  SEQ: foreach my $seq_id (keys %interpro_features) {
        
        my $interpro_count = 0;
        
        my $result = $seq_adp->find_by_query(
                                             $seq_query,
                                             -name   => 'pap_upload_result_sequence',
                                             -values => [ $seq_id ],
                                         );
        
        my $db_seq = $result->next_object();
        
        unless (defined($db_seq)) {
            die "failed to find sequence object for '$seq_id'";
        }
        
        foreach my $feature (@{$interpro_features{$seq_id}}) {
            
            my $display_name = $feature->display_name();
            
            $interpro_count++;
            
            ## Change the name so we don't fetch it back below...
            $feature->display_name(join('.', $display_name, 'InterPro', $interpro_count));
            
            $db_seq->add_SeqFeature($feature);
            
        }
        
        $db_seq->store();
        
    }
    
  FEATURE: foreach my $feature (@features) {
        
        my $display_name = $feature->display_name();
        
        my $result = $feature_adp->find_by_query(
                                                 $feature_query,
                                                 -name   => 'pap_upload_result_feature',
                                                 -values => [ $display_name ],
                                             );
        
        my $db_feature = $result->next_object();
        
        unless (defined($db_feature)) {
            warn "failed to find feature object for '$display_name'";
            next FEATURE;
        }
        
        my $feature_ac    = $feature->annotation();
        my $db_feature_ac = $db_feature->annotation();
        
        foreach my $annotation ($feature_ac->get_Annotations()) {
            
            $db_feature_ac->add_Annotation($annotation);
            
        }
        
        foreach my $tagname (
                             qw(
                                psort_localization 
                                psort_score
                                kegg_evalue
                                kegg_description
                                blastp_bit_score
                                blastp_evalue
                                blastp_percent_identical
                                blastp_query_start
                                blastp_query_end
                                blastp_subject_start
                                blastp_subject_end
                                blastp_hit_name
                                blastp_hit_description
                                blastp_category
                               )
                         ) {
            
            if ($feature->has_tag($tagname)) {
                
                my ($tagvalue) = $feature->each_tag_value($tagname);
                $db_feature->add_tag_value($tagname, $tagvalue);
                
            }
            
        }
        
        $db_feature->store();
        
    }
    
    $feature_adp->commit();
    $seq_adp->commit();
    
    return 1;
    
}


1;
