use strict;
use warnings;

#use lib '/gscuser/mjohnson/bioperl-svn/bioperl-live';
#use lib '/gscuser/mjohnson/bioperl-svn/bioperl-run';

use above 'PAP';
use Workflow;

use Bio::DB::BioDB;
use Bio::DB::Query::BioQuery;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;

use Data::Dumper;
use File::Temp;
use Test::More tests => 4;

BEGIN {
    use_ok('PAP::Command');
    use_ok('PAP::Command::UploadResult');
}

my $biosql_namespace = 'AUTOMATED_TEST';

my $db_adp = connect_db();

my $seq_adp = $db_adp->get_object_adaptor('Bio::Seq');

my $bp_seq = create_seq($biosql_namespace);

my $pseq = $db_adp->create_persistent($bp_seq);

$pseq->store();
$seq_adp->commit();

my $feature_ref = create_feature_ref();

my $command = PAP::Command::UploadResult->create();

$command->dev_flag(1);
$command->biosql_namespace($biosql_namespace);
$command->bio_seq_features($feature_ref);

isa_ok($command, 'PAP::Command::UploadResult');

ok($command->execute());


$seq_adp->rollback();

$pseq->remove();
$seq_adp->commit();

sub connect_db {
    
    return Bio::DB::BioDB->new(
                               -database => 'biosql',
                               -user     => 'sg_user',
                               -pass     => 'sgus3r',
                               -dbname   => 'DWDEV',
                               -driver   => 'Oracle',
                           );
    
}

sub create_seq {

    my ($namespace) = @_;

    
    my $seq = Bio::Seq->new(
                            -id               => 'TST_Contig98.1',
                            -accession_number => '999999999', 
                            -namespace        => $namespace,
                            -seq              => 'GATTACA' x 1000, 
                           );
    
    $seq->add_SeqFeature(
                         Bio::SeqFeature::Generic->new(
                                                       -seq_id       => 'TST_Contig98.1',
                                                       -display_name => 'TST_Contig98.1.GeneMark.1',
                                                       -primary      => 'gene',
                                                       -source       => 'genemark',
                                                       -start        => 1,
                                                       -end          => 1000,
                                                       -strand       => 1,
                                                      )
                        );

    return $seq;
    
}

sub create_feature_ref {

    my $display_name = 'TST_Contig98.1.GeneMark.1';

    my $psort = Bio::SeqFeature::Generic->new(
                                              -display_name => $display_name,
                                              -tag          => {
                                                                'psort_localization' => 'Cytoplasmic',
                                                                'psort_score'        => 8.96,
                                                               },
                                             );


    my $hit_description = 'transposase, putative [Deinococcus radiodurans] R1';                

    my $blastp = Bio::SeqFeature::Generic->new(
                                                -display_name => $display_name,
                                                -tag          => {
                                                                  'blastp_bit_score'         => 260,
                                                                  'blastp_evalue'            => 2.1e-21,
                                                                  'blastp_percent_identical' => 29.4,
                                                                  'blastp_query_start'       => 15,
                                                                  'blastp_query_end'         => 241,
                                                                  'blastp_subject_start'     => 13,
                                                                  'blastp_subject_end'       => 246,
                                                                  'blastp_hit_name'          => 'ref|NP_051595.1|',
                                                                  'blastp_hit_description'   => $hit_description,
                                                                  'blastp_category'          => 'Predicted Protein',
                                                                 },
                                              );

    $blastp->annotation->add_Annotation(
                                        'dblink',
                                         Bio::Annotation::DBLink->new(
                                                                      -database   => 'GenBank',
                                                                      -primary_id => 'NP_051595.1',
                                                                     ),
                                       );


    my $kegg = Bio::SeqFeature::Generic->new(
                                             -display_name     => $display_name,
                                             -kegg_evalue      => 4.2e-197,
                                             -kegg_description => 'hypothetical protein',
                                            );

    $kegg->annotation->add_Annotation(
                                      'dblink',
                                      Bio::Annotation::DBLink->new(
                                                                   -database   => 'KEGG',
                                                                   -primary_id => 'eci:UTI89_C4937',
                                                                  ),
                                     );

    $kegg->annotation->add_Annotation(
                                      'dblink',
                                      Bio::Annotation::DBLink->new(
                                                                   -database   => 'KEGG',
                                                                   -primary_id => 'K01152',
                                                                  ),
                                     );


    my $interpro = Bio::SeqFeature::Generic->new(
                                                 -display_name => $display_name,
                                                 -primary      => 'HMMPfam',
                                                 -source_tag   => 'InterPro',
                                                 -start        => 1,
                                                 -end          => 100,
                                                 -score        => 2.3e-12,
                                                 -tag          => {
                                                                   'interpro_analysis'    => 'HMMPfam',
                                                                   'interpro_evalue'      => 2.3e-12,
                                                                   'interpro_description' => 'ATP-binding region, ATPase-like',
                                                                  },
                                                );

    $interpro->annotation->add_Annotation(
                                          'dblink',
                                          Bio::Annotation::DBLink->new(
                                                                       -database   => 'InterPro',
                                                                       -primary_id => 'IPR003594',
                                                                      ),
                                         );
    
    return [ [ $psort, $blastp, $kegg, $interpro, ] ];

}
