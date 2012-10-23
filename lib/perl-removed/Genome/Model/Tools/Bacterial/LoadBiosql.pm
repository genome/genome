package Genome::Model::Tools::Bacterial::LoadBiosql;

use strict;
use warnings;

use Genome;
use Command;

use BAP::DB::SequenceSet;
use BAP::DB::GeneTag;
use BAP::DB::Tag;

use Bio::Annotation::SimpleValue;
use Bio::DB::BioDB;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::SeqIO;
use Bio::Species;

use Data::Dumper;
use Getopt::Long;
use Carp;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is         => 'Command',
    has        => [
        'sequence_set_id' => {
            is  => 'Integer',
            doc => "sequence set id in mgap",
        },
    ],
    has_optional => [
        'biosql_dev' => {
            is      => 'Boolean',
            doc     => "use development biosql db",
            default => 0,
        },

        'dev' => {
            is      => 'Boolean',
            doc     => "use development mgap db",
            default => 0,
        },
    ],
);

sub help_brief
{
    "tool for loading mgap data to biosql",;
}

sub help_detail
{
    return <<EOS
loads mgap data to biosql
EOS
}

sub help_synopsis
{
    return <<EOS
gmt bacterial load-biosql --sequence-set-id ssid  [--dev] [--biosql-dev]
EOS
}

sub execute
{
    my $self = shift;

    if ( $self->dev ) { $BAP::DB::DBI::db_env = 'dev'; }

    my $biosql_db;

    if ($self->biosql_dev)
    {

        $biosql_db = Bio::DB::BioDB->new(
            -database => 'biosql',
            -user     => 'sg_user',
            -pass     => 'sgus3r',
            -dbname   => 'DWDEV',
            -driver   => 'Oracle',
        );

    }
    else
    {

        $biosql_db = Bio::DB::BioDB->new(
            -database => 'biosql',
            -user     => 'sg_user',
            -pass     => 'sg_us3r',
            -dbname   => 'DWRAC',
            -driver   => 'Oracle',
        );
    }

    my $biosql_oa = $biosql_db->get_object_adaptor('Bio::Seq');

    $biosql_oa->dbh->{'LongTruncOk'} = 0;
    $biosql_oa->dbh->{'LongReadLen'} = 1000000000;

    my $sequence_set   = BAP::DB::SequenceSet->retrieve($self->sequence_set_id);
    my @sequences      = $sequence_set->sequences();
    my @sequence_names = ();

    my $mgap_coding_gene_count = 0;
    my $mgap_trna_gene_count   = 0;
    my $mgap_rna_gene_count    = 0;

SEQUENCE: foreach my $sequence (@sequences)
    {

        my $sequence_id   = $sequence->sequence_id();
        my $sequence_name = $sequence->sequence_name();

        push @sequence_names, $sequence_name;

        my $bp_seq = Bio::Seq->new(
            -id               => $sequence_name,
            -accession_number => $sequence_id,
            -seq              => $sequence->sequence_string(),
        );

        my $sequence_length = $bp_seq->length();

        my @coding_genes = $sequence->coding_genes( 'phase_5' => 1 );
        my @trna_genes   = $sequence->trna_genes();
        my @rna_genes    = $sequence->rna_genes();

        @rna_genes = grep { !$_->redundant() } @rna_genes;

        $mgap_coding_gene_count += @coding_genes;
        $mgap_trna_gene_count   += @trna_genes;
        $mgap_rna_gene_count    += @rna_genes;

        my @genes = ( @coding_genes, @trna_genes, @rna_genes );

        unless (@genes)
        {
            warn
                "no genes found in MGAP for sequence '$sequence_name' ($sequence_length bp)";
            next SEQUENCE;
        }

        foreach my $gene (@genes)
        {

            my $gene_name = $gene->gene_name();
            my $start     = $gene->start();
            my $end       = $gene->end();
            my $strand    = $gene->strand();
            my $source    = $gene->source();
            my $score     = $gene->score();

            if ( $start > $end ) { ( $start, $end ) = ( $end, $start ); }

            my $primary = 'gene';

            my $class = ref($gene);

            if ( $class eq 'BAP::DB::tRNAGene' )
            {
                $primary = 'trna';
            }
            elsif ( $class eq 'BAP::DB::RNAGene' )
            {
                $primary = 'rrna';
            }

            my $feature = Bio::SeqFeature::Generic->new(
                -seq_id       => $sequence_name,
                -display_name => $gene_name,
                -primary      => $primary,
                -source       => $source,
                -start        => $start,
                -end          => $end,
                -strand       => $strand,
                -score        => $score,
            );

            if ( $class eq 'BAP::DB::CodingGene' )
            {

                if ( $gene->missing_start() )
                {
                    $feature->add_tag_value( 'missing_start' => 1 );
                }

                if ( $gene->missing_stop() )
                {
                    $feature->add_tag_value( 'missing_stop' => 1 );
                }

                my @genetags
                    = BAP::DB::GeneTag->search( { gene_id => $gene } );
                foreach my $tag (@genetags)
                {
                    $feature->add_tag_value(
                        $tag->tag_id->tag_name => $tag->tag_id->tag_value );
                }

            }
            elsif ( $class eq 'BAP::DB::tRNAGene' )
            {
                $feature->add_tag_value( 'amino_acid', $gene->aa() );
            }
            elsif ( $class eq 'BAP::DB::RNAGene' )
            {
                if ( $gene->source() eq 'rnammer' )
                {
                    $feature->add_tag_value( 'description' => $gene->desc() );
                }
                elsif ( $gene->source() eq 'rfam' )
                {
                    $feature->add_tag_value( 'description' => $gene->desc() );
                    $feature->add_tag_value( 'rfam_acc'    => $gene->acc() );
                    $feature->add_tag_value(
                        'rfam_product' => $gene->rfam_prod() );

                   # added due to bioperl-db not serializine a feature's score
                    $feature->add_tag_value( 'rfam_score' => $score );
                }
            }

            carp "adding seq feature for $gene_name";
            print STDERR Dumper( \$feature ), "\n";
            $bp_seq->add_SeqFeature($feature);
            carp "added seq feature for $gene_name";
        }
        carp "done adding features.";
        $bp_seq->namespace('MGAP');

        carp "done setting namespace";
        my $lseq = $biosql_oa->find_by_unique_key($bp_seq);

        carp "done checking unique key";
        if ( defined($lseq) )
        {
            die
                "sequence '$sequence_name' already exists in BioSQL, sorry, you'll have to delete it manually";
        }

        my $pseq = $biosql_db->create_persistent($bp_seq);
        carp "creating persistent";
        $pseq->store();
        carp "storing...";

    }

    my $locus_tag = $self->locus_tag(@sequence_names);

    unless ( defined($locus_tag)
        && ( length($locus_tag) > 3 ) )
    {

        die
            "locus_tag '$locus_tag' calculated from sequence names seems wrong";

    }

    $biosql_oa->commit();

    my $biosql_coding_gene_count = 0;
    my $biosql_rna_gene_count    = 0;
    my $biosql_trna_gene_count   = 0;

    my $query = Bio::DB::Query::BioQuery->new();

    $query->datacollections( [ "Bio::PrimarySeqI s", ] );
    $query->where( ["s.display_id like '$locus_tag%'"] );

    my $res = $biosql_oa->find_by_query($query);

    while ( my $seq = $res->next_object() )
    {

        my @features = $seq->get_SeqFeatures();

        foreach my $feature (@features)
        {

            my $primary = $feature->primary_tag();

            if ( $primary eq 'gene' )
            {
                $biosql_coding_gene_count++;
            }
            elsif ( $primary eq 'rrna' )
            {
                $biosql_rna_gene_count++;
            }
            elsif ( $primary eq 'trna' )
            {
                $biosql_trna_gene_count++;
            }

        }

    }

    unless ( ( $mgap_coding_gene_count == $biosql_coding_gene_count )
        && ( $mgap_rna_gene_count == $biosql_rna_gene_count )
        && ( $mgap_trna_gene_count == $biosql_trna_gene_count ) )
    {

        my $expected = join '/',
            $mgap_coding_gene_count,
            $mgap_rna_gene_count,
            $mgap_trna_gene_count;

        my $found = join '/',
            $biosql_coding_gene_count,
            $biosql_rna_gene_count,
            $biosql_trna_gene_count;

        die
            "gene count mismatch expected $expected coding/rna/trna, found $found";

    }

    print "loaded gene count: ",
        join( "\t",
        $biosql_coding_gene_count, $biosql_rna_gene_count,
        $biosql_trna_gene_count ),
        "\n";

    return 1;
}

sub locus_tag
{
    my $self = shift;
    return '' unless @_;

    my $prefix = shift;

    for (@_)
    {
        chop $prefix while ( !/^\Q$prefix\E/ );
    }

    return $prefix;

}

1;
