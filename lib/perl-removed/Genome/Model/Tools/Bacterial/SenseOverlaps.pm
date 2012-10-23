package Genome::Model::Tools::Bacterial::SenseOverlaps;
use strict;
use warnings;
use Genome;
use Command;
use Carp;

# bap mess
use BAP::DB::DBI;
use BAP::DB::Tag;
use BAP::DB::GeneTag;

UR::Object::Type->define(
    class_name => __PACKAGE__,
    is => 'Command',
    has => [
        sequence_set_id => {
            is => 'Integer',
            doc => "sequence set id",
        },

    ],
    has_optional => [
        tag_name => {
            is => 'String',
            doc => "",
        },
        dev => {
            is => 'Boolean',
            doc => "switch to development mgap database",
        },
    ],
);

sub help_brief
{
    "tags genes that have 100 percent overlap";
}

sub help_detail
{
    return <<EOS
tags genes that have 100 percent overlap
EOS
}

sub execute
{
    my $self = shift;

    my $ssid = $self->sequence_set_id;
    
    # connect to mgap
    # one day, we may have to UR-ify.
    if($self->dev) {
        $BAP::DB::DBI::db_env = 'dev';
    }
    my $dbh = BAP::DB::DBI::db_Main();


    # prep query
    my $sth = $dbh->prepare($self->query);
    # execute query
    $sth->execute($ssid);

    # sort thru results.
    my $result_count = 0;
    my @genes2tag;
    my @genes4manual_review;
    while(my $results = $sth->fetchrow_arrayref) {
        # check PERCENTAGE_A and PERCENTAGE_B == 100
        my @columns = @{$results};
        # gene_name, seq_start, seq_end, strand(1/-1), other_gene_name, other_gene_start, other_gene_end, other_gene_strand, overlap,
        # pct_overlap, other_pct_overlap
        $result_count += 1;
        $self->status_message(" overlap ". $columns[0]."/".$columns[4]." pct ".
                             sprintf("%.02f : %.02f",$columns[9],$columns[10]));
        if(($columns[10] == 100.0) && ( $columns[3] == $columns[7])) {
            #print $columns[0],":",$columns[9]," ",$columns[4],":",$columns[10],"\n";
            push(@genes2tag,$columns[4]);
        }
        elsif(($columns[9] == 100.0) && ( $columns[3] == $columns[7])) {
            #print $columns[0],":",$columns[9]," ",$columns[4],":",$columns[10],"\n";
            push(@genes2tag,$columns[0]);
        }
        elsif($columns[10] == 100.0)  
        {
            push(@genes4manual_review,$columns[4]);
        }
        elsif($columns[9] == 100.0)  
        {
            push(@genes4manual_review,$columns[0]);
        }
    }
    $self->status_message("there are ". $result_count ." overlaps found");
    $self->status_message("there are " . scalar(@genes2tag) . " genes (matching strand) found 100% contained in a larger gene.");
    $self->status_message("there are " . scalar(@genes4manual_review) . " genes (diff strand) found 100% contained in a larger gene.");

    # tag 'em and bag 'em.
#    $self->tag_genes(\@genes2tag);
#    $self->tag_manual_review(\@genes4manual_review);

    return 1;
}


sub query {
    my $self = shift;

    my $query = " SELECT a.gene_name gene_name,
       a.seq_start seq_start,
       a.seq_end seq_end,
       b.strand strand,
       b.gene_name other_gene_name,
       b.seq_start other_gene_start,
       b.seq_end other_gene_end,
       b.strand other_gene_strand,
       least(a.seq_end, b.seq_end) - greatest(a.seq_start, b.seq_start) + 1,
       ((least(a.seq_end, b.seq_end) - greatest(a.seq_start, b.seq_start) + 1) / (b.seq_end - b.seq_start + 1)) * 100,
       ((least(a.seq_end, b.seq_end) - greatest(a.seq_start, b.seq_start) + 1) / (a.seq_end - a.seq_start + 1)) * 100
  FROM coding_gene a,
       coding_gene b,
       dna_sequence
 WHERE a.sequence_id = dna_sequence.sequence_id AND
       dna_sequence.sequence_set_id = ? AND
       a.sequence_id = b.sequence_id AND
       a.gene_id != b.gene_id AND
       a.gene_id < b.gene_id AND
       a.pfam_evidence = 1 AND
       b.pfam_evidence = 1 AND
       a.phase_5 = 1 AND
       b.phase_5 = 1 AND
       a.seq_start != b.seq_start AND
       a.seq_end != b.seq_end AND
       a.seq_start != b.seq_end AND
       a.seq_end != b.seq_start AND
       (
        (a.seq_start between b.seq_start and b.seq_end) OR
        (a.seq_end between b.seq_start and b.seq_end) OR
        (b.seq_start between a.seq_start and a.seq_end) OR
        (b.seq_end between a.seq_start and a.seq_end)
       ) ";
    return $query;
}



1;
