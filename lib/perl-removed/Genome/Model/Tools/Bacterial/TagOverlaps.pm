package Genome::Model::Tools::Bacterial::TagOverlaps;

# bdericks: Holy shit. If you needed proof of the need to refactor BAP to 
# use UR, this is it. Take a look at that query.

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

sub execute {
    my $self = shift;

    my $ssid = $self->sequence_set_id;
    
    # Connect to MGAP schema
    if($self->dev) {
        $BAP::DB::DBI::db_env = 'dev';
    }
    my $dbh = BAP::DB::DBI::db_Main();

    # Prepare and execute the query
    my $sth = $dbh->prepare($self->query);
    $sth->execute($ssid);

    # Get the tag that will be applied to these overlaps
    my ($manual_review_tag) = BAP::DB::Tag->search({
        tag_name => 'ManualReview',
        tag_value => 'fully overlapped gene',
    });

    my ($dead_tag) = BAP::DB::Tag->search({
        tag_name => 'Dead',
        tag_value => 'fully overlapped gene',
    });

    while(my $results = $sth->fetchrow_hashref) {
        next unless $results->{pct_overlap} == 100.0 or $results->{other_pct_overlap} == 100.0;
        my $gene_length = abs($results->{seq_end} - $results->{seq_start}) + 1;
        my $other_gene_length = abs($results->{other_gene_start} - $results->{other_gene_end}) + 1;

        my $dead_gene;
        my $manual_review_gene;
        if ($gene_length < $other_gene_length) {
            $dead_gene = $results->{gene_name};
            $manual_review_gene = $results->{other_gene_name};
        }
        else {
            $dead_gene = $results->{other_gene_name};
            $manual_review_gene = $results->{gene_name};
        }

        $self->debug_message("$dead_gene will be tagged as dead.");
        $self->debug_message("$manual_review_gene will be tagged for manual review.");

        # Tagging dead gene
        my $dead_coding_gene_results = BAP::DB::CodingGene->search({gene_name => $dead_gene});
        my $dead_coding_gene = $dead_coding_gene_results->next;
        next unless defined $dead_coding_gene;

        my $dead_gene_tag = BAP::DB::GeneTag->find_or_create({
            gene_id => $dead_coding_gene,
            tag_id => $dead_tag,
        });

        # Tagging manual review gene
        my $review_coding_gene_results = BAP::DB::CodingGene->search({gene_name => $manual_review_gene});
        my $review_coding_gene = $review_coding_gene_results->next;
        next unless defined $review_coding_gene;

        my $review_gene_tag = BAP::DB::GeneTag->find_or_create({
            gene_id => $review_coding_gene,
            tag_id => $manual_review_tag,
        });
    }    

    # Check if any no commit flags are set...
    if(exists($ENV{UR_DBI_NO_COMMIT}) && ($ENV{UR_DBI_NO_COMMIT} == 1)) {
        $self->debug_message("UR_DBI_NO_COMMIT set; not commiting changes");
    }
    else {
        $self->debug_message("Committing changes to MGAP database.");
        BAP::DB::DBI->dbi_commit();
    }

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
       least(a.seq_end, b.seq_end) - greatest(a.seq_start, b.seq_start) + 1 overlap,
       ((least(a.seq_end, b.seq_end) - greatest(a.seq_start, b.seq_start) + 1) / (b.seq_end - b.seq_start + 1)) * 100 pct_overlap,
       ((least(a.seq_end, b.seq_end) - greatest(a.seq_start, b.seq_start) + 1) / (a.seq_end - a.seq_start + 1)) * 100 other_pct_overlap
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
